#!/usr/bin/env python3
"""
Alternative Electron Carriers Experiment.

Explores how different quinone redox potentials affect proton pumping
capacity and ATP yield in the PTR-94 Perfect Coupling Module. Tests
ubiquinone (native), plastoquinone, menaquinone, and engineered
variants with custom redox potentials.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple


from ..thermodynamics import (
    STANDARD_GIBBS_GLUCOSE,
    STANDARD_GIBBS_ATP,
    STANDARD_REDOX_POTENTIALS,
    max_work_from_redox,
    compute_coupling_efficiency,
)
from ..pareto_optimizer import (
    DesignParameters,
    ObjectiveFunction,
)

# ---------------------------------------------------------------------------
# Carrier definitions
# ---------------------------------------------------------------------------

CARRIER_DATABASE: Dict[str, Dict[str, Any]] = {
    "ubiquinone": {
        "e0_prime": 0.040,
        "n_electrons": 2,
        "pmf_yield_per_e": 1.0,
        "source": "Native mitochondrial",
        "description": "Ubiquinone/ubiquinol (CoQ) — Standard mitochondrial carrier",
    },
    "plastoquinone": {
        "e0_prime": 0.110,
        "n_electrons": 2,
        "pmf_yield_per_e": 0.85,
        "source": "Chloroplast",
        "description": "Plastoquinone — Higher potential, lower PMF yield",
    },
    "menaquinone": {
        "e0_prime": -0.074,
        "n_electrons": 2,
        "pmf_yield_per_e": 1.15,
        "source": "Bacterial (e.g., E. coli)",
        "description": "Menaquinone — Lower potential, higher PMF yield",
    },
    "engineered_low": {
        "e0_prime": -0.150,
        "n_electrons": 2,
        "pmf_yield_per_e": 1.30,
        "source": "Hypothetical engineered",
        "description": "Engineered low-potential quinone — Maximises PMF",
    },
    "engineered_high": {
        "e0_prime": 0.200,
        "n_electrons": 2,
        "pmf_yield_per_e": 0.70,
        "source": "Hypothetical engineered",
        "description": "Engineered high-potential quinone — Minimises PMF",
    },
}
"""
ASSUMPTION: Standard reduction potentials at pH 7 (E0') from Nicholls &
Ferguson (2013) and Berg et al. (2015). The PMF yield per electron is a
relative estimate of how much proton pumping capacity the quinone provides
compared to ubiquinone (baseline = 1.0). Engineered quinones are speculative.
"""


@dataclass
class CarrierResult:
    """Results for a single electron carrier evaluation.

    Attributes:
        carrier_name: Name of the carrier.
        e0_prime: Standard reduction potential (V).
        atp_yield: ATP per glucose with this carrier.
        efficiency: Thermodynamic efficiency (percent).
        redox_work: Maximum redox work from NADH to O2 (kJ/mol).
        pmf_per_electron: Relative PMF yield per electron.
        proton_pumping_capacity: Estimated H+ pumped per NADH.
        description: Carrier description.
    """
    carrier_name: str
    e0_prime: float
    atp_yield: float
    efficiency: float
    redox_work: float
    pmf_per_electron: float
    proton_pumping_capacity: float
    description: str


@dataclass
class AlternativeCarrierResult:
    """Container for alternative carrier experiment results.

    Attributes:
        carriers: List of CarrierResult for each tested carrier.
        standard_carrier: Name of the reference (ubiquinone).
        best_carrier: Name of the best-performing carrier.
        improvement: ATP yield improvement of best over standard.
        pmf_range: (min, max) PMF yield across carriers.
        summary: Formatted summary text.
    """
    carriers: List[CarrierResult] = field(default_factory=list)
    standard_carrier: str = "ubiquinone"
    best_carrier: str = ""
    improvement: float = 0.0
    pmf_range: Tuple[float, float] = (0.0, 0.0)
    summary: str = ""


class AlternativeCarrierExperiment:
    """Experiment testing alternative electron carriers.

    Evaluates how the quinone pool composition affects the PCM's proton
    pumping capacity, and consequently the ATP yield.

    Attributes:
        carriers_to_test: List of carrier names from CARRIER_DATABASE.
            If None, tests all.
        base_params: Nominal DesignParameters for the PCM.
        result: AlternativeCarrierResult after run().
    """

    def __init__(
        self,
        carriers_to_test: List[str] | None = None,
        base_params: DesignParameters | None = None,
    ) -> None:
        """
        Parameters
        ----------
        carriers_to_test : list of str, optional
            Names of carriers to test. Defaults to all in CARRIER_DATABASE.
        base_params : DesignParameters, optional
            Baseline PCM parameters.
        """
        if carriers_to_test is None:
            self.carriers_to_test = list(CARRIER_DATABASE.keys())
        else:
            self.carriers_to_test = [c for c in carriers_to_test
                                     if c in CARRIER_DATABASE]
        self.base_params = base_params or DesignParameters()
        self.result: AlternativeCarrierResult | None = None

        # ASSUMPTION: The PCM's proton pumping stoichiometry scales linearly
        # with the redox gap between NADH and the quinone. A wider gap allows
        # more H+ to be pumped per NADH oxidised.

    def run(self) -> AlternativeCarrierResult:
        """Evaluate all selected carriers and return a comparison.

        Returns
        -------
        AlternativeCarrierResult with per-carrier breakdown.
        """
        # Reference: NADH -> O2 redox work
        nadh_potential = STANDARD_REDOX_POTENTIALS["NAD_NADH"][0]
        o2_potential = STANDARD_REDOX_POTENTIALS["O2_H2O"][0]

        results: List[CarrierResult] = []
        for name in self.carriers_to_test:
            info = CARRIER_DATABASE[name]
            e0 = info["e0_prime"]
            pmf_per_e = info["pmf_yield_per_e"]

            # Redox work for NADH -> carrier -> O2
            # ASSUMPTION: The carrier sits between Complex I and III/IV.
            # The total work is the sum of two steps.
            work_step1 = max_work_from_redox(nadh_potential, e0, 2)
            work_step2 = max_work_from_redox(e0, o2_potential, 2)
            total_work = work_step1 + work_step2

            # Estimate proton pumping capacity relative to ubiquinone
            # ASSUMPTION: Baseline ubiquinone gives 30 H+/NADH in PCM.
            base_carrier = CARRIER_DATABASE.get("ubiquinone", {})
            base_pmf = base_carrier.get("pmf_yield_per_e", 1.0)
            relative_pmf = pmf_per_e / base_pmf if base_pmf > 0 else 1.0
            proton_capacity = relative_pmf * self.base_params.h_per_nadh

            # Estimate ATP yield
            # ASSUMPTION: The yield scales directly with the proton pumping
            # capacity of the carrier. Other parameters remain fixed.
            scaled_params = DesignParameters(
                h_per_nadh=min(proton_capacity, 40.0),
                h_per_fadh2=float(
                    self.base_params.h_per_fadh2 * relative_pmf
                ),
                h_per_atp=self.base_params.h_per_atp,
                atp_synthase_efficiency=self.base_params.atp_synthase_efficiency,
                membrane_leak_conductance=self.base_params.membrane_leak_conductance,
                proton_slip_probability=self.base_params.proton_slip_probability,
                roq_quinone_coupling=self.base_params.roq_quinone_coupling,
                scaffold_channeling_efficiency=self.base_params.scaffold_channeling_efficiency,
                ros_bypass_fraction=self.base_params.ros_bypass_fraction,
            )
            atp = ObjectiveFunction.compute_atp_yield(scaled_params)
            eff = compute_coupling_efficiency(
                atp, STANDARD_GIBBS_GLUCOSE, STANDARD_GIBBS_ATP
            )

            results.append(CarrierResult(
                carrier_name=name,
                e0_prime=e0,
                atp_yield=atp,
                efficiency=eff,
                redox_work=total_work,
                pmf_per_electron=pmf_per_e,
                proton_pumping_capacity=proton_capacity,
                description=info.get("description", ""),
            ))

        # Identify best carrier (highest ATP yield)
        best = max(results, key=lambda r: r.atp_yield)
        standard = next(
            (r for r in results if r.carrier_name == self.carriers_to_test[0]),
            results[0],
        )
        improvement = best.atp_yield - standard.atp_yield
        pmf_vals = [r.pmf_per_electron for r in results]
        pmf_range = (min(pmf_vals), max(pmf_vals))

        self.result = AlternativeCarrierResult(
            carriers=results,
            standard_carrier=results[0].carrier_name,
            best_carrier=best.carrier_name,
            improvement=improvement,
            pmf_range=pmf_range,
        )
        self.result.summary = self._generate_summary()
        return self.result

    def _generate_summary(self) -> str:
        """Generate a formatted text summary of results."""
        if self.result is None or not self.result.carriers:
            return ""
        r = self.result

        e0_prime_header = "E0' (V)"
        lines = [
            "=" * 84,
            "  ALTERNATIVE ELECTRON CARRIERS EXPERIMENT — PTR-94",
            "=" * 84,
            f"  {'Carrier':<20} {e0_prime_header:>10} {'ATP':>8} {'Eff%':>7} "
            f"{'Redox(kJ)':>10} {'PMF/e':>7} {'H+/NADH':>8}",
            "-" * 84,
        ]

        for carrier in r.carriers:
            lines.append(
                f"  {carrier.carrier_name:<20} {carrier.e0_prime:>+10.4f} "
                f"{carrier.atp_yield:>8.2f} {carrier.efficiency:>6.2f} "
                f"{carrier.redox_work:>10.1f} {carrier.pmf_per_electron:>7.2f} "
                f"{carrier.proton_pumping_capacity:>8.1f}"
            )

        lines.extend([
            "-" * 84,
            f"  Standard carrier:  {r.standard_carrier}",
            f"  Best carrier:      {r.best_carrier} "
            f"(ΔATP = {r.improvement:+.2f})",
            f"  PMF yield range:   [{r.pmf_range[0]:.2f}, {r.pmf_range[1]:.2f}]",
            "",
            "  Notes:",
            "  * PMF/electron is relative to ubiquinone (baseline = 1.0).",
            "  * Engineered carriers are speculative and assume perfect tuning.",
            "  * Real-world yields may differ due to kinetic constraints.",
            "=" * 84,
        ])
        return "\n".join(lines)

    def data_table(self) -> List[Dict[str, Any]]:
        """Return structured data for export."""
        if self.result is None:
            return []
        return [
            {
                "carrier": c.carrier_name,
                "e0_prime": c.e0_prime,
                "atp_yield": c.atp_yield,
                "efficiency_pct": c.efficiency,
                "redox_work_kj": c.redox_work,
                "pmf_per_electron": c.pmf_per_electron,
                "proton_pumping_capacity": c.proton_pumping_capacity,
                "description": c.description,
            }
            for c in self.result.carriers
        ]


if __name__ == "__main__":
    exp = AlternativeCarrierExperiment()
    result = exp.run()
    print(result.summary)

    print("\n  Data table:")
    for row in exp.data_table():
        print(f"  {row['carrier']:<20} ATP={row['atp_yield']:.2f}  "
              f"E0'={row['e0_prime']:+.4f}V  "
              f"PMF/e={row['pmf_per_electron']:.2f}")
