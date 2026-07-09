#!/usr/bin/env python3
"""
Enzyme Efficiency Experiment.

Explores how enzyme catalytic efficiency (kcat/Km) affects pathway flux
and ATP yield in the PTR-94 network. Scans kcat/Km for each module enzyme
to identify kinetic bottlenecks and determine which reactions most limit
the overall ATP production rate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple

import numpy as np

from ..kinetics import (
    EnzymeKinetics,
)
from ..pareto_optimizer import (
    DesignParameters,
    ObjectiveFunction,
)


# ---------------------------------------------------------------------------
# Enzyme module definitions
# ---------------------------------------------------------------------------

ENZYME_MODULES: Dict[str, Dict[str, Any]] = {
    "hexokinase": {
        "reaction": "Glucose + ATP -> G6P + ADP",
        "kcat_km_nominal": 1e6,
        "substrate": "Glucose",
        "pathway": "glycolysis",
        "km": 1e-4,
        "vmax": 1e-3,
    },
    "phosphofructokinase": {
        "reaction": "F6P + ATP -> F16BP + ADP",
        "kcat_km_nominal": 5e5,
        "substrate": "F6P",
        "pathway": "glycolysis",
        "km": 5e-5,
        "vmax": 1e-3,
    },
    "glyceraldehyde_3p_dehydrogenase": {
        "reaction": "G3P + NAD+ + Pi -> BPG + NADH + H+",
        "kcat_km_nominal": 2e5,
        "substrate": "G3P",
        "pathway": "glycolysis",
        "km": 4e-5,
        "vmax": 1e-3,
    },
    "pyruvate_kinase": {
        "reaction": "PEP + ADP -> Pyr + ATP",
        "kcat_km_nominal": 8e5,
        "substrate": "PEP",
        "pathway": "glycolysis",
        "km": 2e-5,
        "vmax": 1e-3,
    },
    "pyruvate_dehydrogenase": {
        "reaction": "Pyr + CoA + NAD+ -> AcCoA + NADH + CO2",
        "kcat_km_nominal": 3e5,
        "substrate": "Pyr",
        "pathway": "pdh_tca",
        "km": 5e-4,
        "vmax": 1e-3,
    },
    "citrate_synthase": {
        "reaction": "AcCoA + OAA + H2O -> Cit + CoA + H+",
        "kcat_km_nominal": 4e5,
        "substrate": "AcCoA",
        "pathway": "pdh_tca",
        "km": 1e-4,
        "vmax": 1e-3,
    },
    "isocitrate_dehydrogenase": {
        "reaction": "IsoCit + NAD+ -> AKG + NADH + CO2",
        "kcat_km_nominal": 2e5,
        "substrate": "IsoCit",
        "pathway": "pdh_tca",
        "km": 1e-4,
        "vmax": 1e-3,
    },
    "alpha_kg_dehydrogenase": {
        "reaction": "AKG + CoA + NAD+ -> SucCoA + NADH + CO2",
        "kcat_km_nominal": 1e5,
        "substrate": "AKG",
        "pathway": "pdh_tca",
        "km": 2e-4,
        "vmax": 1e-3,
    },
    "pcm_complex_I": {
        "reaction": "NADH + Q -> NAD + QH2 (H+ pumping)",
        "kcat_km_nominal": 5e5,
        "substrate": "NADH",
        "pathway": "pcm",
        "km": 1e-4,
        "vmax": 1e-3,
    },
    "pcm_complex_III": {
        "reaction": "QH2 + 2CytC_ox -> Q + 2CytC_red (H+ pumping)",
        "kcat_km_nominal": 3e5,
        "substrate": "QH2",
        "pathway": "pcm",
        "km": 1e-4,
        "vmax": 1e-3,
    },
    "pcm_complex_IV": {
        "reaction": "2CytC_red + 0.5O2 -> 2CytC_ox + H2O (H+ pumping)",
        "kcat_km_nominal": 4e5,
        "substrate": "CytC_red",
        "pathway": "pcm",
        "km": 1e-4,
        "vmax": 1e-3,
    },
    "pcm_atp_synthase": {
        "reaction": "ADP + Pi + 3H+ -> ATP + H2O",
        "kcat_km_nominal": 6e5,
        "substrate": "ADP",
        "pathway": "pcm",
        "km": 1e-5,
        "vmax": 1e-3,
    },
}
"""
ASSUMPTION: Nominal kcat/Km values (in M^-1 s^-1) are estimated from typical
enzyme kinetic data for the native enzymes (Berg et al., 2015). The PCM
complex enzymes are hypothetical but use values similar to their natural
counterparts. The ranking of sensitivities is more robust than absolute values.
"""


@dataclass
class EnzymeEfficiencyResult:
    """Container for enzyme efficiency experiment results.

    Attributes:
        enzyme_names: Names of tested enzymes.
        efficiency_multipliers: Tested multiplier values (e.g., [0.1, 1.0, 10.0]).
        atp_yields: (n_enzymes, n_multipliers) array of ATP yields.
        flux_values: (n_enzymes, n_multipliers) array of relative flux.
        base_yield: ATP yield at nominal efficiency.
        bottleneck_ranking: List of (enzyme_name, sensitivity_score) sorted.
        summary: Formatted text summary.
    """
    enzyme_names: List[str]
    efficiency_multipliers: np.ndarray
    atp_yields: np.ndarray
    flux_values: np.ndarray
    base_yield: float = 0.0
    bottleneck_ranking: List[Tuple[str, float]] = field(default_factory=list)
    summary: str = ""


class EnzymeEfficiencyExperiment:
    """Experiment scanning enzyme efficiency for each module enzyme.

    For each enzyme, the kcat/Km is varied from 0.1x to 10x the nominal
    value, and the effect on ATP yield is recorded. Enzymes with the
    largest impact are identified as kinetic bottlenecks.

    Attributes:
        multiplier_range: (low, high, n_steps) for efficiency scaling.
        enzymes_to_test: List of enzyme names. If None, tests all.
        base_params: Nominal DesignParameters.
        result: EnzymeEfficiencyResult after run().
    """

    def __init__(
        self,
        multiplier_range: Tuple[float, float, int] = (0.1, 10.0, 21),
        enzymes_to_test: List[str] | None = None,
        base_params: DesignParameters | None = None,
    ) -> None:
        """
        Parameters
        ----------
        multiplier_range : tuple
            (low, high, n_steps) for folding the nominal kcat/Km.
        enzymes_to_test : list of str, optional
            Enzyme names. Defaults to all in ENZYME_MODULES.
        base_params : DesignParameters, optional
            Baseline PCM parameters.
        """
        self.multiplier_range = multiplier_range
        if enzymes_to_test is None:
            self.enzymes_to_test = list(ENZYME_MODULES.keys())
        else:
            self.enzymes_to_test = enzymes_to_test
        self.base_params = base_params or DesignParameters()
        self.result: EnzymeEfficiencyResult | None = None

        # ASSUMPTION: kcat/Km scaling captures both catalytic rate and
        # substrate affinity changes. In reality, kcat and Km can change
        # independently, but for bottleneck analysis the lumped metric
        # is a useful first approximation.

    def _compute_yield_for_efficiency(
        self,
        enzyme_name: str,
        multiplier: float,
    ) -> float:
        """Compute ATP yield at a given efficiency multiplier for one enzyme.

        Parameters
        ----------
        enzyme_name : str
            Enzyme to perturb.
        multiplier : float
            Fold-change in kcat/Km (1.0 = nominal).

        Returns
        -------
        float
            Estimated ATP yield per glucose.
        """
        info = ENZYME_MODULES[enzyme_name]
        # ASSUMPTION: The ATP yield is proportional to the flux through
        # each enzyme. We use a simplified model where flux scales as:
        #   flux = Vmax * [S] / (Km + [S])
        # with Km = kcat / (kcat/Km) approx.
        # Since Vmax = kcat * [E], scaling kcat/Km by multiplier while
        # holding Vmax constant means Km = Km_nominal / multiplier.
        km_nominal = info["km"]
        vmax = info["vmax"]
        substrate_conc = 1e-4  # M, typical

        # Effective Km after scaling
        km_eff = km_nominal / max(multiplier, 0.01)
        rate = EnzymeKinetics.michaelis_menten(
            substrate_conc, vmax, km_eff
        )

        # ASSUMPTION: The pathway yield scales with the minimum flux
        # across all steps (linear pathway approximation). We use the
        # relative flux change to scale the ATP yield.
        nominal_rate = EnzymeKinetics.michaelis_menten(
            substrate_conc, vmax, km_nominal
        )
        relative_flux = rate / max(nominal_rate, 1e-30)

        # Scale the ATP yield by the relative flux
        base_yield = ObjectiveFunction.compute_atp_yield(self.base_params)
        return base_yield * min(relative_flux, 2.0)

    def run(self) -> EnzymeEfficiencyResult:
        """Execute the enzyme efficiency scan.

        Returns
        -------
        EnzymeEfficiencyResult with sensitivity data.
        """
        lo, hi, n = self.multiplier_range
        multipliers = np.logspace(np.log10(lo), np.log10(hi), n)

        base_yield = ObjectiveFunction.compute_atp_yield(self.base_params)
        n_enzymes = len(self.enzymes_to_test)

        atp_yields = np.zeros((n_enzymes, n))
        flux_values = np.zeros((n_enzymes, n))

        for i, enzyme_name in enumerate(self.enzymes_to_test):
            for j, mult in enumerate(multipliers):
                atp_yields[i, j] = self._compute_yield_for_efficiency(
                    enzyme_name, mult
                )

            # Relative flux as fraction of nominal
            flux_values[i, :] = atp_yields[i, :] / max(base_yield, 1e-15)

        # Bottleneck ranking: sensitivity = (yield at 10x - yield at 0.1x)
        # / log10(10/0.1)
        bottleneck_scores: List[Tuple[str, float]] = []
        for i, name in enumerate(self.enzymes_to_test):
            y_low = atp_yields[i, 0]
            y_high = atp_yields[i, -1]
            sensitivity = (y_high - y_low) / (np.log10(hi) - np.log10(lo))
            bottleneck_scores.append((name, sensitivity))

        bottleneck_scores.sort(key=lambda x: x[1])

        self.result = EnzymeEfficiencyResult(
            enzyme_names=self.enzymes_to_test,
            efficiency_multipliers=multipliers,
            atp_yields=atp_yields,
            flux_values=flux_values,
            base_yield=base_yield,
            bottleneck_ranking=bottleneck_scores,
        )
        self.result.summary = self._generate_summary()
        return self.result

    def _generate_summary(self) -> str:
        """Generate a formatted text summary."""
        if self.result is None:
            return ""
        r = self.result

        lines = [
            "=" * 84,
            "  ENZYME EFFICIENCY EXPERIMENT — PTR-94",
            "=" * 84,
            f"  Multiplier range: [{r.efficiency_multipliers[0]:.2f}x, "
            f"{r.efficiency_multipliers[-1]:.2f}x]  "
            f"({len(r.efficiency_multipliers)} steps)",
            f"  Base ATP yield: {r.base_yield:.2f}",
            "",
            "  Bottleneck Ranking (most limiting first):",
            f"  {'Enzyme':<32} {'Sensitivity':>12} {'Pathway':<15}",
            "-" * 60,
        ]

        # Most limiting (lowest sensitivity = bottleneck)
        for name, score in r.bottleneck_ranking[:6]:
            pathway = ENZYME_MODULES.get(name, {}).get("pathway", "unknown")
            lines.append(
                f"  {name:<32} {score:>+12.4f}  {pathway:<15}"
            )

        lines.extend([
            "-" * 60,
            "",
            "  Key points:",
            "  * Positive sensitivity = yield increases with efficiency.",
            "  * Enzymes with near-zero sensitivity are not rate-limiting.",
            "  * Bottlenecks are targets for protein engineering.",
            "",
            "  Yield at extreme multipliers for each enzyme:",
            f"  {'Enzyme':<32} {'0.1x':>10} {'1x':>10} {'10x':>10}",
            "-" * 64,
        ])

        mid_idx = len(r.efficiency_multipliers) // 2
        for i, name in enumerate(r.enzyme_names):
            lines.append(
                f"  {name:<32} {r.atp_yields[i, 0]:>10.2f} "
                f"{r.atp_yields[i, mid_idx]:>10.2f} "
                f"{r.atp_yields[i, -1]:>10.2f}"
            )

        lines.append("=" * 84)
        return "\n".join(lines)

    def data_table(self) -> Dict[str, Any]:
        """Return structured data for export."""
        if self.result is None:
            return {}
        return {
            "enzyme_names": self.enzymes_to_test,
            "multipliers": self.efficiency_multipliers.tolist(),
            "atp_yields": self.result.atp_yields.tolist(),
            "bottleneck_ranking": [
                (name, score)
                for name, score in self.result.bottleneck_ranking
            ],
        }


if __name__ == "__main__":
    exp = EnzymeEfficiencyExperiment()
    result = exp.run()
    print(result.summary)

    print("\n  Top 3 bottlenecks:")
    for name, score in result.bottleneck_ranking[:3]:
        print(f"  {name}: sensitivity = {score:.4f}")
