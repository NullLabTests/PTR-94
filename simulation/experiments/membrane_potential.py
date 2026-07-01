#!/usr/bin/env python3
"""
Membrane Potential Experiment.

Explores how the mitochondrial membrane potential (Δψ) affects the proton
motive force, ATP synthesis rate, and proton leakage in the PTR-94 Perfect
Coupling Module. Scans Δψ from -100 mV to -300 mV (inside negative).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple

import numpy as np

from ..thermodynamics import (
    CONST,
    proton_motive_force,
    gibbs_from_ph_and_potential,
)
from ..kinetics import ProtonLeakKinetics, ATPSynthaseKinetics
from ..pareto_optimizer import (
    DesignParameters,
    ObjectiveFunction,
    NADH_COUNT,
    FADH2_COUNT,
    SUBSTRATE_ATP,
    THEORETICAL_MAX_ATP,
)


@dataclass
class MembranePotentialResult:
    """Container for membrane potential experiment results.

    Attributes:
        delta_psi_values: Scanned Δψ values (mV or V).
        pmf_values: Proton motive force at each Δψ (kJ/mol).
        atp_synthesis_rates: ATP synthase rate at each Δψ (M/s).
        leak_rates: Proton leak rate at each Δψ (M/s).
        net_atp_flux: Net ATP flux (synthesis - leak cost) (M/s).
        atp_yields: Estimated ATP yield per glucose.
        critical_pmf: PMF where leak equals synthesis rate.
        summary: Formatted text summary.
    """
    delta_psi_values: np.ndarray
    pmf_values: np.ndarray
    atp_synthesis_rates: np.ndarray
    leak_rates: np.ndarray
    net_atp_flux: np.ndarray
    atp_yields: np.ndarray
    critical_pmf: float = 0.0
    summary: str = ""


class MembranePotentialExperiment:
    """Experiment scanning membrane potential.

    Evaluates the trade-off between ATP synthase driving force and proton
    leakage as Δψ varies. Identifies the optimal operating potential for
    the PTR-94 design.

    Attributes:
        psi_range: (low, high, n_steps) for Δψ in volts.
            Default: (-0.100, -0.300, 41) i.e. -100 to -300 mV.
        delta_ph: Fixed pH gradient (pH_out - pH_in). Default: 0.5.
        adp_conc: ADP concentration (M). Default: 1e-3.
        pi_conc: Phosphate concentration (M). Default: 5e-3.
        atp_adp_ratio: Cellular ATP/ADP ratio. Default: 10.0.
        base_params: Nominal PCM DesignParameters.
        result: MembranePotentialResult after run().
    """

    def __init__(
        self,
        psi_range: Tuple[float, float, int] = (-0.100, -0.300, 41),
        delta_ph: float = 0.5,
        adp_conc: float = 1e-3,
        pi_conc: float = 5e-3,
        atp_adp_ratio: float = 10.0,
        base_params: DesignParameters | None = None,
    ) -> None:
        """
        Parameters
        ----------
        psi_range : tuple
            (low, high, n_steps) for Δψ in volts.
        delta_ph : float
            Transmembrane pH gradient (default: 0.5).
        adp_conc : float
            ADP concentration in molar (default: 1e-3).
        pi_conc : float
            Phosphate concentration in molar (default: 5e-3).
        atp_adp_ratio : float
            Cellular ATP/ADP ratio (default: 10.0).
        base_params : DesignParameters, optional
            Baseline PCM parameters.
        """
        self.psi_range = psi_range
        self.delta_ph = delta_ph
        self.adp_conc = adp_conc
        self.pi_conc = pi_conc
        self.atp_adp_ratio = atp_adp_ratio
        self.base_params = base_params or DesignParameters()
        self.result: MembranePotentialResult | None = None

        # ASSUMPTION: The pH gradient is held constant as Δψ varies. In real
        # cells, ΔpH and Δψ are interconvertible (Mitchell's chemiosmotic
        # theory), but for this parametric scan we treat ΔpH as fixed.

    def run(self) -> MembranePotentialResult:
        """Execute the membrane potential scan.

        Returns
        -------
        MembranePotentialResult with arrays over the Δψ range.
        """
        lo, hi, n = self.psi_range
        # Make the scan run from low magnitude to high magnitude
        psi_values = np.linspace(lo, hi, n)

        pmf_vals = np.zeros(n)
        atp_rates = np.zeros(n)
        leak_vals = np.zeros(n)
        net_flux = np.zeros(n)
        atp_yields = np.zeros(n)

        # Use absolute membrane potential for proton motive force calculation
        # (Δψ is negative for mitochondria: inside negative)
        # ASSUMPTION: Standard physiological temperature.
        temp = CONST.T_phys

        # Proton leak kinetics
        # ASSUMPTION: PCM membrane has lower basal conductance than natural.
        leak_model = ProtonLeakKinetics(
            l0=1e-9,
            alpha=0.06,
            q10=1.5,
            temp_ref=temp,
        )

        # ATP synthase kinetics
        atp_synthase = ATPSynthaseKinetics(
            h_per_atp=self.base_params.h_per_atp,
            k_rotation=150.0,
            km_atp_site=1e-5,
            km_pi_site=1e-3,
            slip_probability=self.base_params.proton_slip_probability,
        )

        for i, psi in enumerate(psi_values):
            # PMF in kJ/mol (absolute value since Δψ is negative)
            psi_abs = abs(psi)
            pmf = proton_motive_force(psi_abs, self.delta_ph, temp)
            pmf_vals[i] = pmf

            # ATP synthesis rate
            rate = atp_synthase.rotation_rate(
                pmf_kj_per_mol=pmf,
                adp_conc=self.adp_conc,
                pi_conc=self.pi_conc,
                atp_adp_ratio=self.atp_adp_ratio,
            )
            atp_rates[i] = rate

            # Proton leak rate
            leak = leak_model.leak_rate(pmf, temp)
            leak_vals[i] = leak

            # Net ATP flux reduces by the fraction lost to leak
            # ASSUMPTION: Leak H+ flux reduces net ATP synthesis proportionally
            # to the H+/ATP stoichiometry.
            h_per_atp = self.base_params.h_per_atp
            leak_cost_in_atp = leak / max(h_per_atp, 1.0)
            net_flux[i] = max(rate - leak_cost_in_atp, 0.0)

        # Estimate ATP yield per glucose (after loop, arrays fully populated)
        max_net = max(net_flux.max(), 1e-30)
        reference_yield = min(
            SUBSTRATE_ATP
            + (self.base_params.h_per_nadh * NADH_COUNT
               + self.base_params.h_per_fadh2 * FADH2_COUNT)
            / self.base_params.h_per_atp,
            THEORETICAL_MAX_ATP,
        )
        atp_yields = net_flux / max_net * reference_yield

        # Find critical PMF where leak == synthesis rate
        diff = np.abs(atp_rates - leak_vals)
        crit_idx = int(np.argmin(diff))
        critical_pmf = float(pmf_vals[crit_idx])

        idx_best = int(np.argmax(net_flux))

        self.result = MembranePotentialResult(
            delta_psi_values=psi_values,
            pmf_values=pmf_vals,
            atp_synthesis_rates=atp_rates,
            leak_rates=leak_vals,
            net_atp_flux=net_flux,
            atp_yields=atp_yields,
            critical_pmf=critical_pmf,
        )
        self.result.summary = self._generate_summary(idx_best)
        return self.result

    def _generate_summary(self, best_idx: int) -> str:
        """Generate formatted results table."""
        if self.result is None:
            return ""
        r = self.result

        lines = [
            "=" * 84,
            "  MEMBRANE POTENTIAL EXPERIMENT — PTR-94",
            "=" * 84,
            f"  Δψ range: [{r.delta_psi_values[0]*1000:.0f}, "
            f"{r.delta_psi_values[-1]*1000:.0f}] mV",
            f"  ΔpH (fixed): {self.delta_ph:.2f}",
            f"  Temperature: {CONST.T_phys:.1f} K",
            "-" * 84,
            f"  {'Δψ(mV)':>8} {'PMF(kJ)':>9} {'ATP rate':>9} "
            f"{'Leak':>9} {'Net ATP':>9} {'Yield':>7}",
            f"  {'':>8} {'/mol':>9} {'(M/s)':>9} {'(M/s)':>9} "
            f"{'(M/s)':>9} {'':>7}",
            "-" * 84,
        ]

        # Show selected points: first, every 10th, best, last
        n = len(r.delta_psi_values)
        indices = {0, best_idx, n - 1}
        step = max(n // 10, 1)
        indices.update(range(0, n, step))
        for idx in sorted(indices):
            lines.append(
                f"  {r.delta_psi_values[idx]*1000:>8.0f} "
                f"{r.pmf_values[idx]:>9.2f} "
                f"{r.atp_synthesis_rates[idx]:>9.2e} "
                f"{r.leak_rates[idx]:>9.2e} "
                f"{r.net_atp_flux[idx]:>9.2e} "
                f"{r.atp_yields[idx]:>7.1f}"
            )

        lines.append("-" * 84)

        best_psi = r.delta_psi_values[best_idx] * 1000
        lines.append(f"\n  Optimal Δψ:           {best_psi:.0f} mV")
        lines.append(f"  PMF at optimum:       {r.pmf_values[best_idx]:.2f} kJ/mol")
        lines.append(f"  Max net ATP flux:     {r.net_atp_flux[best_idx]:.2e} M/s")
        lines.append(f"  ATP yield at optimum: {r.atp_yields[best_idx]:.1f} per glucose")
        lines.append(f"  Critical PMF (leak=synth): {r.critical_pmf:.1f} kJ/mol")

        lines.append("")
        lines.append("  ASSUMPTION: The pH gradient is fixed at "
                      f"{self.delta_ph}. Real cells")
        lines.append("  interconvert Δψ and ΔpH via ion transport.")
        lines.append("=" * 84)
        return "\n".join(lines)

    def data_table(self) -> Dict[str, List[float]]:
        """Return structured data for export."""
        if self.result is None:
            return {}
        r = self.result
        return {
            "delta_psi_mV": (r.delta_psi_values * 1000).tolist(),
            "pmf_kj_per_mol": r.pmf_values.tolist(),
            "atp_synthesis_rate_M_per_s": r.atp_synthesis_rates.tolist(),
            "leak_rate_M_per_s": r.leak_rates.tolist(),
            "net_atp_flux_M_per_s": r.net_atp_flux.tolist(),
            "atp_yield": r.atp_yields.tolist(),
        }


if __name__ == "__main__":
    exp = MembranePotentialExperiment()
    result = exp.run()
    print(result.summary)

    print("\n  Key data points:")
    table = exp.data_table()
    for i in [0, len(table["delta_psi_mV"]) // 2, -1]:
        psi = table["delta_psi_mV"][i]
        pmf = table["pmf_kj_per_mol"][i]
        atp = table["atp_yield"][i]
        print(f"  Δψ={psi:.0f} mV  PMF={pmf:.2f} kJ/mol  Yield={atp:.1f} ATP")
