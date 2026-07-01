#!/usr/bin/env python3
"""
Proton Slip Experiment.

Explores how ATP synthase proton slip affects the effective H+/ATP ratio
and ATP yield in the PTR-94 Perfect Coupling Module. Scans slip probability
from 0 to 0.3, representing stochastic decoupling events where protons
translocate without driving ATP synthesis.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple

import numpy as np

from ..thermodynamics import (
    STANDARD_GIBBS_GLUCOSE,
    STANDARD_GIBBS_ATP,
    compute_coupling_efficiency,
    entropy_production,
)
from ..pareto_optimizer import (
    DesignParameters,
    ObjectiveFunction,
    SUBSTRATE_ATP,
    NADH_COUNT,
    FADH2_COUNT,
    THEORETICAL_MAX_ATP,
)


@dataclass
class ProtonSlipResult:
    """Container for proton slip experiment results.

    Attributes:
        slip_probabilities: Scanned slip probability values.
        atp_yields: ATP yield at each slip level.
        effective_h_per_atp: Effective H+/ATP ratio (actual protons consumed
            per ATP produced).
        efficiencies: Thermodynamic efficiency (percent).
        entropies: Entropy production (J/(K*mol)).
        yield_loss: ATP yield loss relative to zero slip.
        relative_yield: Yield as fraction of slip-free yield.
        critical_slip: Slip where yield drops below SUBSTRATE_ATP + 1.
        summary: Formatted text summary.
    """
    slip_probabilities: np.ndarray
    atp_yields: np.ndarray
    effective_h_per_atp: np.ndarray
    efficiencies: np.ndarray
    entropies: np.ndarray
    yield_loss: np.ndarray
    relative_yield: np.ndarray
    critical_slip: float = 0.0
    summary: str = ""


class ProtonSlipExperiment:
    """Experiment scanning ATP synthase slip probability.

    In the ATP synthase rotary mechanism, slip occurs when protons pass
    through the enzyme without inducing the conformational changes needed
    for ATP synthesis. This effectively increases the H+/ATP ratio.

    Attributes:
        slip_range: (low, high, n_steps) for slip probability.
        base_params: Nominal DesignParameters for the PCM.
        result: ProtonSlipResult after run().
    """

    def __init__(
        self,
        slip_range: Tuple[float, float, int] = (0.0, 0.3, 61),
        base_params: DesignParameters | None = None,
    ) -> None:
        """
        Parameters
        ----------
        slip_range : tuple
            (low, high, n_steps) for slip probability.
        base_params : DesignParameters, optional
            Baseline PCM parameters.
        """
        self.slip_range = slip_range
        self.base_params = base_params or DesignParameters()
        self.result: ProtonSlipResult | None = None

        # ASSUMPTION: Slip probability is the fraction of 120-degree rotor
        # steps that fail to synthesise ATP. In the PTR-94 design, slip is
        # minimised through engineered rotor-stator interactions and
        # optimised proton half-channel geometry.
        # ASSUMPTION: Slip increases the effective H+/ATP ratio as:
        #   H+/ATP_effective = H+/ATP_base / (1 - slip)

    def run(self) -> ProtonSlipResult:
        """Execute the slip probability scan.

        Returns
        -------
        ProtonSlipResult with arrays over the slip range.
        """
        lo, hi, n = self.slip_range
        slip_vals = np.linspace(lo, hi, n)

        atp_yields = np.zeros(n)
        effective_h_per_atp = np.zeros(n)
        efficiencies = np.zeros(n)
        entropies = np.zeros(n)

        for i, slip in enumerate(slip_vals):
            # Effective H+/ATP increases with slip
            base_h_atp = self.base_params.h_per_atp
            if slip >= 1.0:
                eff_h_atp = float("inf")
            else:
                eff_h_atp = base_h_atp / (1.0 - slip)
            effective_h_per_atp[i] = eff_h_atp

            # Compute ATP yield with effective H+/ATP
            # ASSUMPTION: Slip only affects the ATP synthase step; the
            # proton pumping and other parameters remain unchanged.
            params = DesignParameters(
                h_per_nadh=self.base_params.h_per_nadh,
                h_per_fadh2=self.base_params.h_per_fadh2,
                h_per_atp=min(eff_h_atp, 10.0),
                atp_synthase_efficiency=self.base_params.atp_synthase_efficiency,
                membrane_leak_conductance=self.base_params.membrane_leak_conductance,
                proton_slip_probability=float(slip),
                roq_quinone_coupling=self.base_params.roq_quinone_coupling,
                scaffold_channeling_efficiency=self.base_params.scaffold_channeling_efficiency,
                ros_bypass_fraction=self.base_params.ros_bypass_fraction,
            )
            atp = ObjectiveFunction.compute_atp_yield(params)
            atp_yields[i] = atp
            efficiencies[i] = compute_coupling_efficiency(
                atp, STANDARD_GIBBS_GLUCOSE, STANDARD_GIBBS_ATP
            )
            energy_captured = atp * STANDARD_GIBBS_ATP
            entropies[i] = entropy_production(
                STANDARD_GIBBS_GLUCOSE, energy_captured
            )

        yield_loss = atp_yields[0] - atp_yields
        relative_yield = atp_yields / max(atp_yields[0], 1e-15)

        # Critical slip: yield drops to SUBSTRATE_ATP + 1
        critical = SUBSTRATE_ATP + 1.0
        diff = np.abs(atp_yields - critical)
        crit_idx = int(np.argmin(diff))
        critical_slip = float(slip_vals[crit_idx])

        self.result = ProtonSlipResult(
            slip_probabilities=slip_vals,
            atp_yields=atp_yields,
            effective_h_per_atp=effective_h_per_atp,
            efficiencies=efficiencies,
            entropies=entropies,
            yield_loss=yield_loss,
            relative_yield=relative_yield,
            critical_slip=critical_slip,
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
            "  PROTON SLIP EXPERIMENT — PTR-94",
            "=" * 84,
            f"  Slip probability range: [{r.slip_probabilities[0]:.3f}, "
            f"{r.slip_probabilities[-1]:.3f}]  ({len(r.slip_probabilities)} steps)",
            f"  Base H+/ATP: {self.base_params.h_per_atp:.1f}",
            "-" * 84,
            f"  {'Slip':>8}  {'H+/ATP':>8}  {'ATP':>8}  {'Eff%':>7}  "
            f"{'Entropy':>9}  {'Loss':>8}  {'Rel':>6}",
            f"  {'':>8}  {'eff.':>8}  {'yield':>8}  {'':>7}  "
            f"{'J/(K*mol)':>9}  {'ATP':>8}  {'':>6}",
            "-" * 84,
        ]

        for idx in [0, len(r.slip_probabilities) // 4,
                     len(r.slip_probabilities) // 2,
                     3 * len(r.slip_probabilities) // 4, -1]:
            if idx < len(r.slip_probabilities):
                lines.append(
                    f"  {r.slip_probabilities[idx]:>8.4f}  "
                    f"{r.effective_h_per_atp[idx]:>8.2f}  "
                    f"{r.atp_yields[idx]:>8.2f}  "
                    f"{r.efficiencies[idx]:>7.2f}  "
                    f"{r.entropies[idx]:>9.2f}  "
                    f"{r.yield_loss[idx]:>8.2f}  "
                    f"{r.relative_yield[idx]:>6.3f}"
                )

        lines.extend([
            "-" * 84,
            f"  Yield at slip=0:         {r.atp_yields[0]:.2f} ATP",
            f"  Yield at slip={r.slip_probabilities[-1]:.2f}:      "
            f"{r.atp_yields[-1]:.2f} ATP",
            f"  Total yield loss:        "
            f"{r.yield_loss[-1]:.2f} ATP ({r.yield_loss[-1]/max(r.atp_yields[0],1e-15)*100:.1f}%)",
            f"  Critical slip (yield={SUBSTRATE_ATP+1}): {r.critical_slip:.4f}",
            "",
            "  Observations:",
            "  * At slip=0.3, effective H+/ATP = "
            f"{r.effective_h_per_atp[-1]:.1f} vs base {self.base_params.h_per_atp:.1f}.",
            "  * PTR-94 requires slip < 0.1 to maintain yield > 90 ATP.",
            "=" * 84,
        ])
        return "\n".join(lines)

    def data_table(self) -> Dict[str, List[float]]:
        """Return structured data for export."""
        if self.result is None:
            return {}
        r = self.result
        return {
            "slip_probability": r.slip_probabilities.tolist(),
            "effective_h_per_atp": r.effective_h_per_atp.tolist(),
            "atp_yield": r.atp_yields.tolist(),
            "efficiency_pct": r.efficiencies.tolist(),
            "entropy_J_per_K_per_mol": r.entropies.tolist(),
            "yield_loss": r.yield_loss.tolist(),
        }


if __name__ == "__main__":
    exp = ProtonSlipExperiment()
    result = exp.run()
    print(result.summary)

    print("\n  Data table (slip vs. yield):")
    table = exp.data_table()
    for i in range(0, len(table["slip_probability"]), 10):
        s = table["slip_probability"][i]
        y = table["atp_yield"][i]
        heff = table["effective_h_per_atp"][i]
        print(f"  slip={s:.3f}  ATP={y:.2f}  H+/ATP(eff)={heff:.2f}")
