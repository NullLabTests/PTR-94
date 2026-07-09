#!/usr/bin/env python3
"""
ATP Free Energy Experiment.

Explores how the Gibbs free energy of ATP synthesis (ΔG_ATP) affects
the theoretical maximum ATP yield from glucose oxidation. Scans ΔG_ATP
from 30 to 60 kJ/mol, covering standard biochemical to physiological
conditions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

from ..thermodynamics import (
    STANDARD_GIBBS_GLUCOSE,
    STANDARD_GIBBS_ATP,
    PHYSIOLOGICAL_GIBBS_ATP,
    compute_theoretical_max_atp,
    compute_coupling_efficiency,
    thermodynamic_ceiling,
)


@dataclass
class ATPFreeEnergyResult:
    """Container for ATP free energy experiment results.

    Attributes:
        delta_g_atp_values: Scanned ΔG_ATP values (kJ/mol).
        theoretical_max_atp: Continuous thermodynamic max at each value.
        integer_ceiling: Integer floor of max ATP at each value.
        efficiencies: Thermodynamic efficiency (percent) assuming 94 ATP.
        ptr94_achievable: Whether 94 ATP is thermodynamically feasible.
        breakpoint_dg: ΔG_ATP where max falls below 94.
        summary: Text summary.
    """
    delta_g_atp_values: np.ndarray
    theoretical_max_atp: np.ndarray
    integer_ceiling: np.ndarray
    efficiencies: np.ndarray
    ptr94_achievable: np.ndarray
    breakpoint_dg: float = 0.0
    summary: str = ""


class ATPFreeEnergyExperiment:
    """Experiment scanning ATP synthesis free energy.

    The theoretical maximum ATP yield equals the free energy of glucose
    oxidation divided by the free energy required per ATP. As ΔG_ATP
    increases (more physiological), the theoretical yield drops.

    Attributes:
        dg_range: (low, high, n_steps) for ΔG_ATP in kJ/mol.
        result: ATPFreeEnergyResult after run().
    """

    def __init__(
        self,
        dg_range: Tuple[float, float, int] = (30.0, 60.0, 101),
    ) -> None:
        """
        Parameters
        ----------
        dg_range : tuple
            (low, high, n_steps) for ΔG_ATP scan in kJ/mol.
            30 kJ/mol = standard biochemical.
            55-60 kJ/mol = physiological.
        """
        self.dg_range = dg_range

        # ASSUMPTION: Standard glucose oxidation free energy (2870 kJ/mol)
        # is used throughout. In vivo values may vary slightly.
        # ASSUMPTION: We use the standard value ΔG_glucose = 2870 kJ/mol.

    def run(self) -> ATPFreeEnergyResult:
        """Execute the ΔG_ATP scan.

        Returns
        -------
        ATPFreeEnergyResult with arrays over the free energy range.
        """
        lo, hi, n = self.dg_range
        dg_values = np.linspace(lo, hi, n)

        theoretical_max = np.array([
            compute_theoretical_max_atp(
                STANDARD_GIBBS_GLUCOSE, float(dg)
            )
            for dg in dg_values
        ])
        integer_ceilings = np.array([
            thermodynamic_ceiling(STANDARD_GIBBS_GLUCOSE, float(dg))
            for dg in dg_values
        ])
        efficiencies = np.array([
            compute_coupling_efficiency(94.0, STANDARD_GIBBS_GLUCOSE, float(dg))
            for dg in dg_values
        ])
        ptr94_possible = integer_ceilings >= 94

        # ASSUMPTION: The PTR-94 target of 94 ATP is achievable when
        # the integer ceiling >= 94.
        breakpoint_idx = int(np.argmax(integer_ceilings < 94))
        if integer_ceilings[breakpoint_idx] < 94:
            breakpoint_dg = float(dg_values[breakpoint_idx])
        else:
            breakpoint_dg = float(dg_values[-1])

        self.result = ATPFreeEnergyResult(
            delta_g_atp_values=dg_values,
            theoretical_max_atp=theoretical_max,
            integer_ceiling=integer_ceilings,
            efficiencies=efficiencies,
            ptr94_achievable=ptr94_possible,
            breakpoint_dg=breakpoint_dg,
        )
        self.result.summary = self._generate_summary()
        return self.result

    def _generate_summary(self) -> str:
        """Generate a formatted text summary."""
        if self.result is None:
            return ""
        r = self.result

        lines = [
            "=" * 74,
            "  ATP FREE ENERGY EXPERIMENT — PTR-94",
            "=" * 74,
            f"  ΔG_ATP range: [{r.delta_g_atp_values[0]:.1f}, "
            f"{r.delta_g_atp_values[-1]:.1f}] kJ/mol  "
            f"({len(r.delta_g_atp_values)} steps)",
            f"  ΔG_glucose (fixed): {STANDARD_GIBBS_GLUCOSE:.0f} kJ/mol",
            "-" * 74,
            f"  {'ΔG_ATP':>8}  {'Max ATP':>10}  {'Floor':>6}  {'Eff(94)':>8}  {'PTR-94':>7}",
            f"  {'(kJ/mol)':>8}  {'':>10}  {'':>6}  {'(%)':>8}  {'feas?':>7}",
            "-" * 74,
        ]

        # Show key points: standard, physiological, breakpoint, ends
        dg = r.delta_g_atp_values
        key_indices = [0]
        # Find closest to 30.5 (standard), 45, 55 (physiological)
        for target in [STANDARD_GIBBS_ATP, 45.0, PHYSIOLOGICAL_GIBBS_ATP]:
            idx = int(np.argmin(np.abs(dg - target)))
            key_indices.append(idx)
        key_indices.append(len(dg) - 1)
        key_indices = sorted(set(key_indices))

        for idx in key_indices:
            lines.append(
                f"  {r.delta_g_atp_values[idx]:>8.1f}  "
                f"{r.theoretical_max_atp[idx]:>10.2f}  "
                f"{r.integer_ceiling[idx]:>6d}  "
                f"{r.efficiencies[idx]:>8.2f}  "
                f"{'Yes' if r.ptr94_achievable[idx] else 'No':>7}"
            )

        lines.append("-" * 74)

        # Standard conditions
        std_idx = int(np.argmin(np.abs(dg - STANDARD_GIBBS_ATP)))
        lines.append(f"\n  Standard (ΔG_ATP = {STANDARD_GIBBS_ATP:.1f} kJ/mol):")
        lines.append(f"    Max ATP = {r.theoretical_max_atp[std_idx]:.2f}")
        lines.append(f"    Floor  = {r.integer_ceiling[std_idx]} ATP")

        # Physiological
        phys_idx = int(np.argmin(np.abs(dg - PHYSIOLOGICAL_GIBBS_ATP)))
        lines.append(f"\n  Physiological (ΔG_ATP = {PHYSIOLOGICAL_GIBBS_ATP:.1f} kJ/mol):")
        lines.append(f"    Max ATP = {r.theoretical_max_atp[phys_idx]:.2f}")
        lines.append(f"    Floor  = {r.integer_ceiling[phys_idx]} ATP")
        lines.append(f"    Efficiency (at 94 ATP): {r.efficiencies[phys_idx]:.2f}%")
        # ASSUMPTION: Under physiological conditions, 94 ATP would represent
        # >100% efficiency, indicating it is not thermodynamically feasible.

        lines.append(f"\n  Breakpoint ΔG_ATP (ceiling drops below 94): {r.breakpoint_dg:.1f} kJ/mol")
        lines.append("  Below this, PTR-94 target becomes thermodynamically infeasible.")
        lines.append("=" * 74)
        return "\n".join(lines)

    def data_table(self) -> Dict[str, List[float]]:
        """Return structured data for export."""
        if self.result is None:
            return {}
        r = self.result
        return {
            "delta_g_atp": r.delta_g_atp_values.tolist(),
            "theoretical_max_atp": r.theoretical_max_atp.tolist(),
            "integer_ceiling": [int(x) for x in r.integer_ceiling],
            "efficiency_pct": r.efficiencies.tolist(),
            "ptr94_achievable": [bool(x) for x in r.ptr94_achievable],
        }


if __name__ == "__main__":
    exp = ATPFreeEnergyExperiment()
    result = exp.run()
    print(result.summary)

    print("\n  Data table (key points):")
    table = exp.data_table()
    dg = table["delta_g_atp"]
    for target in [30.5, 45.0, 55.0]:
        idx = min(range(len(dg)), key=lambda i: abs(dg[i] - target))
        print(f"  ΔG={dg[idx]:.1f}:  max={table['theoretical_max_atp'][idx]:.2f}  "
              f"floor={table['integer_ceiling'][idx]}  "
              f"eff={table['efficiency_pct'][idx]:.2f}%  "
              f"ptr94={'yes' if table['ptr94_achievable'][idx] else 'no'}")
