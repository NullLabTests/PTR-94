#!/usr/bin/env python3
"""
Leakage Sensitivity Experiment.

Explores how membrane proton leakage affects ATP yield in the PTR-94
Perfect Coupling Module. Scans leakage conductance from 0 to 0.5
(dimensionless), computing ATP yield, thermodynamic efficiency, and
entropy production at each level.
"""

from __future__ import annotations

import math
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
)


@dataclass
class LeakageSensitivityResult:
    """Container for results from a leakage sensitivity scan.

    Attributes:
        leakage_conductances: Array of scanned leakage values.
        atp_yields: ATP yield at each leakage level.
        efficiencies: Thermodynamic efficiency (percent) at each level.
        entropies: Entropy production (J/(K*mol)) at each level.
        proton_fluxes: Effective proton flux (H+ per glucose).
        redox_atps: Redox-derived ATP per glucose.
        summary: Text summary of findings.
    """
    leakage_conductances: np.ndarray
    atp_yields: np.ndarray
    efficiencies: np.ndarray
    entropies: np.ndarray
    proton_fluxes: np.ndarray
    redox_atps: np.ndarray
    summary: str = ""


class LeakageSensitivityExperiment:
    """Experiment scanning membrane proton leakage conductance.

    Evaluates the impact of passive proton leak on the ATP yield,
    efficiency, and entropy production of the PTR-94 design.

    Attributes:
        base_params: Nominal DesignParameters (leak will be overridden).
        conductance_range: (low, high, n_steps) for the scan.
        result: LeakageSensitivityResult after run().
    """

    def __init__(
        self,
        base_params: DesignParameters | None = None,
        conductance_range: Tuple[float, float, int] = (0.0, 0.5, 51),
    ) -> None:
        """
        Parameters
        ----------
        base_params : DesignParameters, optional
            Baseline parameters. Defaults to nominal PTR-94 design.
        conductance_range : tuple
            (low, high, n_steps) for leakage conductance scan.
        """
        self.base_params = base_params or DesignParameters()
        self.conductance_range = conductance_range
        self.result: LeakageSensitivityResult | None = None

        # ASSUMPTION: Nominal PTR-94 design has membrane_leak_conductance = 0.05.
        # The scan range covers ideal (0) through extreme (0.5) leakage.
        # ASSUMPTION: Leak reduces effective proton flux linearly.

    def run(self) -> LeakageSensitivityResult:
        """Execute the leakage sensitivity scan.

        Returns
        -------
        LeakageSensitivityResult with arrays over the conductance range.
        """
        lo, hi, n = self.conductance_range
        conductances = np.linspace(lo, hi, n)

        atp_yields = np.zeros(n)
        efficiencies = np.zeros(n)
        entropies = np.zeros(n)
        proton_fluxes = np.zeros(n)
        redox_atps = np.zeros(n)

        for i, leak in enumerate(conductances):
            params = DesignParameters(
                h_per_nadh=self.base_params.h_per_nadh,
                h_per_fadh2=self.base_params.h_per_fadh2,
                h_per_atp=self.base_params.h_per_atp,
                atp_synthase_efficiency=self.base_params.atp_synthase_efficiency,
                membrane_leak_conductance=float(leak),
                proton_slip_probability=self.base_params.proton_slip_probability,
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

            # ASSUMPTION: Effective proton flux is total pumped H+ scaled by
            # the same loss factors used in the yield calculation.
            protons_nadh = params.h_per_nadh * NADH_COUNT
            protons_fadh2 = params.h_per_fadh2 * FADH2_COUNT
            total_pumped = protons_nadh + protons_fadh2
            effective = (
                total_pumped
                * (1.0 - leak)
                * (1.0 - params.proton_slip_probability)
                * params.roq_quinone_coupling
            )
            proton_fluxes[i] = effective
            redox_atps[i] = atp - SUBSTRATE_ATP

        self.result = LeakageSensitivityResult(
            leakage_conductances=conductances,
            atp_yields=atp_yields,
            efficiencies=efficiencies,
            entropies=entropies,
            proton_fluxes=proton_fluxes,
            redox_atps=redox_atps,
        )
        self.result.summary = self._generate_summary()
        return self.result

    def _generate_summary(self) -> str:
        """Generate a text summary of the experiment results."""
        if self.result is None:
            return "No results available."
        r = self.result
        i_zero = 0
        i_mid = len(r.leakage_conductances) // 2
        i_end = -1

        lines = [
            "=" * 74,
            "  LEAKAGE SENSITIVITY EXPERIMENT — PTR-94",
            "=" * 74,
            f"  Conductance range: [{r.leakage_conductances[0]:.3f}, "
            f"{r.leakage_conductances[-1]:.3f}]  ({len(r.leakage_conductances)} steps)",
            "-" * 74,
            f"  {'Leak':>8}  {'ATP yield':>10}  {'Efficiency':>10}  {'Entropy':>10}",
            f"  {'':>8}  {'':>10}  {'(%)':>10}  {'J/(K*mol)':>10}",
            "-" * 74,
        ]

        for idx in [i_zero, i_mid, i_end]:
            lines.append(
                f"  {r.leakage_conductances[idx]:>8.4f}  "
                f"{r.atp_yields[idx]:>10.2f}  "
                f"{r.efficiencies[idx]:>10.2f}  "
                f"{r.entropies[idx]:>10.2f}"
            )

        lines.append("-" * 74)

        # Key findings
        loss_zero = r.atp_yields[i_zero]
        loss_end = r.atp_yields[i_end]
        pct_loss = (loss_zero - loss_end) / loss_zero * 100.0

        lines.append(f"\n  Nominal yield (leak=0):          {loss_zero:.2f} ATP")
        lines.append(f"  Yield at max leak ({r.leakage_conductances[i_end]:.2f}):     {loss_end:.2f} ATP")
        lines.append(f"  Total loss:                      {loss_zero - loss_end:.2f} ATP ({pct_loss:.1f}%)")
        lines.append(f"  Yield halved at conductance:     {self._find_halving_point():.4f}")

        crit = self._critical_leak()
        lines.append(f"  Critical leak (yield <= {SUBSTRATE_ATP}):      {crit:.4f}")

        lines.append("=" * 74)
        return "\n".join(lines)

    def _find_halving_point(self) -> float:
        """Find the leakage conductance where yield drops by half."""
        if self.result is None:
            return float("nan")
        r = self.result
        half = r.atp_yields[0] / 2.0
        idx = int(np.argmin(np.abs(r.atp_yields - half)))
        return float(r.leakage_conductances[idx])

    def _critical_leak(self) -> float:
        """Find the leakage where yield falls to substrate-level ATP only."""
        if self.result is None:
            return float("nan")
        r = self.result
        critical = SUBSTRATE_ATP + 0.5
        idx = int(np.argmin(np.abs(r.atp_yields - critical)))
        return float(r.leakage_conductances[idx])

    def data_table(self) -> Dict[str, List[float]]:
        """Return structured data for export or further analysis.

        Returns
        -------
        dict with keys: 'leakage', 'atp_yield', 'efficiency',
                         'entropy', 'proton_flux', 'redox_atp'.
        """
        if self.result is None:
            return {}
        r = self.result
        return {
            "leakage": r.leakage_conductances.tolist(),
            "atp_yield": r.atp_yields.tolist(),
            "efficiency": r.efficiencies.tolist(),
            "entropy": r.entropies.tolist(),
            "proton_flux": r.proton_fluxes.tolist(),
            "redox_atp": r.redox_atps.tolist(),
        }


if __name__ == "__main__":
    exp = LeakageSensitivityExperiment()
    result = exp.run()
    print(result.summary)

    print("\n  Data table (first 5, last 2 rows):")
    table = exp.data_table()
    keys = list(table.keys())
    print(f"  {'  '.join(f'{k:>12}' for k in keys)}")
    for i in list(range(5)) + [len(table["leakage"]) - 2, len(table["leakage"]) - 1]:
        if i < len(table["leakage"]):
            print(f"  {'  '.join(f'{table[k][i]:>12.4f}' for k in keys)}")
