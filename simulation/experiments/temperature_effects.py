#!/usr/bin/env python3
"""
Temperature Effects Experiment.

Explores how temperature affects reaction rates and ATP yield in the
PTR-94 pathway. Scans temperature from 280 K to 330 K, computing Q10
effects on enzyme kinetics, proton leak, and overall ATP production.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple

import numpy as np

from ..thermodynamics import (
    CONST,
    STANDARD_GIBBS_GLUCOSE,
    STANDARD_GIBBS_ATP,
    compute_coupling_efficiency,
    compute_theoretical_max_atp,
    proton_motive_force,
)
from ..kinetics import (
    ProtonLeakKinetics,
    ATPSynthaseKinetics,
    EnzymeKinetics,
)
from ..pareto_optimizer import (
    DesignParameters,
    ObjectiveFunction,
    NADH_COUNT,
    FADH2_COUNT,
    SUBSTRATE_ATP,
    THEORETICAL_MAX_ATP,
)


# Standard Q10 values for different process types
# ASSUMPTION: Q10 values are from biological rate-temperature dependence.
# Enzymatic reactions typically have Q10 ~ 2-3, membrane processes ~ 1.5.
Q10_ENZYME: float = 2.0
Q10_LEAK: float = 1.5
Q10_DIFFUSION: float = 1.3


def q10_factor(q10: float, temp: float, temp_ref: float = CONST.T_phys) -> float:
    """Compute the Q10 temperature adjustment factor.

    Parameters
    ----------
    q10 : float
        Q10 coefficient (rate ratio per 10 K rise).
    temp : float
        Current temperature (K).
    temp_ref : float
        Reference temperature (K).

    Returns
    -------
    float
        Multiplicative factor for the rate.
    """
    return q10 ** ((temp - temp_ref) / 10.0)


@dataclass
class TemperatureEffectsResult:
    """Container for temperature effects experiment results.

    Attributes:
        temperatures: Scanned temperature values (K).
        atp_yields: ATP yield at each temperature.
        efficiencies: Thermodynamic efficiency (percent).
        theoretical_max: Thermodynamic maximum ATP at each temperature.
        q10_enzyme_factors: Q10 multiplier for enzyme rates.
        q10_leak_factors: Q10 multiplier for leak rates.
        atp_synthase_rates: ATP synthase rate (M/s) at each temperature.
        leak_rates: Proton leak rate (M/s) at each temperature.
        optimum_temp: Temperature with highest ATP yield (K).
        yield_at_optimum: ATP yield at optimal temperature.
        summary: Formatted text summary.
    """
    temperatures: np.ndarray
    atp_yields: np.ndarray
    efficiencies: np.ndarray
    theoretical_max: np.ndarray
    q10_enzyme_factors: np.ndarray
    q10_leak_factors: np.ndarray
    atp_synthase_rates: np.ndarray
    leak_rates: np.ndarray
    optimum_temp: float = 0.0
    yield_at_optimum: float = 0.0
    summary: str = ""


class TemperatureEffectsExperiment:
    """Experiment scanning temperature.

    Evaluates the effect of temperature on every component of the PTR-94
    pathway: enzyme kinetics, ATP synthase rotation, proton leakage, and
    the thermodynamic driving force (PMF).

    Attributes:
        temp_range: (low, high, n_steps) for temperature in Kelvin.
        base_params: Nominal DesignParameters.
        result: TemperatureEffectsResult after run().
    """

    def __init__(
        self,
        temp_range: Tuple[float, float, int] = (280.0, 330.0, 51),
        base_params: DesignParameters | None = None,
    ) -> None:
        """
        Parameters
        ----------
        temp_range : tuple
            (low, high, n_steps) for temperature in Kelvin.
        base_params : DesignParameters, optional
            Baseline PCM parameters.
        """
        self.temp_range = temp_range
        self.base_params = base_params or DesignParameters()
        self.result: TemperatureEffectsResult | None = None

        # ASSUMPTION: The thermodynamic ceiling changes slightly with
        # temperature due to changes in ΔG of glucose oxidation and ATP
        # synthesis. We approximate this as a linear correction.
        # ASSUMPTION: All Q10 factors remain constant across the tested
        # range. In reality, Q10 itself is temperature-dependent.

    def run(self) -> TemperatureEffectsResult:
        """Execute the temperature scan.

        Returns
        -------
        TemperatureEffectsResult with arrays over temperature range.
        """
        lo, hi, n = self.temp_range
        temps = np.linspace(lo, hi, n)

        atp_yields = np.zeros(n)
        efficiencies = np.zeros(n)
        theoretical_max = np.zeros(n)
        q10_enz = np.zeros(n)
        q10_lk = np.zeros(n)
        atp_rates = np.zeros(n)
        leak_vals = np.zeros(n)

        # Reference conditions
        delta_psi = 0.180  # V (absolute mitochondrial Δψ)
        delta_ph = 0.5
        adp_conc = 1e-3
        pi_conc = 5e-3
        atp_adp_ratio = 10.0

        for i, temp in enumerate(temps):
            # Q10 factors
            q10_enz[i] = q10_factor(Q10_ENZYME, temp, CONST.T_phys)
            q10_lk[i] = q10_factor(Q10_LEAK, temp, CONST.T_phys)

            # PMF at this temperature
            pmf = proton_motive_force(delta_psi, delta_ph, temp)

            # ATP synthase rate with temperature scaling
            atp_synthase = ATPSynthaseKinetics(
                h_per_atp=self.base_params.h_per_atp,
                k_rotation=100.0 * q10_enz[i],
                km_atp_site=1e-5,
                km_pi_site=1e-3,
                slip_probability=self.base_params.proton_slip_probability,
            )
            rate = atp_synthase.rotation_rate(
                pmf_kj_per_mol=pmf,
                adp_conc=adp_conc,
                pi_conc=pi_conc,
                atp_adp_ratio=atp_adp_ratio,
            )
            atp_rates[i] = rate

            # Proton leak with temperature scaling
            leak_model = ProtonLeakKinetics(
                l0=1e-9 * q10_lk[i],
                alpha=0.06,
                q10=Q10_LEAK,
                temp_ref=CONST.T_phys,
            )
            leak = leak_model.leak_rate(pmf, temp)
            leak_vals[i] = leak

            # Temperature correction for ΔG
            # ASSUMPTION: The free energy of glucose oxidation changes by
            # approximately -0.2 kJ/(mol*K) with temperature (from Kirchhoff's
            # law). The effect on ATP synthesis is smaller.
            dg_glucose_at_t = STANDARD_GIBBS_GLUCOSE - 0.2 * (temp - CONST.T_std)
            dg_atp_at_t = STANDARD_GIBBS_ATP + 0.02 * (temp - CONST.T_std)

            theoretical_max[i] = compute_theoretical_max_atp(
                dg_glucose_at_t, dg_atp_at_t
            )

            # Estimate ATP yield
            # ASSUMPTION: Yield is limited by the ratio of net synthase rate
            # to reference rate, scaled by the thermodynamic ceiling.
            h_per_atp = self.base_params.h_per_atp
            leak_cost = leak / max(h_per_atp, 1.0)
            net_rate = max(rate - leak_cost, 0.0)

            # Normalise to yield
            ref_rate = 1e-6  # M/s, reference
            yield_ = (
                net_rate / max(ref_rate, 1e-30)
                * min(theoretical_max[i], 94.0)
            )
            yield_ = max(min(yield_, theoretical_max[i]), 0.0)
            atp_yields[i] = yield_
            efficiencies[i] = compute_coupling_efficiency(
                yield_, dg_glucose_at_t, dg_atp_at_t
            )

        # Find optimum
        opt_idx = int(np.argmax(atp_yields))
        optimum_temp = float(temps[opt_idx])
        yield_at_opt = float(atp_yields[opt_idx])

        self.result = TemperatureEffectsResult(
            temperatures=temps,
            atp_yields=atp_yields,
            efficiencies=efficiencies,
            theoretical_max=theoretical_max,
            q10_enzyme_factors=q10_enz,
            q10_leak_factors=q10_lk,
            atp_synthase_rates=atp_rates,
            leak_rates=leak_vals,
            optimum_temp=optimum_temp,
            yield_at_optimum=yield_at_opt,
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
            "  TEMPERATURE EFFECTS EXPERIMENT — PTR-94",
            "=" * 84,
            f"  Temperature range: [{r.temperatures[0]:.1f}, "
            f"{r.temperatures[-1]:.1f}] K  "
            f"({len(r.temperatures)} steps)",
            f"  Reference temperature: {CONST.T_phys:.1f} K "
            f"({CONST.T_phys - 273.15:.1f} C)",
            f"  Q10 (enzyme): {Q10_ENZYME}  |  Q10 (leak): {Q10_LEAK}  |  "
            f"Q10 (diffusion): {Q10_DIFFUSION}",
            "-" * 84,
            f"  {'T(K)':>7} {'T(C)':>6} {'ATP':>8} {'Eff%':>7} "
            f"{'Max':>8} {'Synth':>10} {'Leak':>10}",
            f"  {'':>7} {'':>6} {'yield':>8} {'':>7} "
            f"{'theor':>8} {'(M/s)':>10} {'(M/s)':>10}",
            "-" * 84,
        ]

        n = len(r.temperatures)
        opt_idx = int(np.argmax(r.atp_yields))
        indices = {0, n // 4, n // 2, 3 * n // 4, n - 1, opt_idx}
        for idx in sorted(indices):
            tc = r.temperatures[idx] - 273.15
            lines.append(
                f"  {r.temperatures[idx]:>7.1f} {tc:>6.1f} "
                f"{r.atp_yields[idx]:>8.2f} {r.efficiencies[idx]:>7.2f} "
                f"{r.theoretical_max[idx]:>8.2f} "
                f"{r.atp_synthase_rates[idx]:>10.2e} "
                f"{r.leak_rates[idx]:>10.2e}"
            )

        lines.extend([
            "-" * 84,
            f"  Optimum temperature:         {r.optimum_temp:.1f} K "
            f"({r.optimum_temp - 273.15:.1f} C)",
            f"  ATP yield at optimum:        {r.yield_at_optimum:.2f}",
            f"  Yield at 310.15 K (37 C):    "
            f"{float(np.interp(CONST.T_phys, r.temperatures, r.atp_yields)):.2f}",
            f"  Yield at 298.15 K (25 C):    "
            f"{float(np.interp(CONST.T_std, r.temperatures, r.atp_yields)):.2f}",
            "",
            "  Observations:",
            "  * Below optimum: enzyme rates are too slow.",
            "  * Above optimum: proton leak dominates, reducing yield.",
            "  * The optimal temperature is a trade-off between kinetic"
            " rate and membrane integrity.",
            "=" * 84,
        ])
        return "\n".join(lines)

    def data_table(self) -> Dict[str, List[float]]:
        """Return structured data for export."""
        if self.result is None:
            return {}
        r = self.result
        return {
            "temperature_K": r.temperatures.tolist(),
            "temperature_C": (r.temperatures - 273.15).tolist(),
            "atp_yield": r.atp_yields.tolist(),
            "efficiency_pct": r.efficiencies.tolist(),
            "theoretical_max_atp": r.theoretical_max.tolist(),
            "atp_synthase_rate_M_per_s": r.atp_synthase_rates.tolist(),
            "leak_rate_M_per_s": r.leak_rates.tolist(),
        }


if __name__ == "__main__":
    exp = TemperatureEffectsExperiment()
    result = exp.run()
    print(result.summary)

    print(f"\n  Optimal temperature: {result.optimum_temp:.1f} K "
          f"({result.optimum_temp - 273.15:.1f} C)")
    print(f"  ATP yield at optimum: {result.yield_at_optimum:.2f}")
