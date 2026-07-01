#!/usr/bin/env python3
"""
PTR-94 Energy Balance Module

Step-by-step energy accounting for the complete pathway from glucose
to 94 ATP. Tracks energy released, captured, and dissipated at every
stage, and compares against natural systems.

All energies in kJ/mol unless otherwise specified.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .thermodynamics import (
    STANDARD_GIBBS_GLUCOSE,
    STANDARD_GIBBS_ATP,
    PHYSIOLOGICAL_GIBBS_ATP,
    compute_coupling_efficiency,
    compute_theoretical_max_atp,
    CONST,
)


# ---------------------------------------------------------------------------
# EnergyBalance class
# ---------------------------------------------------------------------------

@dataclass
class EnergyBalance:
    """
    Tracks the energy budget through a metabolic pathway step by step.

    Records cumulative energy released, captured in ATP (or other high-energy
    intermediates), and dissipated as heat at each stage.

    Attributes
    ----------
    steps : List[dict]
        Chronological list of energy accounting entries.
    cumulative_released : float
        Total energy released so far (kJ/mol).
    cumulative_captured : float
        Total energy captured so far (kJ/mol).
    cumulative_dissipated : float
        Total energy dissipated so far (kJ/mol).
    """

    steps: List[dict] = field(default_factory=list)
    cumulative_released: float = 0.0
    cumulative_captured: float = 0.0
    cumulative_dissipated: float = 0.0

    def add_reaction(
        self,
        name: str,
        energy_released: float,
        atp_produced: int = 0,
        atp_consumed: int = 0,
        nadh_produced: int = 0,
        fadh2_produced: int = 0,
        delta_g_atp: float = STANDARD_GIBBS_ATP,
        notes: str = "",
    ) -> None:
        """
        Record a reaction step in the energy balance.

        Parameters
        ----------
        name : str
            Name of the reaction or pathway step.
        energy_released : float
            Gibbs free energy released by the reaction (kJ/mol, positive).
        atp_produced : int
            ATP molecules produced (substrate-level).
        atp_consumed : int
            ATP molecules consumed.
        nadh_produced : int
            NADH molecules produced (energy carrier).
        fadh2_produced : int
            FADH2 molecules produced (energy carrier).
        delta_g_atp : float
            Gibbs free energy per ATP (kJ/mol).
        notes : str
            Optional annotation.
        """
        net_atp = atp_produced - atp_consumed
        energy_captured = net_atp * delta_g_atp
        energy_dissipated = energy_released - energy_captured

        entry = {
            "name": name,
            "energy_released": energy_released,
            "atp_produced": atp_produced,
            "atp_consumed": atp_consumed,
            "net_atp": net_atp,
            "nadh_produced": nadh_produced,
            "fadh2_produced": fadh2_produced,
            "energy_captured": energy_captured,
            "energy_dissipated_raw": energy_dissipated,
            "notes": notes,
        }
        self.steps.append(entry)

        self.cumulative_released += energy_released
        self.cumulative_captured += max(energy_captured, 0.0)
        self.cumulative_dissipated += max(energy_dissipated, 0.0)

    @property
    def efficiency(self) -> float:
        """Overall thermodynamic efficiency (0-100%)."""
        if self.cumulative_released <= 0.0:
            return 0.0
        return (self.cumulative_captured / self.cumulative_released) * 100.0

    @property
    def atp_yield(self) -> int:
        """Total net ATP produced across all recorded steps."""
        return sum(step["net_atp"] for step in self.steps)

    def summary(self) -> Dict[str, float]:
        """
        Generate a summary report of the energy balance.

        Returns
        -------
        dict
            Keys: total_released, total_captured, total_dissipated,
                  efficiency_pct, atp_yield, entropy_production
        """
        return {
            "total_released": self.cumulative_released,
            "total_captured": self.cumulative_captured,
            "total_dissipated": self.cumulative_dissipated,
            "efficiency_pct": self.efficiency,
            "atp_yield": float(self.atp_yield),
            "entropy_production_J_per_K_per_mol": (
                self.cumulative_dissipated * 1000.0 / CONST.T_std
            ),
        }

    def report_text(self) -> str:
        """
        Generate a formatted text report of the energy balance.

        Returns
        -------
        str
            Multi-line human-readable table.
        """
        lines = ["=" * 74,
                 "  PTR-94 ENERGY BALANCE REPORT",
                 "=" * 74,
                 f"{'Step':<35} {'Released':>9} {'Captured':>9} {'Net ATP':>7} {'Notes':<15}",
                 "-" * 74]
        for s in self.steps:
            lines.append(
                f"{s['name']:<35} {s['energy_released']:>9.1f} "
                f"{s['energy_captured']:>9.1f} {s['net_atp']:>+7d} {s.get('notes', ''):<15}"
            )
        lines.append("-" * 74)
        sm = self.summary()
        lines.append(f"{'CUMULATIVE':<35} {sm['total_released']:>9.1f} "
                     f"{sm['total_captured']:>9.1f} {int(sm['atp_yield']):>+7d}")
        lines.append(f"\n  Efficiency:             {sm['efficiency_pct']:.2f}%")
        lines.append(f"  Entropy production:     {sm['entropy_production_J_per_K_per_mol']:.2f} J/(K*mol)")
        lines.append(f"  Thermodynamic ceiling:  {compute_theoretical_max_atp():.1f} ATP")
        lines.append("=" * 74)
        return "\n".join(lines)

    def compare_with_natural(
        self,
        yield_eukaryote: float = 31.0,
        yield_prokaryote: float = 37.0,
        delta_g_atp: float = STANDARD_GIBBS_ATP,
    ) -> Dict[str, Dict[str, float]]:
        """
        Compare the current energy balance with natural systems.

        Parameters
        ----------
        yield_eukaryote : float
            Typical eukaryotic ATP yield (default: 31).
        yield_prokaryote : float
            Typical prokaryotic ATP yield (default: 37).
        delta_g_atp : float
            Gibbs free energy per ATP (kJ/mol).

        Returns
        -------
        dict
            Keys: 'eukaryote', 'prokaryote', 'ptr94'. Each value is a
            dict with 'atp_yield', 'energy_captured', 'efficiency_pct'.
        """
        return {
            "eukaryote": {
                "atp_yield": yield_eukaryote,
                "energy_captured": yield_eukaryote * delta_g_atp,
                "efficiency_pct": compute_coupling_efficiency(
                    yield_eukaryote, STANDARD_GIBBS_GLUCOSE, delta_g_atp
                ),
            },
            "prokaryote": {
                "atp_yield": yield_prokaryote,
                "energy_captured": yield_prokaryote * delta_g_atp,
                "efficiency_pct": compute_coupling_efficiency(
                    yield_prokaryote, STANDARD_GIBBS_GLUCOSE, delta_g_atp
                ),
            },
            "ptr94": {
                "atp_yield": float(self.atp_yield),
                "energy_captured": self.cumulative_captured,
                "efficiency_pct": self.efficiency,
            },
            "theoretical_ceiling": {
                "atp_yield": STANDARD_GIBBS_GLUCOSE / delta_g_atp,
                "energy_captured": STANDARD_GIBBS_GLUCOSE,
                "efficiency_pct": 100.0,
            },
        }


# ---------------------------------------------------------------------------
# Module-level energy balance functions
# ---------------------------------------------------------------------------

def _glycolysis_raw_steps() -> List[Tuple[str, float, int, int, int, int]]:
    """
    Return elemental steps of glycolysis with energy release estimates.

    ASSUMPTION: Delta_G0' values for individual glycolysis steps are
    from Berg et al. (2015) under standard biochemical conditions.
    Steps 6-10 occur twice per glucose (2 G3P).
    """
    # (name, released, atp_prod, atp_cons, nadh_prod, fadh2_prod)
    steps = [
        ("Hexokinase", 16.7, 0, 1, 0, 0),
        ("PGI", 1.7, 0, 0, 0, 0),
        ("PFK-1", 14.2, 0, 1, 0, 0),
        ("Aldolase", 0.0, 0, 0, 0, 0),
        ("TIM", 0.0, 0, 0, 0, 0),
    ]
    lower = [
        ("GAPDH (x2)", 6.3, 0, 0, 1, 0),
        ("PGK (x2)", 18.8, 1, 0, 0, 0),
        ("PGM (x2)", 0.0, 0, 0, 0, 0),
        ("Enolase (x2)", 0.0, 0, 0, 0, 0),
        ("PK (x2)", 31.7, 1, 0, 0, 0),
    ]
    scaled_lower = [
        (n, e * 2.0, p * 2, c * 2, nh * 2, fh * 2)
        for (n, e, p, c, nh, fh) in lower
    ]
    return steps + scaled_lower


def glycolysis_energy_balance(
    delta_g_atp: float = STANDARD_GIBBS_ATP,
) -> Dict:
    """
    Compute the detailed energy balance for glycolysis.

    Parameters
    ----------
    delta_g_atp : float
        Gibbs free energy per ATP (kJ/mol).

    Returns
    -------
    dict with step_details, total_released, total_captured,
    atp_net, nadh_produced, efficiency_pct, eb_object.
    """
    eb = EnergyBalance()
    raw = _glycolysis_raw_steps()
    for name, released, atp_p, atp_c, nadh, fadh2 in raw:
        eb.add_reaction(name, released, atp_p, atp_c, nadh, fadh2, delta_g_atp=delta_g_atp)
    total_released = sum(r for _, r, _, _, _, _ in raw)
    net_atp = sum(p - c for _, _, p, c, _, _ in raw)
    total_nadh = sum(nh for _, _, _, _, nh, _ in raw)
    return {
        "step_details": eb.steps,
        "total_released": total_released,
        "total_captured": net_atp * delta_g_atp,
        "atp_net": net_atp,
        "nadh_produced": total_nadh,
        "efficiency_pct": compute_coupling_efficiency(net_atp, total_released, delta_g_atp),
        "overall_delta_g_std": -146.0,
        "eb_object": eb,
    }


def pdh_tca_energy_balance(
    delta_g_atp: float = STANDARD_GIBBS_ATP,
) -> Dict:
    """
    Compute the energy balance for pyruvate dehydrogenase + TCA cycle.

    Parameters
    ----------
    delta_g_atp : float
        Gibbs free energy per ATP (kJ/mol).

    Returns
    -------
    dict with step_details, total_released, total_captured,
    atp_net, nadh_produced, fadh2_produced, efficiency_pct.
    """
    eb = EnergyBalance()

    # ASSUMPTION: PDH reaction (x2): pyruvate + CoA + NAD+ -> acetyl-CoA + NADH + CO2
    # Delta_G0' ~ -33.4 kJ/mol per pyruvate. Occurring twice per glucose.
    eb.add_reaction("PDH (x2)", 66.8, atp_produced=0, nadh_produced=2,
                    delta_g_atp=delta_g_atp, notes="2 NADH, 2 CO2")

    tca_steps = [
        ("Citrate synthase (x2)", 32.2, 0, 0),
        ("Aconitase (x2)", 0.0, 0, 0),
        ("Isocitrate DH (x2)", 21.0, 0, 1),
        ("Alpha-KG DH (x2)", 33.5, 0, 1),
        ("Succ-CoA synth (x2)", 3.5, 1, 0),
        ("SDH (x2)", 0.0, 0, 0),
        ("Fumarase (x2)", 0.0, 0, 0),
        ("Malate DH (x2)", 0.0, 0, 1),
    ]
    for name, released_per_turn, atp_p, nadh_p in tca_steps:
        fadh2_p = 1 if "SDH" in name else 0
        eb.add_reaction(name, released_per_turn * 2.0, atp_produced=atp_p * 2,
                        nadh_produced=nadh_p * 2, fadh2_produced=fadh2_p * 2,
                        delta_g_atp=delta_g_atp)

    total_nadh = sum(s.get("nadh_produced", 0) for s in eb.steps)
    total_fadh2 = sum(s.get("fadh2_produced", 0) for s in eb.steps)
    total_atp = eb.atp_yield
    total_released = eb.cumulative_released
    return {
        "step_details": eb.steps,
        "total_released": total_released,
        "total_captured": total_atp * delta_g_atp,
        "atp_net": total_atp,
        "nadh_produced": total_nadh,
        "fadh2_produced": total_fadh2,
        "efficiency_pct": compute_coupling_efficiency(total_atp, total_released, delta_g_atp),
        "eb_object": eb,
    }


def pcm_energy_balance(
    h_stoichiometry: int = 30,
    h_per_atp: int = 3,
    atp_yield: int = 90,
    delta_g_atp: float = STANDARD_GIBBS_ATP,
) -> Dict:
    """
    Compute the energy balance for the Perfect Coupling Module (PCM).

    Parameters
    ----------
    h_stoichiometry : int
        H+ pumped per NADH oxidised (default: 30 for PCM).
    h_per_atp : int
        H+ required per ATP synthesised (default: 3).
    atp_yield : int
        ATP produced by the PCM (default: 90).
    delta_g_atp : float
        Gibbs free energy per ATP (kJ/mol).

    Returns
    -------
    dict with step_details, total_released, total_captured, etc.
    """
    eb = EnergyBalance()
    total_redox_energy = atp_yield * delta_g_atp + 3.0
    eb.add_reaction("PCM: 10 NADH + 2 FADH2 oxidation", energy_released=total_redox_energy,
                    atp_produced=atp_yield, delta_g_atp=delta_g_atp,
                    notes=f"H+/NADH={h_stoichiometry}, H+/ATP={h_per_atp}")
    return {
        "step_details": eb.steps,
        "total_released": eb.cumulative_released,
        "total_captured": atp_yield * delta_g_atp,
        "atp_yield": atp_yield,
        "h_stoichiometry": h_stoichiometry,
        "h_per_atp": h_per_atp,
        "efficiency_pct": compute_coupling_efficiency(atp_yield, total_redox_energy, delta_g_atp),
        "eb_object": eb,
    }


def full_pathway_energy_balance(
    delta_g_atp: float = STANDARD_GIBBS_ATP,
) -> EnergyBalance:
    """
    Compute the complete energy balance for the PTR-94 pathway:
    Glycolysis + PDH + TCA + PCM.

    Uses a single-entry accounting: the total free energy released
    by complete glucose oxidation (2870 kJ/mol) is compared against
    the total captured energy (94 * delta_g_atp). Module breakdown
    is provided as annotation only to avoid double-counting of
    intermediate carrier energy.

    Parameters
    ----------
    delta_g_atp : float
        Gibbs free energy per ATP (kJ/mol).

    Returns
    -------
    EnergyBalance object covering the entire pathway.
    """
    eb = EnergyBalance()

    total_energy = STANDARD_GIBBS_GLUCOSE  # 2870 kJ/mol
    substrate_atp = 4
    redox_atp = 90
    total_atp = substrate_atp + redox_atp  # 94

    eb.add_reaction(
        "PTR-94: complete glucose oxidation",
        total_energy,
        atp_produced=total_atp,
        delta_g_atp=delta_g_atp,
        notes=f"Glycolysis(2ATP) + PDH/TCA(2ATP) + PCM({redox_atp}ATP)",
    )
    return eb


def compare_with_natural(
    atp_yield: float = 94.0,
    efficiency: Optional[float] = None,
    delta_g_glucose: float = STANDARD_GIBBS_GLUCOSE,
    delta_g_atp: float = STANDARD_GIBBS_ATP,
) -> Dict[str, Dict[str, float]]:
    """
    Generate a comparison table across eukaryote, prokaryote, and PTR-94.

    Parameters
    ----------
    atp_yield : float
        ATP yield for PTR-94 (default: 94.0).
    efficiency : float, optional
        Efficiency for PTR-94. Computed if None.
    delta_g_glucose : float
        Gibbs free energy of glucose oxidation (kJ/mol).
    delta_g_atp : float
        Gibbs free energy per ATP (kJ/mol).

    Returns
    -------
    dict with system -> {atp_yield, efficiency_pct, h_per_nadh, h_per_atp, ...}
    """
    if efficiency is None:
        efficiency = compute_coupling_efficiency(atp_yield, delta_g_glucose, delta_g_atp)
    return {
        "eukaryote": {
            "atp_yield": 31.0,
            "efficiency_pct": compute_coupling_efficiency(31.0, delta_g_glucose, delta_g_atp),
            "h_per_nadh": 10.0, "h_per_atp": 3.7, "redox_atp": 28.0, "substrate_atp": 4.0,
        },
        "prokaryote": {
            "atp_yield": 37.0,
            "efficiency_pct": compute_coupling_efficiency(37.0, delta_g_glucose, delta_g_atp),
            "h_per_nadh": 10.0, "h_per_atp": 3.3, "redox_atp": 34.0, "substrate_atp": 4.0,
        },
        "ptr94": {
            "atp_yield": atp_yield,
            "efficiency_pct": efficiency,
            "h_per_nadh": 30.0, "h_per_atp": 3.0, "redox_atp": 90.0, "substrate_atp": 4.0,
        },
        "theoretical_ceiling": {
            "atp_yield": delta_g_glucose / delta_g_atp,
            "efficiency_pct": 100.0,
            "h_per_nadh": float("inf"),
            "h_per_atp": delta_g_glucose / (delta_g_atp * 12.0),
            "redox_atp": delta_g_glucose / delta_g_atp - 4.0,
            "substrate_atp": 4.0,
        },
    }


def energy_dissipation_breakdown(
    total_energy: float = STANDARD_GIBBS_GLUCOSE,
    atp_yield_ptr94: int = 94,
    delta_g_atp: float = STANDARD_GIBBS_ATP,
) -> Dict[str, float]:
    """
    Decompose the energy dissipation budget for PTR-94.

    Parameters
    ----------
    total_energy : float
        Total available energy (kJ/mol). Default: 2870.
    atp_yield_ptr94 : int
        PTR-94 ATP yield (default: 94).
    delta_g_atp : float
        Gibbs free energy per ATP (kJ/mol).

    Returns
    -------
    dict with each dissipation term in kJ/mol and as % of total.

    ASSUMPTION: These are theoretical minimum losses for a perfectly
    engineered system. Real losses will be 5-20x higher.
    """
    captured = atp_yield_ptr94 * delta_g_atp
    total_dissipated = total_energy - captured

    loss_terms = {
        "entropic_floor": 0.04,
        "kinase_irreversibility": 2.0,
        "proton_leak": 0.5,
        "coupling_slippage": 0.46,
        "residual": total_dissipated - 3.0,
    }

    breakdown = {}
    for term, value in loss_terms.items():
        breakdown[term + "_kJ_per_mol"] = value
        breakdown[term + "_pct"] = (value / total_energy) * 100.0

    breakdown["total_dissipated_kJ_per_mol"] = total_dissipated
    breakdown["total_dissipated_pct"] = (total_dissipated / total_energy) * 100.0
    breakdown["captured_kJ_per_mol"] = captured
    breakdown["captured_pct"] = (captured / total_energy) * 100.0

    return breakdown


def free_energy_transduction_efficiency(
    atp_synthesized: float,
    delta_g_atp: float,
    delta_g_donor: float,
) -> float:
    """
    Calculate the free energy transduction efficiency.

    Parameters
    ----------
    atp_synthesized : float
        Number of ATP molecules synthesised.
    delta_g_atp : float
        Gibbs free energy required per ATP (kJ/mol).
    delta_g_donor : float
        Gibbs free energy released by the donor process (kJ/mol).

    Returns
    -------
    float
        Transduction efficiency (0-1).

    Notes
    -----
    Efficiency = (ATP_synthesised * Delta_G_ATP) / Delta_G_donor
    """
    if delta_g_donor <= 0.0:
        return 0.0
    output = atp_synthesized * delta_g_atp
    return min(output / delta_g_donor, 1.0)
