#!/usr/bin/env python3
"""
Natural Mitochondria Benchmark.

Benchmarks the PTR-94 pathway against eukaryotic mitochondrial oxidative
phosphorylation. Typical ATP yield: 30-32 per glucose with ~32-34%
thermodynamic efficiency. Provides detailed energy accounting and comparison.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple

import numpy as np

from ..thermodynamics import (
    STANDARD_GIBBS_GLUCOSE,
    STANDARD_GIBBS_ATP,
    PHYSIOLOGICAL_GIBBS_ATP,
    compute_theoretical_max_atp,
    compute_coupling_efficiency,
    entropy_production,
    NATURAL_COMPARISON,
)
from ..energy_balance import EnergyBalance, full_pathway_energy_balance


@dataclass
class NaturalMitochondriaResult:
    """Container for the natural mitochondria benchmark results.

    Attributes:
        atp_yield: ATP per glucose (eukaryotic, ~30-32).
        efficiency: Thermodynamic efficiency (percent, ~32-34).
        h_per_nadh: H+ pumped per NADH oxidised (~10).
        h_per_atp: H+ required per ATP synthesised (~3.7).
        redox_atp: ATP from oxidative phosphorylation (~28).
        substrate_atp: ATP from substrate-level phosphorylation (4).
        entropy: Entropy production (J/(K*mol)).
        energy_balance: Full energy accounting.
        ptr94_comparison: Dict comparing features with PTR-94.
        summary: Formatted text report.
    """
    atp_yield: float = 31.0
    efficiency: float = 32.4
    h_per_nadh: float = 10.0
    h_per_atp: float = 3.7
    redox_atp: float = 28.0
    substrate_atp: int = 4
    entropy: float = 0.0
    energy_balance: EnergyBalance | None = None
    ptr94_comparison: Dict[str, Any] = field(default_factory=dict)
    summary: str = ""


class NaturalMitochondriaBenchmark:
    """Benchmark against eukaryotic mitochondrial respiration.

    Provides detailed energy accounting for the natural system and compares
    it with the PTR-94 design. Eukaryotic yields reflect the cost of
    metabolite shuttles (malate-aspartate, glycerol-3-phosphate) across
    the mitochondrial membrane.

    Attributes:
        delta_g_atp: Gibbs free energy per ATP for the benchmark.
        result: NaturalMitochondriaResult after run().
    """

    def __init__(self, delta_g_atp: float = STANDARD_GIBBS_ATP) -> None:
        """
        Parameters
        ----------
        delta_g_atp : float
            Gibbs free energy per ATP (kJ/mol). Default: standard (30.5).
        """
        self.delta_g_atp = delta_g_atp
        self.result: NaturalMitochondriaResult | None = None

        # ASSUMPTION: Eukaryotic yield of 30-32 ATP/glucose includes:
        #   - 2 ATP (glycolysis, net after hexokinase & PFK)
        #   - 2 ATP (TCA, succinyl-CoA synthetase)
        #   - ~28 ATP from oxidative phosphorylation (10 NADH * 2.5 + 2 FADH2 * 1.5)
        #   - minus 1-2 ATP for shuttle costs (malate-aspartate: 1 ATP-equivalent/NADH)
        #   - minus ~0.5 ATP for ADP/ATP translocase and Pi carrier
        #   Net: 30-32 ATP
        # ASSUMPTION: H+/ATP for natural ATP synthase is ~3.7 (includes
        # the cost of ATP export via the adenine nucleotide translocase).

    def run(self) -> NaturalMitochondriaResult:
        """Execute the benchmark calculation.

        Returns
        -------
        NaturalMitochondriaResult with detailed energy accounting.
        """
        atp_yield = 31.0  # canonical eukaryotic yield
        efficiency = compute_coupling_efficiency(
            atp_yield, STANDARD_GIBBS_GLUCOSE, self.delta_g_atp
        )
        energy_captured = atp_yield * self.delta_g_atp
        entropy = entropy_production(
            STANDARD_GIBBS_GLUCOSE, energy_captured
        )

        # Detailed energy balance
        eb = EnergyBalance()
        eb.add_reaction(
            "Glycolysis (substrate-level)",
            146.0, atp_produced=2, atp_consumed=0,
            delta_g_atp=self.delta_g_atp,
            notes="2 ATP net, 2 NADH",
        )
        eb.add_reaction(
            "PDH + TCA (substrate-level)",
            90.0, atp_produced=2, atp_consumed=0,
            delta_g_atp=self.delta_g_atp,
            notes="2 ATP net, 8 NADH, 2 FADH2",
        )
        # Add mitochondrial shuttles cost
        shuttle_cost = 1.5 * self.delta_g_atp
        eb.add_reaction(
            "Mitochondrial shuttles (malate-aspartate + G3P)",
            1.5 * self.delta_g_atp, atp_produced=0, atp_consumed=0,
            delta_g_atp=self.delta_g_atp,
            notes=f"Cost: {1.5:.1f} ATP equiv.",
        )
        # Natural OXPHOS
        # ASSUMPTION: 10 H+/NADH, 6 H+/FADH2, 3.7 H+/ATP
        h_total = 10.0 * 10 + 6.0 * 2
        redox_atp = h_total / 3.7
        redox_energy = redox_atp * self.delta_g_atp
        eb.add_reaction(
            "Mitochondrial OXPHOS",
            redox_energy + 20.0,  # ~20 kJ dissipated
            atp_produced=int(round(redox_atp)),
            delta_g_atp=self.delta_g_atp,
            notes=f"10 H+/NADH, 6 H+/FADH2, 3.7 H+/ATP",
        )

        # PTR-94 comparison
        ptr94 = NATURAL_COMPARISON.get("ptr94", {})
        ptr94_comparison = {
            "ptr94_atp_yield": ptr94.get("atp_per_glucose", 94.0),
            "ptr94_efficiency": ptr94.get("efficiency_pct", 99.9),
            "ptr94_h_per_nadh": ptr94.get("h_per_nadh", 30.0),
            "ptr94_h_per_atp": ptr94.get("h_per_atp", 3.0),
            "improvement_factor": ptr94.get("atp_per_glucose", 94.0) / atp_yield,
            "efficiency_gain": ptr94.get("efficiency_pct", 99.9) - efficiency,
        }

        self.result = NaturalMitochondriaResult(
            atp_yield=atp_yield,
            efficiency=efficiency,
            h_per_nadh=10.0,
            h_per_atp=3.7,
            redox_atp=redox_atp,
            substrate_atp=4,
            entropy=entropy,
            energy_balance=eb,
            ptr94_comparison=ptr94_comparison,
        )
        self.result.summary = self._generate_summary()
        return self.result

    def _generate_summary(self) -> str:
        """Generate the benchmark report text."""
        if self.result is None:
            return ""
        r = self.result

        lines = [
            "=" * 84,
            "  NATURAL MITOCHONDRIA BENCHMARK — Eukaryotic OXPHOS",
            "=" * 84,
            "",
            "  OVERVIEW",
            "  --------",
            f"  ATP yield:              {r.atp_yield:.1f} per glucose",
            f"  Thermodynamic efficiency: {r.efficiency:.2f}%",
            f"  Entropy production:      {r.entropy:.2f} J/(K*mol)",
            f"  Theoretical ceiling:     {compute_theoretical_max_atp():.1f} ATP",
            "",
            "  ENERGY ACCOUNTING",
            "  -----------------",
            f"  Substrate-level ATP:    {r.substrate_atp} ATP",
            f"  Redox-derived ATP:      {r.redox_atp:.1f} ATP",
            f"  H+/NADH:                {r.h_per_nadh:.0f}",
            f"  H+/FADH2:               6.0",
            f"  H+/ATP synthase:        {r.h_per_atp:.1f}",
            f"  Total H+ pumped:        {10*10 + 6*2:.0f} per glucose",
            "",
            "  SHUTTLE & TRANSPORT COSTS",
            "  -------------------------",
            "  Malate-aspartate shuttle:   1 ATP equiv. per NADH",
            "  Glycerol-3-P shuttle:       0.5 ATP equiv. per NADH",
            "  ADP/ATP translocase:        ~0.3 ATP equiv. per ATP exported",
            "  Pi carrier:                 ~0.1 ATP equiv. per Pi imported",
            "",
            "  PTR-94 COMPARISON",
            "  -----------------",
            f"  PTR-94 ATP yield:          {r.ptr94_comparison['ptr94_atp_yield']:.0f}",
            f"  PTR-94 efficiency:         {r.ptr94_comparison['ptr94_efficiency']:.1f}%",
            f"  Improvement factor:        {r.ptr94_comparison['improvement_factor']:.1f}x",
            f"  Efficiency gain:           {r.ptr94_comparison['efficiency_gain']:.1f} pp",
            "",
            "  LIMITATIONS OF NATURAL SYSTEM",
            "  -----------------------------",
            "  * Proton slip in ATP synthase reduces effective H+/ATP.",
            "  * Membrane leak dissipates ~10-20% of PMF as heat.",
            "  * Metabolite shuttles cost 1-2 ATP per NADH.",
            "  * Superoxide production at Complex I/III wastes ~1-2% of flux.",
            "  * The cristae structure limits diffusion at high rates.",
            "  * H+/ATP ratio of 3.7 vs. theoretical minimum of 3.0.",
            "=" * 84,
        ]
        return "\n".join(lines)

    def report(self) -> str:
        """Generate a detailed benchmark report.

        Returns
        -------
        str
            Multi-line report with energy accounting and comparison.
        """
        if self.result is None:
            self.run()
        if self.result is None or self.result.energy_balance is None:
            return ""
        header = "\nNATURAL MITOCHONDRIA – DETAILED ENERGY ACCOUNTING\n"
        return header + self.result.energy_balance.report_text() + "\n" + self.result.summary

    def data_table(self) -> Dict[str, Any]:
        """Return structured data for export."""
        if self.result is None:
            return {}
        r = self.result
        return {
            "system": "natural_mitochondria",
            "atp_yield": r.atp_yield,
            "efficiency_pct": r.efficiency,
            "h_per_nadh": r.h_per_nadh,
            "h_per_atp": r.h_per_atp,
            "redox_atp": r.redox_atp,
            "substrate_atp": r.substrate_atp,
            "entropy_J_per_K_per_mol": r.entropy,
            "ptr94_comparison": r.ptr94_comparison,
        }


if __name__ == "__main__":
    bench = NaturalMitochondriaBenchmark()
    result = bench.run()
    print(result.summary)
    print("\n" + bench.report())
