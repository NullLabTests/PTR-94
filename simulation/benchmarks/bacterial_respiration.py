#!/usr/bin/env python3
"""
Bacterial Respiration Benchmark.

Benchmarks the PTR-94 pathway against prokaryotic aerobic respiration.
Typical ATP yield: 36-38 per glucose with ~37-40% thermodynamic efficiency.
Prokaryotes lack mitochondrial compartments, avoiding shuttle costs.
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
from ..energy_balance import EnergyBalance


@dataclass
class BacterialRespirationResult:
    """Container for the bacterial respiration benchmark results.

    Attributes:
        atp_yield: ATP per glucose (~36-38).
        efficiency: Thermodynamic efficiency (percent, ~37-40).
        h_per_nadh: H+ pumped per NADH oxidised (~8-10).
        h_per_atp: H+ required per ATP synthesised (~3.3).
        redox_atp: ATP from oxidative phosphorylation (~34).
        substrate_atp: ATP from substrate-level phosphorylation (4).
        entropy: Entropy production (J/(K*mol)).
        energy_balance: Full energy accounting.
        ptr94_comparison: Comparison dict with PTR-94.
        summary: Formatted text report.
    """
    atp_yield: float = 37.0
    efficiency: float = 39.4
    h_per_nadh: float = 10.0
    h_per_atp: float = 3.3
    redox_atp: float = 34.0
    substrate_atp: int = 4
    entropy: float = 0.0
    energy_balance: EnergyBalance | None = None
    ptr94_comparison: Dict[str, Any] = field(default_factory=dict)
    summary: str = ""


class BacterialRespirationBenchmark:
    """Benchmark against prokaryotic aerobic respiration.

    Bacteria such as E. coli perform aerobic respiration without
    mitochondrial compartmentalisation. This avoids shuttle costs but
    typically has a lower H+/ATP ratio due to a different ATP synthase
    c-ring stoichiometry.

    Attributes:
        delta_g_atp: Gibbs free energy per ATP for the benchmark.
        result: BacterialRespirationResult after run().
    """

    def __init__(self, delta_g_atp: float = STANDARD_GIBBS_ATP) -> None:
        """
        Parameters
        ----------
        delta_g_atp : float
            Gibbs free energy per ATP (kJ/mol). Default: standard (30.5).
        """
        self.delta_g_atp = delta_g_atp
        self.result: BacterialRespirationResult | None = None

        # ASSUMPTION: Prokaryotic yield of 36-38 ATP/glucose includes:
        #   - 2 ATP (glycolysis, net after HK & PFK)
        #   - 2 ATP (TCA, succinyl-CoA synthetase)
        #   - ~32-34 ATP from oxidative phosphorylation
        #   - No shuttle costs (direct oxidation in cytoplasm/periplasm)
        #   - H+/ATP = 3.3 (E. coli has c10-12 ring vs. mammalian c8)
        # ASSUMPTION: Some bacteria (e.g., Paracoccus denitrificans)
        # achieve higher yields due to more efficient electron transport
        # chain organisation.

    def run(self) -> BacterialRespirationResult:
        """Execute the benchmark calculation.

        Returns
        -------
        BacterialRespirationResult with detailed energy accounting.
        """
        atp_yield = 37.0  # canonical prokaryotic yield
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
            146.0, atp_produced=2,
            delta_g_atp=self.delta_g_atp,
            notes="2 ATP net, 2 NADH",
        )
        eb.add_reaction(
            "PDH + TCA (substrate-level)",
            90.0, atp_produced=2,
            delta_g_atp=self.delta_g_atp,
            notes="2 ATP net, 8 NADH, 2 FADH2",
        )

        # Bacterial OXPHOS
        # ASSUMPTION: E. coli uses NADH dehydrogenase I (NDH-1) which
        # pumps 4 H+/NADH, plus cytochrome bo3 which pumps 2 H+/e-.
        # Some NADH may use NDH-2 (non-pumping). Net: ~8-10 H+/NADH.
        # H+/ATP = 3.3 (c10 ring).
        h_total = 8.0 * 10 + 4.0 * 2
        redox_atp = h_total / 3.3
        eb.add_reaction(
            "Bacterial OXPHOS (no shuttles)",
            redox_atp * self.delta_g_atp + 15.0,
            atp_produced=int(round(redox_atp)),
            delta_g_atp=self.delta_g_atp,
            notes="8-10 H+/NADH, 3.3 H+/ATP, no shuttle cost",
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

        self.result = BacterialRespirationResult(
            atp_yield=atp_yield,
            efficiency=efficiency,
            h_per_nadh=10.0,
            h_per_atp=3.3,
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

        efficiency_natural = compute_coupling_efficiency(
            self.result.atp_yield, STANDARD_GIBBS_GLUCOSE, self.delta_g_atp
        )

        lines = [
            "=" * 84,
            "  BACTERIAL RESPIRATION BENCHMARK — Prokaryotic OXPHOS",
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
            f"  H+/NADH (net):          {r.h_per_nadh:.0f}",
            f"  H+/FADH2 (net):         4.0",
            f"  H+/ATP synthase:        {r.h_per_atp:.1f}",
            f"  Total H+ pumped:        {int(8*10 + 4*2)} per glucose (approx.)",
            "",
            "  KEY ADVANTAGES OVER EUKARYOTES",
            "  ------------------------------",
            "  * No mitochondrial shuttles (saves 1-2 ATP/NADH).",
            "  * Direct access to periplasmic proton gradient.",
            "  * Often more efficient ATP synthase (c10 vs c8 ring).",
            "  * No metabolite transport costs across membranes.",
            "",
            "  PTR-94 COMPARISON",
            "  -----------------",
            f"  PTR-94 ATP yield:          {r.ptr94_comparison['ptr94_atp_yield']:.0f}",
            f"  PTR-94 efficiency:         {r.ptr94_comparison['ptr94_efficiency']:.1f}%",
            f"  Improvement over bacteria:  {r.ptr94_comparison['improvement_factor']:.1f}x",
            f"  Efficiency gain:           {r.ptr94_comparison['efficiency_gain']:.1f} pp",
            "",
            "  LIMITATIONS OF NATURAL BACTERIAL SYSTEM",
            "  ----------------------------------------",
            "  * Some NADH uses non-pumping NDH-2, wasting energy.",
            "  * Branched electron transport chains分流 flux.",
            "  * Membrane leak is higher in single-membrane systems.",
            "  * Lower H+/NADH ratio (8-10) vs. theoretical potential.",
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
        header = "\nBACTERIAL RESPIRATION – DETAILED ENERGY ACCOUNTING\n"
        return header + self.result.energy_balance.report_text() + "\n" + self.result.summary

    def data_table(self) -> Dict[str, Any]:
        """Return structured data for export."""
        if self.result is None:
            return {}
        r = self.result
        return {
            "system": "bacterial_respiration",
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
    bench = BacterialRespirationBenchmark()
    result = bench.run()
    print(result.summary)
