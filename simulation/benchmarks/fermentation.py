#!/usr/bin/env python3
"""
Fermentation Benchmark.

Benchmarks the PTR-94 pathway against anaerobic fermentation (glycolysis
to lactate or ethanol). Typical ATP yield: 2 per glucose with only ~2-3%
thermodynamic efficiency. Shows the dramatic improvement from anaerobic
to aerobic to the PTR-94 design.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict


from ..thermodynamics import (
    STANDARD_GIBBS_GLUCOSE,
    STANDARD_GIBBS_ATP,
    compute_coupling_efficiency,
    entropy_production,
)
from ..energy_balance import EnergyBalance


@dataclass
class FermentationResult:
    """Container for the fermentation benchmark results.

    Attributes:
        atp_yield: ATP per glucose (2 for homolactic).
        efficiency: Thermodynamic efficiency (percent, ~2-3).
        atp_from_glycolysis: ATP from glycolysis alone (2).
        atp_from_tca: ATP from TCA (0, anaerobic).
        atp_from_oxphos: ATP from OXPHOS (0, anaerobic).
        end_product: Fermentation end product ('lactate' or 'ethanol').
        energy_balance: Full energy accounting.
        anaerobic_gap: Energy difference from glucose to 2 ATP (kJ).
        aerobic_yield: Typical aerobic yield for comparison.
        ptr94_yield: PTR-94 yield for comparison.
        summary: Formatted text report.
    """
    atp_yield: float = 2.0
    efficiency: float = 2.1
    atp_from_glycolysis: int = 2
    atp_from_tca: int = 0
    atp_from_oxphos: int = 0
    end_product: str = "lactate"
    entropy: float = 0.0
    energy_balance: EnergyBalance | None = None
    anaerobic_gap: float = 0.0
    aerobic_yield: float = 31.0
    ptr94_yield: float = 94.0
    summary: str = ""


class FermentationBenchmark:
    """Benchmark against anaerobic fermentation.

    In the absence of oxygen, cells rely on glycolysis alone for ATP
    production, fermenting pyruvate to lactate or ethanol to regenerate
    NAD+. This captures only ~2% of the available energy in glucose.

    Attributes:
        delta_g_atp: Gibbs free energy per ATP for the benchmark.
        end_product: 'lactate' or 'ethanol'.
        result: FermentationResult after run().
    """

    def __init__(
        self,
        delta_g_atp: float = STANDARD_GIBBS_ATP,
        end_product: str = "lactate",
    ) -> None:
        """
        Parameters
        ----------
        delta_g_atp : float
            Gibbs free energy per ATP (kJ/mol).
        end_product : str
            'lactate' or 'ethanol'. Affects total energy released.
        """
        self.delta_g_atp = delta_g_atp
        self.end_product = end_product
        self.result: FermentationResult | None = None

        # ASSUMPTION: Homolactic fermentation (glucose -> 2 lactate):
        #   Glucose + 2 ADP + 2 Pi -> 2 Lactate + 2 ATP + 2 H2O
        #   ΔG0' = -196 kJ/mol (much less than complete oxidation)
        # ASSUMPTION: Alcoholic fermentation (glucose -> 2 ethanol + 2 CO2):
        #   ΔG0' = -235 kJ/mol

    def run(self) -> FermentationResult:
        """Execute the benchmark calculation.

        Returns
        -------
        FermentationResult with energy accounting.
        """
        atp_yield = 2.0

        # ASSUMPTION: The free energy released by fermentation depends
        # on the end product.
        if self.end_product == "lactate":
            total_energy = 196.0  # kJ/mol for glucose -> 2 lactate
        elif self.end_product == "ethanol":
            total_energy = 235.0  # kJ/mol for glucose -> 2 ethanol + 2 CO2
        else:
            raise ValueError(f"Unknown end product: {self.end_product}")

        efficiency = compute_coupling_efficiency(
            atp_yield, total_energy, self.delta_g_atp
        )
        energy_captured = atp_yield * self.delta_g_atp
        entropy = entropy_production(total_energy, energy_captured)

        # Energy balance
        eb = EnergyBalance()
        eb.add_reaction(
            "Glycolysis (substrate-level)",
            total_energy, atp_produced=2,
            delta_g_atp=self.delta_g_atp,
            notes=f"2 ATP, 2 NADH -> 2 {self.end_product}",
        )

        # Energy not captured
        self.anaerobic_gap = total_energy - energy_captured

        # Comparison to more advanced systems
        aerobic_yield = 31.0
        ptr94_yield = 94.0

        self.result = FermentationResult(
            atp_yield=atp_yield,
            efficiency=efficiency,
            entropy=entropy,
            atp_from_glycolysis=2,
            end_product=self.end_product,
            energy_balance=eb,
            anaerobic_gap=self.anaerobic_gap,
            aerobic_yield=aerobic_yield,
            ptr94_yield=ptr94_yield,
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
            f"  FERMENTATION BENCHMARK — Anaerobic ({self.end_product})",
            "=" * 84,
            "",
            "  OVERVIEW",
            "  --------",
            f"  ATP yield:              {r.atp_yield:.0f} per glucose",
            f"  Thermodynamic efficiency: {r.efficiency:.2f}%",
            f"  Entropy production:      {r.entropy:.2f} J/(K*mol)",
            f"  Energy not captured:     {r.anaerobic_gap:.0f} kJ/mol "
            f"({r.anaerobic_gap / STANDARD_GIBBS_GLUCOSE * 100:.1f}% of total)",
            "",
            "  WHY SO INEFFICIENT?",
            "  -------------------",
            "  * Only 2/4 substrate-level phosphorylation steps used.",
            "  * NADH cannot be recycled without oxygen.",
            "  * Most glucose carbon is excreted as organic acids/alcohol.",
            "  * The electron transport chain is completely inactive.",
            f"  * Only ~{r.anaerobic_gap:.0f} of {STANDARD_GIBBS_GLUCOSE:.0f} kJ/mol is captured.",
            "",
            "  EVOLUTIONARY PROGRESSION",
            "  ------------------------",
            f"  {'System':<35} {'ATP/glucose':>12} {'Efficiency':>12}",
            "-" * 60,
            f"  {'Anaerobic fermentation':<35} {r.atp_yield:>12.0f} {r.efficiency:>11.2f}%",
            f"  {'Aerobic (eukaryote)':<35} {r.aerobic_yield:>12.0f} "
            f"{compute_coupling_efficiency(r.aerobic_yield, STANDARD_GIBBS_GLUCOSE, self.delta_g_atp):>11.2f}%",
            f"  {'Aerobic (prokaryote)':<35} 37.0 "
            f"{compute_coupling_efficiency(37.0, STANDARD_GIBBS_GLUCOSE, self.delta_g_atp):>11.2f}%",
            f"  {'PTR-94 (theoretical max)':<35} {r.ptr94_yield:>12.0f} "
            f"{compute_coupling_efficiency(r.ptr94_yield, STANDARD_GIBBS_GLUCOSE, self.delta_g_atp):>11.2f}%",
            "-" * 60,
            "",
            "  IMPROVEMENT RATIOS",
            "  ------------------",
            f"  PTR-94 / Fermentation:        {r.ptr94_yield / r.atp_yield:.1f}x",
            f"  Eukaryote / Fermentation:     {r.aerobic_yield / r.atp_yield:.1f}x",
            f"  PTR-94 / Eukaryote:           {r.ptr94_yield / r.aerobic_yield:.1f}x",
            "",
            "  KEY INSIGHT",
            "  -----------",
            "  Fermentation captures only the substrate-level ATP from",
            "  glycolysis. The complete PTR-94 pathway captures redox",
            "  energy through the electron transport chain, increasing",
            "  yield by 47x over anaerobic metabolism.",
            "=" * 84,
        ]
        return "\n".join(lines)

    def report(self) -> str:
        """Generate a detailed benchmark report."""
        if self.result is None:
            self.run()
        if self.result is None or self.result.energy_balance is None:
            return ""
        header = f"\nFERMENTATION ({self.end_product}) – ENERGY ACCOUNTING\n"
        return header + self.result.energy_balance.report_text() + "\n" + self.result.summary

    def data_table(self) -> Dict[str, Any]:
        """Return structured data for export."""
        if self.result is None:
            return {}
        r = self.result
        return {
            "system": f"fermentation_{r.end_product}",
            "atp_yield": r.atp_yield,
            "efficiency_pct": r.efficiency,
            "substrate_atp": r.atp_from_glycolysis,
            "tca_atp": r.atp_from_tca,
            "oxphos_atp": r.atp_from_oxphos,
            "anaerobic_gap_kJ": r.anaerobic_gap,
            "aerobic_atp_yield": r.aerobic_yield,
            "ptr94_atp_yield": r.ptr94_yield,
            "improvement_over_fermentation": r.ptr94_yield / r.atp_yield,
        }


if __name__ == "__main__":
    bench = FermentationBenchmark()
    result = bench.run()
    print(result.summary)
