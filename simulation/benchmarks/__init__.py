#!/usr/bin/env python3
"""
PTR-94 Benchmark Suite.

Quantitative benchmarks comparing the PTR-94 design against natural and
synthetic energy transduction systems. Each benchmark calculates ATP yield,
thermodynamic efficiency, and energy accounting for a specific metabolic
configuration.
"""

from .natural_mitochondria import NaturalMitochondriaBenchmark
from .bacterial_respiration import BacterialRespirationBenchmark
from .fermentation import FermentationBenchmark
from .artificial_pathways import ArtificialPathwaysBenchmark, BENCHMARK_COMPARISON

__all__ = [
    "NaturalMitochondriaBenchmark",
    "BacterialRespirationBenchmark",
    "FermentationBenchmark",
    "ArtificialPathwaysBenchmark",
    "BENCHMARK_COMPARISON",
]
