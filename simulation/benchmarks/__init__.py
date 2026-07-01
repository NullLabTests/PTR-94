#!/usr/bin/env python3
"""
PTR-94 Benchmark Suite.

Quantitative benchmarks comparing the PTR-94 design against natural and
synthetic energy transduction systems. Each benchmark calculates ATP yield,
thermodynamic efficiency, and energy accounting for a specific metabolic
configuration.
"""

from .natural_mitochondria import *
from .bacterial_respiration import *
from .fermentation import *
from .artificial_pathways import *

__all__ = [
    "NaturalMitochondriaBenchmark",
    "BacterialRespirationBenchmark",
    "FermentationBenchmark",
    "ArtificialPathwaysBenchmark",
    "BENCHMARK_COMPARISON",
]
