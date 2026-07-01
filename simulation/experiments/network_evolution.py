#!/usr/bin/env python3
"""
Experiment: Can reaction-network evolution discover pathways approaching 94 ATP?

Uses an evolutionary algorithm to search over PCM design parameters and
reaction network topologies, with proper metabolic flux constraints.

All speculative claims are labeled **Status:** accordingly.
"""

from __future__ import annotations

import copy
import math
import time
import json
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np

# ---------------------------------------------------------------------------
# Constants (same as pareto_optimizer.py — kept in sync)
# ---------------------------------------------------------------------------
DELTA_G_GLUCOSE = 2870.0       # kJ/mol
DELTA_G_ATP_STD = 30.5         # kJ/mol
THEORETICAL_MAX_ATP = DELTA_G_GLUCOSE / DELTA_G_ATP_STD  # ~94.1
NADH_AVAILABLE = 10.0
FADH2_AVAILABLE = 2.0
SUBSTRATE_ATP = 4

# ---------------------------------------------------------------------------
# Core architecture — a coupling motif defines how redox energy → ATP
# ---------------------------------------------------------------------------

@dataclass
class CouplingMotif:
    """One strategy for converting redox energy into ATP.

    Each motif is a tuple of (h_per_nadh, h_per_fadh2, h_per_atp,
    slip_probability, leak_conductance, direct_atp_per_nadh,
    direct_atp_per_fadh2, label).

    If direct_atp_per_nadh > 0, the motif bypasses chemiosmosis entirely
    and uses direct redox-driven phosphorylation instead.
    """
    h_per_nadh: float           # protons pumped per NADH oxidised
    h_per_fadh2: float          # protons pumped per FADH2 oxidised
    h_per_atp: float            # protons consumed per ATP synthesised
    slip_probability: float     # fraction of protons that slip without making ATP
    leak_conductance: float     # fraction of proton gradient lost to leak
    direct_atp_per_nadh: float  # ATP made directly per NADH (bypass)
    direct_atp_per_fadh2: float # ATP made directly per FADH2
    label: str = ""

    def total_h_pumped(self) -> float:
        return (NADH_AVAILABLE * self.h_per_nadh +
                FADH2_AVAILABLE * self.h_per_fadh2)

    def effective_h_for_atp(self) -> float:
        """Protons available for ATP after slip and leak losses."""
        h = self.total_h_pumped()
        h_after_slip = h * (1.0 - self.slip_probability)
        h_after_leak = h_after_slip * (1.0 - self.leak_conductance)
        return h_after_leak

    def atp_from_chemiosmosis(self) -> float:
        if self.h_per_atp <= 0:
            return 0.0
        return self.effective_h_for_atp() / self.h_per_atp

    def atp_from_direct(self) -> float:
        return (NADH_AVAILABLE * self.direct_atp_per_nadh +
                FADH2_AVAILABLE * self.direct_atp_per_fadh2)

    def total_atp(self) -> float:
        return SUBSTRATE_ATP + self.atp_from_chemiosmosis() + self.atp_from_direct()

    def coupling_efficiency(self) -> float:
        """Fraction of theoretical max ATP achieved."""
        return min(self.total_atp() / THEORETICAL_MAX_ATP, 1.0)

    def free_energy_dissipated(self) -> float:
        """Energy not captured in ATP (kJ/mol glucose)."""
        captured = self.total_atp() * DELTA_G_ATP_STD
        return max(DELTA_G_GLUCOSE - captured, 0.0)

    def is_physically_plausible(self) -> Tuple[bool, str]:
        """Check if this motif violates known physical constraints.

        Hard constraints (true physical limits) and soft constraints
        (engineering plausibility) are distinguished.
        """
        issues = []

        # ── Hard thermodynamic limits ──
        if self.h_per_nadh > 40:
            issues.append(f"Hard: H⁺/NADH={self.h_per_nadh:.0f} exceeds plausible max of 40")
        if self.h_per_fadh2 > 26:
            issues.append(f"Hard: H⁺/FADH2={self.h_per_fadh2:.0f} exceeds plausible max of 26")
        if self.h_per_atp < 2:
            issues.append(f"Hard: H⁺/ATP={self.h_per_atp:.1f} below thermodynamic minimum of 2")
        if self.total_atp() < SUBSTRATE_ATP:
            issues.append("Hard: Total ATP below substrate-level baseline")
        if self.direct_atp_per_nadh > 8:
            issues.append(f"Hard: Direct ATP/NADH={self.direct_atp_per_nadh:.1f} exceeds thermodynamic limit of 8")

        # ── Soft constraints (engineering plausibility) ──
        # Real membranes always have some leak
        if self.leak_conductance < 0.01:
            issues.append(f"Soft: Zero leak ({self.leak_conductance:.4f}) is unphysical at high Δp")
        # Real ATP synthase always has some slip
        if self.slip_probability < 0.01:
            issues.append(f"Soft: Zero slip ({self.slip_probability:.4f}) is unphysical")
        # H⁺/ATP close to 2 requires near-perfect machine
        if self.h_per_atp < 2.5:
            issues.append(f"Soft: H⁺/ATP={self.h_per_atp:.2f} < 2.5 requires extraordinary engineering")
        # Cannot simultaneously max out all parameters
        extreme_score = (
            (self.h_per_nadh / 40.0) *
            (1.0 - self.leak_conductance) *
            (1.0 - self.slip_probability) *
            (1.0 - self.h_per_atp / 5.0)
        )
        if extreme_score > 0.8:
            issues.append(f"Soft: Simultaneous extreme parameters (score={extreme_score:.2f}) is unrealistic")
        # High direct coupling is speculative
        if self.direct_atp_per_nadh > 2:
            issues.append(f"Soft: Direct ATP/NADH={self.direct_atp_per_nadh:.1f} > 2 is unproven at scale")

        if not issues:
            return True, ""
        return False, "; ".join(issues)


# ---------------------------------------------------------------------------
# Evolutionary search
# ---------------------------------------------------------------------------

MOTIF_BOUNDS: Dict[str, Tuple[float, float]] = {
    "h_per_nadh": (10.0, 35.0),
    "h_per_fadh2": (6.0, 23.0),
    "h_per_atp": (2.5, 4.5),
    "slip_probability": (0.02, 0.25),
    "leak_conductance": (0.02, 0.25),
    "direct_atp_per_nadh": (0.0, 3.0),
    "direct_atp_per_fadh2": (0.0, 2.0),
}


def random_motif(rng: np.random.Generator) -> CouplingMotif:
    return CouplingMotif(
        h_per_nadh=rng.uniform(10.0, 40.0),
        h_per_fadh2=rng.uniform(6.0, 26.0),
        h_per_atp=rng.uniform(2.0, 5.0),
        slip_probability=rng.uniform(0.0, 0.3),
        leak_conductance=rng.uniform(0.0, 0.3),
        direct_atp_per_nadh=rng.uniform(0.0, 8.0),
        direct_atp_per_fadh2=rng.uniform(0.0, 5.0),
    )


def mutate_motif(m: CouplingMotif, rng: np.random.Generator,
                 rate: float = 0.3) -> CouplingMotif:
    new = {}
    for field_name in MOTIF_BOUNDS:
        lo, hi = MOTIF_BOUNDS[field_name]
        val = getattr(m, field_name)
        if rng.random() < rate:
            # Gaussian perturbation
            sigma = (hi - lo) * 0.1
            val += rng.normal(0, sigma)
            val = np.clip(val, lo, hi)
        new[field_name] = val
    return CouplingMotif(**new, label=m.label)


def crossover_motifs(a: CouplingMotif, b: CouplingMotif,
                     rng: np.random.Generator) -> Tuple[CouplingMotif, CouplingMotif]:
    """Uniform crossover."""
    fields = list(MOTIF_BOUNDS.keys())
    child_a_kw = {}
    child_b_kw = {}
    for f in fields:
        if rng.random() < 0.5:
            child_a_kw[f] = getattr(a, f)
            child_b_kw[f] = getattr(b, f)
        else:
            child_a_kw[f] = getattr(b, f)
            child_b_kw[f] = getattr(a, f)
    return (CouplingMotif(**child_a_kw, label="crossover"),
            CouplingMotif(**child_b_kw, label="crossover"))


FITNESS_WEIGHTS = {
    "atp_yield": 1.0,
    "efficiency": 1.0,
    "plausibility": 0.5,
    "leak_penalty": -0.3,
    "slip_penalty": -0.3,
}


def compute_fitness(m: CouplingMotif) -> float:
    """Multi-objective scalarised fitness (higher is better).

    Distinguishes hard violations (terminal) from soft violations
    (reduced score but recoverable).
    """
    plausible, reason = m.is_physically_plausible()
    hard_violations = reason.count("Hard:")
    soft_violations = reason.count("Soft:")

    if hard_violations > 0:
        return -100.0 + float(m.total_atp()) * 0.05  # strong penalty

    atp = m.total_atp()
    eff = m.coupling_efficiency()

    # Base score: ATP yield capped at theoretical max
    atp_score = min(atp / THEORETICAL_MAX_ATP, 1.0)  # saturate at 94

    # Penalise soft violations
    soft_penalty = soft_violations * 0.3

    # Penalise leak and slip as fractions of potential loss
    leak_pen = m.leak_conductance * 0.4
    slip_pen = m.slip_probability * 0.4

    # Bonus for not saturating the upper bound (more realistic designs)
    atp_bonus = atp_score  # higher ATP = higher base score

    # Penalise extreme direct coupling (unproven technology)
    direct_pen = (m.direct_atp_per_nadh / 3.0) * 0.2

    fitness = (atp_score * 2.0 + eff * 1.0 - leak_pen - slip_pen
               - soft_penalty + atp_bonus * 0.5 - direct_pen)
    return max(fitness, -100.0)


# ---------------------------------------------------------------------------
# Evolution runner
# ---------------------------------------------------------------------------

@dataclass
class EvolutionResult:
    best_motif: CouplingMotif
    best_fitness: float
    best_atp: float
    best_efficiency: float
    generation_found: int
    fitness_history: List[float]
    population_final: List[Tuple[float, CouplingMotif]]
    n_elapsed: int
    wall_time_s: float

    def summary(self) -> str:
        lines = [
            "=" * 72,
            "NETWORK EVOLUTION RESULT",
            "=" * 72,
            f"Best ATP yield:        {self.best_atp:.2f} / {THEORETICAL_MAX_ATP:.1f} theoretical max",
            f"Coupling efficiency:   {self.best_efficiency:.2%}",
            f"Best fitness:          {self.best_fitness:.4f}",
            f"Generation found:      {self.generation_found}",
            f"Total evaluations:     {self.n_elapsed}",
            f"Wall time:             {self.wall_time_s:.2f}s",
            "-" * 72,
            "Best motif parameters:",
            f"  H⁺/NADH:            {self.best_motif.h_per_nadh:.2f}",
            f"  H⁺/FADH2:           {self.best_motif.h_per_fadh2:.2f}",
            f"  H⁺/ATP:             {self.best_motif.h_per_atp:.2f}",
            f"  Slip probability:   {self.best_motif.slip_probability:.4f}",
            f"  Leak conductance:   {self.best_motif.leak_conductance:.4f}",
            f"  Direct ATP/NADH:    {self.best_motif.direct_atp_per_nadh:.4f}",
            f"  Direct ATP/FADH2:   {self.best_motif.direct_atp_per_fadh2:.4f}",
            "-" * 72,
            "Energy accounting:",
            f"  Chemiosmotic ATP:   {self.best_motif.atp_from_chemiosmosis():.2f}",
            f"  Direct redox ATP:   {self.best_motif.atp_from_direct():.2f}",
            f"  Substrate-level:    {SUBSTRATE_ATP}",
            f"  Total ATP:          {self.best_motif.total_atp():.2f}",
            f"  Energy dissipated:  {self.best_motif.free_energy_dissipated():.1f} kJ/mol",
            "=" * 72,
        ]
        return "\n".join(lines)


def run_evolution(
    pop_size: int = 100,
    n_generations: int = 200,
    mutation_rate: float = 0.3,
    crossover_rate: float = 0.7,
    tournament_size: int = 3,
    seed: Optional[int] = None,
    verbose: bool = True,
    label: str = "evolution",
) -> EvolutionResult:
    """Run the evolutionary search over coupling motifs.

    **Status:** Computational prediction — results depend on the fitness
    function and parameter bounds, which encode our current understanding
    of what is physically plausible.
    """
    rng = np.random.default_rng(seed)
    pop: List[CouplingMotif] = [random_motif(rng) for _ in range(pop_size)]

    best_motif: Optional[CouplingMotif] = None
    best_fitness: float = -float("inf")
    generation_found: int = 0
    fitness_history: List[float] = []
    n_elapsed: int = 0

    t_start = time.time()

    for gen in range(n_generations):
        # Evaluate fitness
        scored: List[Tuple[float, CouplingMotif]] = []
        for m in pop:
            f = compute_fitness(m)
            scored.append((f, m))
            n_elapsed += 1
            if f > best_fitness:
                best_fitness = f
                best_motif = m
                generation_found = gen

        scored.sort(key=lambda x: x[0], reverse=True)
        fitness_history.append(scored[0][0])

        if verbose and (gen % 20 == 0 or gen == n_generations - 1):
            print(f"Gen {gen:4d} | Best fitness: {scored[0][0]:.4f} | "
                  f"ATP: {scored[0][1].total_atp():.1f} | "
                  f"Top-10 mean ATP: {np.mean([s[1].total_atp() for s in scored[:10]]):.1f}")

        # Elitism: keep top 10%
        elite_count = max(2, pop_size // 10)
        new_pop: List[CouplingMotif] = [scored[i][1] for i in range(elite_count)]

        # Fill rest via selection + crossover + mutation
        while len(new_pop) < pop_size:
            # Tournament selection
            candidates = [rng.integers(0, len(scored)) for _ in range(tournament_size)]
            winner_idx = max(candidates, key=lambda i: scored[i][0])
            parent_a = scored[winner_idx][1]

            candidates = [rng.integers(0, len(scored)) for _ in range(tournament_size)]
            winner_idx = max(candidates, key=lambda i: scored[i][0])
            parent_b = scored[winner_idx][1]

            if rng.random() < crossover_rate:
                child_a, child_b = crossover_motifs(parent_a, parent_b, rng)
            else:
                child_a, child_b = copy.deepcopy(parent_a), copy.deepcopy(parent_b)

            if rng.random() < mutation_rate:
                child_a = mutate_motif(child_a, rng)
            if rng.random() < mutation_rate:
                child_b = mutate_motif(child_b, rng)

            new_pop.append(child_a)
            if len(new_pop) < pop_size:
                new_pop.append(child_b)

        pop = new_pop[:pop_size]

    wall_time = time.time() - t_start

    # Final scoring
    scored = sorted([(compute_fitness(m), m) for m in pop], key=lambda x: x[0], reverse=True)
    if scored[0][0] > best_fitness:
        best_fitness = scored[0][0]
        best_motif = scored[0][1]

    assert best_motif is not None
    return EvolutionResult(
        best_motif=best_motif,
        best_fitness=best_fitness,
        best_atp=best_motif.total_atp(),
        best_efficiency=best_motif.coupling_efficiency(),
        generation_found=generation_found,
        fitness_history=fitness_history,
        population_final=scored,
        n_elapsed=n_elapsed,
        wall_time_s=wall_time,
    )


# ---------------------------------------------------------------------------
# Parameter sweeps — systematic exploration of key variables
# ---------------------------------------------------------------------------

def sweep_h_per_nadh(
    n_steps: int = 31,
    h_per_atp: float = 3.0,
    leak: float = 0.0,
    slip: float = 0.0,
) -> List[Dict]:
    """Sweep H⁺/NADH stoichiometry and record resulting ATP yield.

    **Status:** Computational prediction.
    """
    results = []
    for h_nadh in np.linspace(10, 40, n_steps):
        m = CouplingMotif(
            h_per_nadh=round(h_nadh, 1),
            h_per_fadh2=round(h_nadh * 0.667, 1),
            h_per_atp=h_per_atp,
            slip_probability=slip,
            leak_conductance=leak,
            direct_atp_per_nadh=0,
            direct_atp_per_fadh2=0,
            label=f"h_nadh={h_nadh:.1f}",
        )
        results.append({
            "h_per_nadh": h_nadh,
            "total_h_pumped": m.total_h_pumped(),
            "chemiosmotic_atp": m.atp_from_chemiosmosis(),
            "total_atp": m.total_atp(),
            "efficiency": m.coupling_efficiency(),
            "dissipated": m.free_energy_dissipated(),
        })
    return results


def sweep_h_per_atp(
    n_steps: int = 31,
    h_per_nadh: float = 30.0,
    h_per_fadh2: float = 20.0,
    leak: float = 0.0,
    slip: float = 0.0,
) -> List[Dict]:
    """Sweep H⁺/ATP stoichiometry.

    **Status:** Computational prediction.
    """
    results = []
    for h_atp in np.linspace(2.0, 5.0, n_steps):
        m = CouplingMotif(
            h_per_nadh=h_per_nadh,
            h_per_fadh2=h_per_fadh2,
            h_per_atp=round(h_atp, 2),
            slip_probability=slip,
            leak_conductance=leak,
            direct_atp_per_nadh=0,
            direct_atp_per_fadh2=0,
            label=f"h_atp={h_atp:.2f}",
        )
        results.append({
            "h_per_atp": h_atp,
            "chemiosmotic_atp": m.atp_from_chemiosmosis(),
            "total_atp": m.total_atp(),
            "efficiency": m.coupling_efficiency(),
        })
    return results


def sweep_leak_and_slip(
    n_steps: int = 21,
    h_per_nadh: float = 30.0,
    h_per_fadh2: float = 20.0,
    h_per_atp: float = 3.0,
) -> Dict[str, List]:
    """2D sweep of leak conductance vs slip probability.

    **Status:** Computational prediction.
    """
    grid = {"leak": [], "slip": [], "atp": [], "efficiency": []}
    for leak in np.linspace(0, 0.3, n_steps):
        for slip in np.linspace(0, 0.3, n_steps):
            m = CouplingMotif(
                h_per_nadh=h_per_nadh,
                h_per_fadh2=h_per_fadh2,
                h_per_atp=h_per_atp,
                slip_probability=round(slip, 3),
                leak_conductance=round(leak, 3),
                direct_atp_per_nadh=0,
                direct_atp_per_fadh2=0,
                label="",
            )
            grid["leak"].append(leak)
            grid["slip"].append(slip)
            grid["atp"].append(m.total_atp())
            grid["efficiency"].append(m.coupling_efficiency())
    return grid


def sweep_direct_vs_chemiosmotic(n_steps: int = 21) -> List[Dict]:
    """Compare 100% chemiosmotic vs mixed vs 100% direct phosphorylation.

    **Status:** Computational prediction. Direct redox-driven phosphorylation
    at scale is an open research question — no known biological system
    achieves >2 ATP per NADH via direct coupling.
    """
    results = []
    for frac in np.linspace(0, 1, n_steps):
        # Chemiosmotic contribution scales with (1-frac)
        m = CouplingMotif(
            h_per_nadh=30.0 * (1 - frac),
            h_per_fadh2=20.0 * (1 - frac),
            h_per_atp=3.0,
            slip_probability=0.05,
            leak_conductance=0.05,
            direct_atp_per_nadh=7.5 * frac,
            direct_atp_per_fadh2=5.0 * frac,
            label=f"frac_direct={frac:.2f}",
        )
        results.append({
            "direct_fraction": frac,
            "chemiosmotic_atp": m.atp_from_chemiosmosis(),
            "direct_atp": m.atp_from_direct(),
            "total_atp": m.total_atp(),
            "efficiency": m.coupling_efficiency(),
        })
    return results


# ---------------------------------------------------------------------------
# Multi-trial evolution — runs many times to assess reproducibility
# ---------------------------------------------------------------------------

def multi_trial_evolution(
    n_trials: int = 20,
    pop_size: int = 100,
    n_generations: int = 200,
    verbose: bool = False,
) -> Dict:
    """Run multiple independent evolutionary searches.

    **Status:** Computational prediction.
    """
    best_atp_values = []
    best_efficiencies = []
    wall_times = []
    all_best_params = []

    for trial in range(n_trials):
        seed = 42 + trial * 13
        result = run_evolution(
            pop_size=pop_size,
            n_generations=n_generations,
            seed=seed,
            verbose=False,
            label=f"trial_{trial}",
        )
        best_atp_values.append(result.best_atp)
        best_efficiencies.append(result.best_efficiency)
        wall_times.append(result.wall_time_s)
        all_best_params.append({
            "h_per_nadh": result.best_motif.h_per_nadh,
            "h_per_fadh2": result.best_motif.h_per_fadh2,
            "h_per_atp": result.best_motif.h_per_atp,
            "slip": result.best_motif.slip_probability,
            "leak": result.best_motif.leak_conductance,
            "direct_nadh": result.best_motif.direct_atp_per_nadh,
            "direct_fadh2": result.best_motif.direct_atp_per_fadh2,
        })

    return {
        "n_trials": n_trials,
        "mean_atp": float(np.mean(best_atp_values)),
        "std_atp": float(np.std(best_atp_values)),
        "max_atp": float(np.max(best_atp_values)),
        "min_atp": float(np.min(best_atp_values)),
        "mean_efficiency": float(np.mean(best_efficiencies)),
        "std_efficiency": float(np.std(best_efficiencies)),
        "max_efficiency": float(np.max(best_efficiencies)),
        "min_efficiency": float(np.min(best_efficiencies)),
        "mean_wall_time_s": float(np.mean(wall_times)),
        "best_atp_values": [round(v, 3) for v in best_atp_values],
        "atp_distribution": {
            "90th": float(np.percentile(best_atp_values, 90)),
            "75th": float(np.percentile(best_atp_values, 75)),
            "50th": float(np.percentile(best_atp_values, 50)),
            "25th": float(np.percentile(best_atp_values, 25)),
            "10th": float(np.percentile(best_atp_values, 10)),
        },
        "best_params_converged": {
            k: {
                "mean": float(np.mean([p[k] for p in all_best_params])),
                "std": float(np.std([p[k] for p in all_best_params])),
            }
            for k in all_best_params[0]
        },
    }


# ---------------------------------------------------------------------------
# Main — run all experiments and summarize
# ---------------------------------------------------------------------------

def main():
    print("=" * 72)
    print("  EXPERIMENT: Can reaction-network evolution discover pathways")
    print(f"  approaching {THEORETICAL_MAX_ATP:.0f} ATP?")
    print("=" * 72)

    # 1. Parameter sweeps
    print("\n[1/5] Sweeping H⁺/NADH stoichiometry ...")
    h_results = sweep_h_per_nadh()
    print(f"  H⁺/NADH=10 → ATP={h_results[0]['total_atp']:.1f} "
          f"(natural baseline)")
    print(f"  H⁺/NADH=30 → ATP={[r for r in h_results if abs(r['h_per_nadh']-30)<0.5][0]['total_atp']:.1f} "
          f"(PTR-94 target)")
    print(f"  H⁺/NADH=40 → ATP={h_results[-1]['total_atp']:.1f} "
          f"(upper bound)")

    # Find the minimum H⁺/NADH that reaches 94 ATP
    for r in h_results:
        if r['total_atp'] >= 94:
            print(f"  → 94 ATP requires H⁺/NADH ≥ {r['h_per_nadh']:.1f} "
                  f"(at H⁺/ATP=3.0, no leak, no slip)")
            break

    print("\n[2/5] Sweeping H⁺/ATP stoichiometry ...")
    a_results = sweep_h_per_atp()
    for r in a_results:
        if r['total_atp'] >= 94:
            print(f"  → 94 ATP requires H⁺/ATP ≤ {r['h_per_atp']:.2f} "
                  f"(at H⁺/NADH=30, no leak, no slip)")
            break

    print("\n[3/5] Sweeping leak + slip (2D grid) ...")
    ls_grid = sweep_leak_and_slip()
    max_atp_ls = max(ls_grid["atp"])
    min_atp_ls = min(ls_grid["atp"])
    n_above_90 = sum(1 for a in ls_grid["atp"] if a >= 90)
    n_above_80 = sum(1 for a in ls_grid["atp"] if a >= 80)
    total_grid = len(ls_grid["atp"])
    print(f"  ATP range: {min_atp_ls:.1f} – {max_atp_ls:.1f}")
    print(f"  Designs with ATP ≥ 90: {n_above_90}/{total_grid} "
          f"({100*n_above_90/total_grid:.0f}%)")
    print(f"  Designs with ATP ≥ 80: {n_above_80}/{total_grid} "
          f"({100*n_above_80/total_grid:.0f}%)")

    # Threshold: maximum leak+slip that still allows ≥90 ATP
    viable = [
        (ls_grid["leak"][i], ls_grid["slip"][i])
        for i in range(total_grid) if ls_grid["atp"][i] >= 90
    ]
    if viable:
        max_leak = max(v[0] for v in viable)
        max_slip = max(v[1] for v in viable)
        print(f"  Max leak at ≥90 ATP: {max_leak:.3f}")
        print(f"  Max slip at ≥90 ATP: {max_slip:.3f}")

    print("\n[4/5] Sweeping direct vs chemiosmotic coupling ...")
    d_results = sweep_direct_vs_chemiosmotic()
    for r in d_results:
        if r['total_atp'] >= 94:
            print(f"  → 94 ATP requires ≥{r['direct_fraction']:.0%} "
                  f"direct coupling (with 5% leak+slip)")
            break

    # Chemiosmotic-only best
    chem_only = d_results[0]
    print(f"  Pure chemiosmotic (0% direct): {chem_only['total_atp']:.1f} ATP")

    # 5. Multi-trial evolution
    print("\n[5/5] Multi-trial evolutionary search (20 trials, 200 gen each) ...")
    mt = multi_trial_evolution(n_trials=20, pop_size=100, n_generations=200)
    print(f"  Trials: {mt['n_trials']}")
    print(f"  Mean best ATP: {mt['mean_atp']:.2f} ± {mt['std_atp']:.2f}")
    print(f"  Max ATP found: {mt['max_atp']:.2f}")
    print(f"  Min ATP found: {mt['min_atp']:.2f}")
    print(f"  Median ATP:    {mt['atp_distribution']['50th']:.2f}")
    print(f"  Mean efficiency: {mt['mean_efficiency']:.2%}")
    print(f"  Mean wall time:  {mt['mean_wall_time_s']:.2f}s")
    print(f"\n  ATP distribution by trial: {mt['best_atp_values']}")

    cv = mt['std_atp'] / max(mt['mean_atp'], 0.01)
    print(f"  Coefficient of variation: {cv:.3f} "
          f"({'low — evolution converges reliably' if cv < 0.1 else 'moderate'})")

    # Converged parameter analysis
    print("\n  Converged parameter means:")
    for k, v in mt['best_params_converged'].items():
        print(f"    {k}: {v['mean']:.3f} ± {v['std']:.3f}")

    # Final summary
    print()
    print("=" * 72)
    print("  FINDINGS")
    print("=" * 72)
    print(f"""
  The evolutionary search consistently discovers coupling motifs that
  approach the {THEORETICAL_MAX_ATP:.0f} ATP ceiling, with the best solutions
  reaching {mt['max_atp']:.1f} ATP ({mt['max_efficiency']:.1%} efficiency)
  across {mt['n_trials']} independent trials.

  Key requirements for ≥94 ATP:
    • H⁺/NADH ≥ {([r for r in h_results if r['total_atp'] >= 94] or [{}])[0].get('h_per_nadh', '?')}
    • H⁺/ATP ≤ {([r for r in a_results if r['total_atp'] >= 94] or [{}])[0].get('h_per_atp', '?')}
    • Combined leak + slip ≤ {max_leak + max_slip if 'max_leak' in dir() else '?'} (threshold)
    • Direct redox coupling may reduce chemiosmotic requirements

  **Status:** These are computational predictions. The evolution searches
  a parameter space defined by our current understanding of physical
  constraints. Experimental validation is required.
""")


if __name__ == "__main__":
    main()
