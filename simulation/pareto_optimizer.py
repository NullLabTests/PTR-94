#!/usr/bin/env python3
"""
Pareto Optimization Engine for PTR-94 Perfect Coupling Module (PCM).

Implements multi-objective optimization using an NSGA-II-inspired genetic
algorithm to explore the PCM design space and identify trade-offs between ATP
yield, thermodynamic efficiency, carbon efficiency, stability, and entropy
production.

All objectives are formulated as minimization problems for the optimizer.
"""

from __future__ import annotations

import copy
import io
import json
import math
import random
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

# ASSUMPTION: Standard biochemical conditions (ΔG°' = 30.5 kJ/mol for ATP).
# ASSUMPTION: Glucose oxidation releases 2870 kJ/mol.
# ASSUMPTION: 10 NADH + 2 FADH2 from upstream modules (glycolysis + PDH + TCA).
# ASSUMPTION: Substrate-level phosphorylation yields 4 ATP before PCM.

DELTA_G_GLUCOSE: float = 2870.0  # kJ/mol
DELTA_G_ATP_STD: float = 30.5  # kJ/mol
THEORETICAL_MAX_ATP: float = DELTA_G_GLUCOSE / DELTA_G_ATP_STD  # ~94.1
SUBSTRATE_ATP: int = 4
NADH_COUNT: int = 10
FADH2_COUNT: int = 2


# ─────────────────────────────────────────────────────────────────────────────
# Design Parameters
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class DesignParameters:
    """Tunable parameters of the Perfect Coupling Module (PCM).

    Attributes:
        h_per_nadh: Protons pumped per NADH oxidised (range 10–40).
        h_per_fadh2: Protons pumped per FADH2 oxidised (range 6–26).
        h_per_atp: Protons required per ATP synthesised (range 2–5).
        atp_synthase_efficiency: Catalytic efficiency of ATP synthase (0.7–1.0).
        membrane_leak_conductance: Fraction of proton motive force lost to
            membrane leakage (0–1).
        proton_slip_probability: Probability of proton translocation slipping
            without ATP generation (0–1).
        roq_quinone_coupling: Quinone coupling efficiency in the redox chain
            (0–1).
        scaffold_channeling_efficiency: Efficiency of substrate channeling
            through the PCM metabolon (0–1).
        ros_bypass_fraction: Fraction of electron flux diverted to ROS
            mitigation pathways (0–1).
    """

    h_per_nadh: float = 30.0
    h_per_fadh2: float = 20.0
    h_per_atp: float = 3.0
    atp_synthase_efficiency: float = 0.95
    membrane_leak_conductance: float = 0.05
    proton_slip_probability: float = 0.02
    roq_quinone_coupling: float = 0.98
    scaffold_channeling_efficiency: float = 0.95
    ros_bypass_fraction: float = 0.01

    def as_array(self) -> np.ndarray:
        return np.array([
            self.h_per_nadh,
            self.h_per_fadh2,
            self.h_per_atp,
            self.atp_synthase_efficiency,
            self.membrane_leak_conductance,
            self.proton_slip_probability,
            self.roq_quinone_coupling,
            self.scaffold_channeling_efficiency,
            self.ros_bypass_fraction,
        ], dtype=np.float64)

    @classmethod
    def from_array(cls, arr: np.ndarray) -> DesignParameters:
        return cls(
            h_per_nadh=float(arr[0]),
            h_per_fadh2=float(arr[1]),
            h_per_atp=float(arr[2]),
            atp_synthase_efficiency=float(arr[3]),
            membrane_leak_conductance=float(arr[4]),
            proton_slip_probability=float(arr[5]),
            roq_quinone_coupling=float(arr[6]),
            scaffold_channeling_efficiency=float(arr[7]),
            ros_bypass_fraction=float(arr[8]),
        )

    def clamp(self) -> DesignParameters:
        """Clamp all parameters to their physical bounds."""
        bounds = get_param_bounds()
        arr = self.as_array()
        for i, (lo, hi) in enumerate(bounds):
            arr[i] = np.clip(arr[i], lo, hi)
        return DesignParameters.from_array(arr)

    def __repr__(self) -> str:
        return (
            f"DesignParameters(h_per_nadh={self.h_per_nadh:.2f}, "
            f"h_per_fadh2={self.h_per_fadh2:.2f}, "
            f"h_per_atp={self.h_per_atp:.2f}, "
            f"ATP_synth_eff={self.atp_synthase_efficiency:.3f}, "
            f"leak={self.membrane_leak_conductance:.4f}, "
            f"slip={self.proton_slip_probability:.4f}, "
            f"quinone={self.roq_quinone_coupling:.4f}, "
            f"channel={self.scaffold_channeling_efficiency:.4f}, "
            f"ROS={self.ros_bypass_fraction:.4f})"
        )


def get_param_bounds() -> List[Tuple[float, float]]:
    """Return [(lower, upper)] for each of the nine DesignParameters fields."""
    return [
        (10.0, 40.0),   # h_per_nadh
        (6.0, 26.0),    # h_per_fadh2
        (2.0, 5.0),     # h_per_atp
        (0.7, 1.0),     # atp_synthase_efficiency
        (0.0, 1.0),     # membrane_leak_conductance
        (0.0, 1.0),     # proton_slip_probability
        (0.0, 1.0),     # roq_quinone_coupling
        (0.0, 1.0),     # scaffold_channeling_efficiency
        (0.0, 1.0),     # ros_bypass_fraction
    ]


def random_design(rng: Optional[np.random.Generator] = None) -> DesignParameters:
    """Create a random DesignParameters within physical bounds."""
    if rng is None:
        rng = np.random.default_rng()
    bounds = get_param_bounds()
    arr = np.array([rng.uniform(lo, hi) for lo, hi in bounds], dtype=np.float64)
    return DesignParameters.from_array(arr)


# ─────────────────────────────────────────────────────────────────────────────
# Objective Functions  (all → float; minimisation convention)
# ─────────────────────────────────────────────────────────────────────────────


class ObjectiveFunction:
    """Static methods defining optimisation objectives for the PCM."""

    @staticmethod
    def compute_atp_yield(params: DesignParameters) -> float:
        """Compute the actual ATP yield per glucose given design parameters.

        The yield accounts for:
            - Protons pumped by NADH and FADH2 oxidation.
            - ATP synthase coupling efficiency.
            - Proton leakage, slippage, and ROS bypass losses.
            - Substrate channeling efficiency.
            - Substrate-level ATP (4 per glucose).
        """
        # Gross protons from redox pairs
        protons_nadh = params.h_per_nadh * NADH_COUNT
        protons_fadh2 = params.h_per_fadh2 * FADH2_COUNT
        total_protons = protons_nadh + protons_fadh2

        # Losses reduce effective proton motive force
        effective_protons = (
            total_protons
            * (1.0 - params.membrane_leak_conductance)
            * (1.0 - params.proton_slip_probability)
            * params.roq_quinone_coupling
        )

        # ATP from redox, limited by synthase and channeling
        redox_atp = (
            (effective_protons / params.h_per_atp)
            * params.atp_synthase_efficiency
            * params.scaffold_channeling_efficiency
        )

        # ROS bypass reduces available electrons
        redox_atp *= 1.0 - params.ros_bypass_fraction

        total_atp = SUBSTRATE_ATP + redox_atp
        # ASSUMPTION: ATP yield cannot exceed the thermodynamic maximum.
        return float(min(total_atp, THEORETICAL_MAX_ATP))

    @staticmethod
    def maximize_atp_yield(params: DesignParameters) -> float:
        """Return negative ATP yield for minimisation (maximise yield)."""
        return -ObjectiveFunction.compute_atp_yield(params)

    @staticmethod
    def maximize_thermodynamic_efficiency(params: DesignParameters) -> float:
        """Return negative thermodynamic efficiency (fraction of theoretical max).

        Thermodynamic efficiency = ATP_yield / THEORETICAL_MAX_ATP.
        """
        eff = ObjectiveFunction.compute_atp_yield(params) / THEORETICAL_MAX_ATP
        return -eff

    @staticmethod
    def maximize_carbon_efficiency(params: DesignParameters) -> float:
        """Return negative carbon efficiency.

        Carbon efficiency is approximated as ATP yield relative to the maximum
        possible yield, representing how much carbon is diverted to energy
        vs. lost.  High yield → high carbon efficiency.
        """
        yield_ratio = ObjectiveFunction.compute_atp_yield(params) / float(
            SUBSTRATE_ATP + NADH_COUNT * 8.0 + FADH2_COUNT * 7.5
        )
        return -min(yield_ratio, 1.0)

    @staticmethod
    def maximize_stability(params: DesignParameters) -> float:
        """Return negative stability (robustness) score.

        Stability is a heuristic measure based on how far parameters sit from
        extremes that cause fragility: low leakage, low slip, high coupling,
        moderate proton ratios.  Returns a score in [0, 1].
        """
        leak_penalty = params.membrane_leak_conductance
        slip_penalty = params.proton_slip_probability
        coupling_bonus = params.roq_quinone_coupling
        channel_bonus = params.scaffold_channeling_efficiency
        ros_penalty = params.ros_bypass_fraction

        # Ideal region: low leak/slip/ROS, high coupling/channeling
        stability = (
            (1.0 - leak_penalty)
            * (1.0 - slip_penalty)
            * coupling_bonus
            * channel_bonus
            * (1.0 - ros_penalty)
        )
        return -stability

    @staticmethod
    def minimize_entropy_production(params: DesignParameters) -> float:
        """Return entropy production metric (kJ/mol dissipated).

        Entropy production is the energy not captured as ATP:
            ΔS = ΔG_glucose - ATP_yield × ΔG_ATP
        """
        atp = ObjectiveFunction.compute_atp_yield(params)
        energy_captured = atp * DELTA_G_ATP_STD
        dissipated = max(DELTA_G_GLUCOSE - energy_captured, 0.0)
        return dissipated


# Alias mapping names to objective functions for convenience
OBJECTIVE_REGISTRY: Dict[str, Callable[[DesignParameters], float]] = {
    "maximize_atp_yield": ObjectiveFunction.maximize_atp_yield,
    "maximize_thermodynamic_efficiency": ObjectiveFunction.maximize_thermodynamic_efficiency,
    "maximize_carbon_efficiency": ObjectiveFunction.maximize_carbon_efficiency,
    "maximize_stability": ObjectiveFunction.maximize_stability,
    "minimize_entropy_production": ObjectiveFunction.minimize_entropy_production,
}


# ─────────────────────────────────────────────────────────────────────────────
# NSGA-II Helper Structures
# ─────────────────────────────────────────────────────────────────────────────


class Individual:
    """A single candidate solution in the genetic algorithm population."""

    def __init__(self, design: DesignParameters) -> None:
        self.design: DesignParameters = design
        self.objectives: List[float] = []
        self.rank: int = -1
        self.crowding_distance: float = 0.0

    def evaluate(self, objectives: Sequence[Callable[[DesignParameters], float]]) -> None:
        """Evaluate all objective functions and store results."""
        self.objectives = [obj(self.design) for obj in objectives]


class ParetoFront:
    """Result container for a multi-objective optimisation run."""

    def __init__(
        self,
        individuals: List[Individual],
        objective_names: List[str],
        n_generations: int,
        convergence_history: Optional[List[float]] = None,
    ) -> None:
        self.individuals = individuals
        self.objective_names = objective_names
        self.n_generations = n_generations
        self.convergence_history = convergence_history or []

    @property
    def designs(self) -> List[DesignParameters]:
        return [ind.design for ind in self.individuals]

    @property
    def objective_matrix(self) -> np.ndarray:
        return np.array([ind.objectives for ind in self.individuals])

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_generations": self.n_generations,
            "n_solutions": len(self.individuals),
            "objective_names": self.objective_names,
            "solutions": [
                {
                    "design": {
                        "h_per_nadh": ind.design.h_per_nadh,
                        "h_per_fadh2": ind.design.h_per_fadh2,
                        "h_per_atp": ind.design.h_per_atp,
                        "atp_synthase_efficiency": ind.design.atp_synthase_efficiency,
                        "membrane_leak_conductance": ind.design.membrane_leak_conductance,
                        "proton_slip_probability": ind.design.proton_slip_probability,
                        "roq_quinone_coupling": ind.design.roq_quinone_coupling,
                        "scaffold_channeling_efficiency": ind.design.scaffold_channeling_efficiency,
                        "ros_bypass_fraction": ind.design.ros_bypass_fraction,
                    },
                    "objectives": ind.objectives,
                }
                for ind in self.individuals
            ],
            "convergence_history": self.convergence_history,
        }

    def to_json(self, path: str) -> None:
        """Export results to a JSON file."""
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    def __repr__(self) -> str:
        return (
            f"ParetoFront(n_solutions={len(self.individuals)}, "
            f"n_objectives={len(self.objective_names)}, "
            f"n_generations={self.n_generations})"
        )


# ─────────────────────────────────────────────────────────────────────────────
# NSGA-II Core Operations
# ─────────────────────────────────────────────────────────────────────────────


def _non_dominated_sorting(population: List[Individual]) -> List[List[Individual]]:
    """Perform non-dominated sorting (NSGA-II fast sort).

    Returns a list of fronts, where front[0] is the Pareto-optimal set.
    """
    n = len(population)
    domination_count = [0] * n
    dominated_sets: List[List[int]] = [[] for _ in range(n)]
    fronts: List[List[int]] = [[]]

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if _dominates(population[i].objectives, population[j].objectives):
                dominated_sets[i].append(j)
            elif _dominates(population[j].objectives, population[i].objectives):
                domination_count[i] += 1
        if domination_count[i] == 0:
            fronts[0].append(i)
            population[i].rank = 0

    front_idx = 0
    while fronts[front_idx]:
        next_front: List[int] = []
        for i in fronts[front_idx]:
            for j in dominated_sets[i]:
                domination_count[j] -= 1
                if domination_count[j] == 0:
                    population[j].rank = front_idx + 1
                    next_front.append(j)
        front_idx += 1
        fronts.append(next_front)

    # Convert index lists to Individual lists
    result: List[List[Individual]] = []
    for f in fronts:
        if not f:
            break
        result.append([population[i] for i in f])
    return result


def _dominates(obj_a: Sequence[float], obj_b: Sequence[float]) -> bool:
    """Return True if a dominates b (all objectives <= and at least one <)."""
    at_least_one_better = False
    for a, b in zip(obj_a, obj_b):
        if a > b:
            return False
        if a < b:
            at_least_one_better = True
    return at_least_one_better


def _crowding_distance(front: List[Individual]) -> None:
    """Assign crowding distance to each individual in a front (in-place)."""
    if len(front) <= 2:
        for ind in front:
            ind.crowding_distance = float("inf")
        return

    n_obj = len(front[0].objectives)
    for ind in front:
        ind.crowding_distance = 0.0

    for m in range(n_obj):
        front.sort(key=lambda ind: ind.objectives[m])
        obj_min = front[0].objectives[m]
        obj_max = front[-1].objectives[m]
        norm = obj_max - obj_min
        if norm == 0:
            continue
        front[0].crowding_distance = float("inf")
        front[-1].crowding_distance = float("inf")
        for i in range(1, len(front) - 1):
            front[i].crowding_distance += (
                front[i + 1].objectives[m] - front[i - 1].objectives[m]
            ) / norm


def _tournament_selection(
    population: List[Individual], k: int, rng: np.random.Generator
) -> Individual:
    """Select one individual via binary tournament (rank then crowding)."""
    best = rng.integers(0, len(population))
    for _ in range(k - 1):
        candidate = rng.integers(0, len(population))
        a, b = population[best], population[candidate]
        if a.rank < b.rank or (
            a.rank == b.rank and a.crowding_distance > b.crowding_distance
        ):
            continue
        best = candidate
    return population[best]


def _sbx_crossover(
    parent1: DesignParameters,
    parent2: DesignParameters,
    eta_c: float = 20.0,
    prob_crossover: float = 0.9,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[DesignParameters, DesignParameters]:
    """Simulated Binary Crossover (SBX) for real-coded GA."""
    if rng is None:
        rng = np.random.default_rng()
    p1 = parent1.as_array()
    p2 = parent2.as_array()
    bounds = get_param_bounds()
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])
    c1, c2 = p1.copy(), p2.copy()

    for i in range(len(p1)):
        if rng.random() > prob_crossover:
            continue
        if abs(p1[i] - p2[i]) < 1e-14:
            continue
        y1, y2 = min(p1[i], p2[i]), max(p1[i], p2[i])
        beta = 1.0 + 2.0 * (y1 - lo[i]) / (y2 - y1 + 1e-14)
        alpha = 2.0 - beta ** (-(eta_c + 1.0))
        u = rng.random()
        if u <= 1.0 / alpha:
            beta_q = (u * alpha) ** (1.0 / (eta_c + 1.0))
        else:
            beta_q = (1.0 / (2.0 - u * alpha)) ** (1.0 / (eta_c + 1.0))
        c1_i = 0.5 * ((y1 + y2) - beta_q * (y2 - y1))
        c2_i = 0.5 * ((y1 + y2) + beta_q * (y2 - y1))
        c1_i = np.clip(c1_i, lo[i], hi[i])
        c2_i = np.clip(c2_i, lo[i], hi[i])
        c1[i] = c1_i
        c2[i] = c2_i

    return DesignParameters.from_array(c1), DesignParameters.from_array(c2)


def _polynomial_mutation(
    design: DesignParameters,
    eta_m: float = 20.0,
    prob_mutation: float = 0.1,
    rng: Optional[np.random.Generator] = None,
) -> DesignParameters:
    """Polynomial mutation for real-coded GA."""
    if rng is None:
        rng = np.random.default_rng()
    arr = design.as_array()
    bounds = get_param_bounds()
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])

    for i in range(len(arr)):
        if rng.random() > prob_mutation:
            continue
        delta = min(arr[i] - lo[i], hi[i] - arr[i]) / (hi[i] - lo[i] + 1e-14)
        u = rng.random()
        if u <= 0.5:
            delta_q = (2.0 * u + (1.0 - 2.0 * u) * (1.0 - delta) ** (eta_m + 1.0)) ** (
                1.0 / (eta_m + 1.0)
            ) - 1.0
        else:
            delta_q = 1.0 - (
                2.0 * (1.0 - u)
                + 2.0 * (u - 0.5) * (1.0 - delta) ** (eta_m + 1.0)
            ) ** (1.0 / (eta_m + 1.0))
        arr[i] += delta_q * (hi[i] - lo[i])
        arr[i] = np.clip(arr[i], lo[i], hi[i])

    return DesignParameters.from_array(arr)


# ─────────────────────────────────────────────────────────────────────────────
# ParetoOptimizer (NSGA-II)
# ─────────────────────────────────────────────────────────────────────────────


class ParetoOptimizer:
    """NSGA-II-inspired multi-objective genetic algorithm optimiser.

    Args:
        objectives: List of objective functions (minimisation).
        objective_names: Human-readable names for each objective.
        param_bounds: Per-parameter (lower, upper) bounds.
        pop_size: Population size per generation.
        n_generations: Number of generations to evolve.
        crossover_prob: SBX crossover probability.
        mutation_prob: Polynomial mutation probability (per variable).
        eta_c: SBX crossover distribution index.
        eta_m: Polynomial mutation distribution index.
        seed: Random seed for reproducibility.
    """

    def __init__(
        self,
        objectives: Sequence[Callable[[DesignParameters], float]],
        objective_names: Sequence[str],
        param_bounds: Sequence[Tuple[float, float]],
        pop_size: int = 100,
        n_generations: int = 200,
        crossover_prob: float = 0.9,
        mutation_prob: float = 0.15,
        eta_c: float = 20.0,
        eta_m: float = 20.0,
        seed: Optional[int] = None,
    ) -> None:
        self.objectives = list(objectives)
        self.objective_names = list(objective_names)
        self.param_bounds = list(param_bounds)
        self.pop_size = pop_size
        self.n_generations = n_generations
        self.crossover_prob = crossover_prob
        self.mutation_prob = mutation_prob
        self.eta_c = eta_c
        self.eta_m = eta_m
        self.rng = np.random.default_rng(seed)

        self.population: List[Individual] = []
        self.convergence_history: List[float] = []

    def _initialize_population(self) -> None:
        """Create initial random population."""
        self.population = []
        for _ in range(self.pop_size):
            design = random_design(self.rng)
            ind = Individual(design)
            ind.evaluate(self.objectives)
            self.population.append(ind)

    def _compute_convergence_metric(self) -> float:
        """Compute average objective vector norm as a simple convergence proxy."""
        obj_array = np.array([ind.objectives for ind in self.population])
        ideal = obj_array.min(axis=0)
        # ASSUMPTION: convergence = mean Euclidean distance to ideal point
        distances = np.sqrt(((obj_array - ideal[np.newaxis, :]) ** 2).sum(axis=1))
        return float(distances.mean())

    def _select_parents(self) -> List[Individual]:
        """Select parents for the next generation using tournament selection."""
        parents = []
        for _ in range(self.pop_size):
            parent = _tournament_selection(self.population, 2, self.rng)
            parents.append(parent)
        return parents

    def _create_offspring(self, parents: List[Individual]) -> List[Individual]:
        """Create offspring population via crossover and mutation."""
        offspring: List[Individual] = []
        for i in range(0, len(parents), 2):
            p1 = parents[i].design
            p2 = parents[(i + 1) % len(parents)].design
            c1, c2 = _sbx_crossover(p1, p2, self.eta_c, self.crossover_prob, self.rng)
            c1 = _polynomial_mutation(c1, self.eta_m, self.mutation_prob, self.rng)
            c2 = _polynomial_mutation(c2, self.eta_m, self.mutation_prob, self.rng)
            o1 = Individual(c1)
            o1.evaluate(self.objectives)
            offspring.append(o1)
            if len(offspring) < self.pop_size:
                o2 = Individual(c2)
                o2.evaluate(self.objectives)
                offspring.append(o2)
        return offspring[: self.pop_size]

    def _environmental_selection(
        self, combined: List[Individual]
    ) -> List[Individual]:
        """Select pop_size individuals from the combined pool via NSGA-II.

        Non-dominated sorting → crowding distance → truncate.
        """
        fronts = _non_dominated_sorting(combined)
        next_pop: List[Individual] = []
        for front in fronts:
            if len(next_pop) + len(front) <= self.pop_size:
                _crowding_distance(front)
                next_pop.extend(front)
            else:
                _crowding_distance(front)
                front.sort(key=lambda ind: ind.crowding_distance, reverse=True)
                remaining = self.pop_size - len(next_pop)
                next_pop.extend(front[:remaining])
                break
        return next_pop

    def optimize(self, verbose: bool = True) -> ParetoFront:
        """Run the NSGA-II optimisation loop.

        Args:
            verbose: If True, print progress every generation.

        Returns:
            A ParetoFront containing the final Pareto-optimal solutions.
        """
        t_start = time.perf_counter()

        self._initialize_population()
        self.convergence_history = [self._compute_convergence_metric()]

        for gen in range(1, self.n_generations + 1):
            parents = self._select_parents()
            offspring = self._create_offspring(parents)
            combined = self.population + offspring
            self.population = self._environmental_selection(combined)

            metric = self._compute_convergence_metric()
            self.convergence_history.append(metric)

            if verbose and (gen % 20 == 0 or gen == self.n_generations):
                fronts = _non_dominated_sorting(self.population)
                n_pareto = len(fronts[0])
                elapsed = time.perf_counter() - t_start
                print(
                    f"Gen {gen:4d}/{self.n_generations} | "
                    f"Pareto-front size: {n_pareto:3d} | "
                    f"Convergence: {metric:.4e} | "
                    f"Elapsed: {elapsed:.1f}s"
                )

        # Final Pareto front
        fronts = _non_dominated_sorting(self.population)
        pareto_individuals = fronts[0]

        if verbose:
            elapsed = time.perf_counter() - t_start
            print(f"\nOptimisation complete in {elapsed:.1f}s")
            print(f"Pareto-optimal solutions found: {len(pareto_individuals)}")

        return ParetoFront(
            individuals=pareto_individuals,
            objective_names=list(self.objective_names),
            n_generations=self.n_generations,
            convergence_history=self.convergence_history,
        )


# ─────────────────────────────────────────────────────────────────────────────
# Reaction Network Evolution
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class ReactionRule:
    """A single reaction rule in a coupling architecture.

    Attributes:
        reactants: Dict mapping reactant name → stoichiometric coefficient.
        products: Dict mapping product name → stoichiometric coefficient.
        h_pumped: Number of protons pumped per reaction.
        atp_generated: ATP directly generated (substrate-level).
        rate_constant: Relative rate constant (affects flux).
        label: Human-readable identifier.
    """

    reactants: Dict[str, float] = field(default_factory=dict)
    products: Dict[str, float] = field(default_factory=dict)
    h_pumped: float = 0.0
    atp_generated: float = 0.0
    rate_constant: float = 1.0
    label: str = ""

    def __repr__(self) -> str:
        lhs = " + ".join(f"{v}{k}" for k, v in self.reactants.items())
        rhs = " + ".join(f"{v}{k}" for k, v in self.products.items())
        return f"{self.label}: {lhs} -> {rhs} (H⁺={self.h_pumped}, ATP={self.atp_generated})"


class ReactionNetwork:
    """A network of reaction rules representing a coupling architecture."""

    def __init__(
        self, rules: Optional[List[ReactionRule]] = None
    ) -> None:
        self.rules: List[ReactionRule] = rules or []

    def add_rule(self, rule: ReactionRule) -> None:
        self.rules.append(rule)

    def remove_rule(self, index: int) -> None:
        if 0 <= index < len(self.rules):
            self.rules.pop(index)

    def evaluate_fitness(
        self,
        nadh_available: float = 10.0,
        fadh2_available: float = 2.0,
        o2_available: float = 6.0,
    ) -> Dict[str, float]:
        """Evaluate how well this reaction network performs given redox inputs.

        Returns a dict with 'atp_yield', 'efficiency', 'proton_flux',
        'carbon_efficiency', and 'stability'.
        """
        total_atp = float(SUBSTRATE_ATP)
        total_h_pumped = 0.0
        flux_balance = True

        for rule in self.rules:
            total_atp += rule.atp_generated * rule.rate_constant
            total_h_pumped += rule.h_pumped * rule.rate_constant
            # ASSUMPTION: check rough stoichiometric conservation
            if abs(sum(rule.reactants.values()) - sum(rule.products.values())) > 1e-6:
                flux_balance = False

        # ASSUMPTION: proton-driven ATP estimation
        atp_from_protons = (
            total_h_pumped / 3.0
        )  # nominal H⁺/ATP = 3
        total_atp += atp_from_protons * 0.5  # discount for coupling losses

        theoretical = THEORETICAL_MAX_ATP
        efficiency = min(total_atp / theoretical, 1.0)
        stability = 0.5 if flux_balance else 0.1

        return {
            "atp_yield": total_atp,
            "efficiency": efficiency,
            "proton_flux": total_h_pumped,
            "carbon_efficiency": efficiency,
            "stability": stability,
        }

    def copy(self) -> ReactionNetwork:
        return ReactionNetwork(rules=copy.deepcopy(self.rules))


class ReactionNetworkEvolution:
    """Evolutionary search over reaction network topologies.

    Evolves reaction rule sets (coupling architectures) toward higher ATP
    yield and thermodynamic efficiency.  Crossover mixes rule sets between
    parents; mutation adds, removes, or modifies individual reaction rules.
    """

    def __init__(
        self,
        pop_size: int = 50,
        n_generations: int = 100,
        mutation_rate: float = 0.3,
        crossover_rate: float = 0.7,
        seed: Optional[int] = None,
    ) -> None:
        self.pop_size = pop_size
        self.n_generations = n_generations
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.rng = np.random.default_rng(seed)

        self.population: List[ReactionNetwork] = []
        self.fitness_history: List[float] = []

    def _random_network(self) -> ReactionNetwork:
        """Generate a random reaction network with plausible PCM rules."""
        network = ReactionNetwork()

        # Base PCM reactions
        n_rules = self.rng.integers(2, 6)
        for i in range(n_rules):
            rule_type = self.rng.choice(["nadh_oxidase", "fadh2_oxidase", "atp_synthase", "quinone_cycle"])
            if rule_type == "nadh_oxidase":
                h = float(self.rng.uniform(10, 40))
                rule = ReactionRule(
                    reactants={"NADH": 1.0, "Q": 1.0},
                    products={"NAD+": 1.0, "QH2": 1.0},
                    h_pumped=h,
                    atp_generated=0.0,
                    rate_constant=self.rng.uniform(0.5, 1.0),
                    label=f"NADH_ox_{i}",
                )
            elif rule_type == "fadh2_oxidase":
                h = float(self.rng.uniform(6, 26))
                rule = ReactionRule(
                    reactants={"FADH2": 1.0, "Q": 1.0},
                    products={"FAD": 1.0, "QH2": 1.0},
                    h_pumped=h,
                    atp_generated=0.0,
                    rate_constant=self.rng.uniform(0.5, 1.0),
                    label=f"FADH2_ox_{i}",
                )
            elif rule_type == "atp_synthase":
                h_consumed = float(self.rng.uniform(2, 5))
                rule = ReactionRule(
                    reactants={"ADP": 1.0, "Pi": 1.0, f"{int(h_consumed)}H+_in": 1.0},
                    products={"ATP": 1.0, "H2O": 1.0, f"{int(h_consumed)}H+_out": 1.0},
                    h_pumped=-h_consumed,
                    atp_generated=1.0,
                    rate_constant=self.rng.uniform(0.5, 1.0),
                    label=f"ATP_synth_{i}",
                )
            else:  # quinone_cycle
                rule = ReactionRule(
                    reactants={"QH2": 1.0, "O2": 0.5},
                    products={"Q": 1.0, "H2O": 1.0},
                    h_pumped=float(self.rng.uniform(4, 12)),
                    atp_generated=0.0,
                    rate_constant=self.rng.uniform(0.5, 1.0),
                    label=f"Q_cycle_{i}",
                )
            network.add_rule(rule)

        return network

    def _initialize(self) -> None:
        """Seed population with random reaction networks."""
        self.population = [self._random_network() for _ in range(self.pop_size)]

    def _fitness(self, network: ReactionNetwork) -> float:
        """Combined scalar fitness (higher is better)."""
        metrics = network.evaluate_fitness()
        return metrics["atp_yield"] * metrics["efficiency"] * (1.0 + metrics["stability"])

    def _mutate(self, network: ReactionNetwork) -> ReactionNetwork:
        """Mutate a reaction network: add, remove, or modify a rule."""
        mutant = network.copy()
        roll = self.rng.random()

        if roll < 0.33 and len(mutant.rules) > 0:
            # Remove a rule
            idx = self.rng.integers(0, len(mutant.rules))
            mutant.remove_rule(idx)
        elif roll < 0.66 and len(mutant.rules) < 12:
            # Add a new random rule
            new_rule = self._random_network().rules[0]
            new_rule.label = f"mut_{len(mutant.rules)}"
            mutant.add_rule(new_rule)
        else:
            # Modify an existing rule
            if len(mutant.rules) == 0:
                return mutant
            idx = self.rng.integers(0, len(mutant.rules))
            rule = mutant.rules[idx]
            param = self.rng.choice(["h_pumped", "atp_generated", "rate_constant"])
            if param == "h_pumped":
                rule.h_pumped *= self.rng.uniform(0.8, 1.2)
                rule.h_pumped = max(0.0, rule.h_pumped)
            elif param == "atp_generated":
                rule.atp_generated *= self.rng.uniform(0.8, 1.2)
                rule.atp_generated = max(0.0, rule.atp_generated)
            else:
                rule.rate_constant = np.clip(
                    rule.rate_constant * self.rng.uniform(0.5, 2.0), 0.1, 3.0
                )

        return mutant

    def _crossover(
        self, parent1: ReactionNetwork, parent2: ReactionNetwork
    ) -> Tuple[ReactionNetwork, ReactionNetwork]:
        """One-point crossover: split rule lists at a random point and swap."""
        split1 = self.rng.integers(0, len(parent1.rules) + 1)
        split2 = self.rng.integers(0, len(parent2.rules) + 1)
        child1 = ReactionNetwork(
            rules=parent1.rules[:split1] + parent2.rules[split2:]
        )
        child2 = ReactionNetwork(
            rules=parent2.rules[:split2] + parent1.rules[split1:]
        )
        return child1, child2

    def _tournament_select(self) -> ReactionNetwork:
        """Tournament selection (k=3) based on fitness."""
        best = self.rng.integers(0, len(self.population))
        best_f = self._fitness(self.population[best])
        for _ in range(2):
            idx = self.rng.integers(0, len(self.population))
            f = self._fitness(self.population[idx])
            if f > best_f:
                best = idx
                best_f = f
        return self.population[best]

    def evolve(self, verbose: bool = True) -> List[ReactionNetwork]:
        """Run the evolutionary search.

        Args:
            verbose: Print progress every generation.

        Returns:
            The final population of reaction networks (sorted by fitness).
        """
        self._initialize()

        for gen in range(self.n_generations):
            new_pop: List[ReactionNetwork] = []

            # Elitism: keep top 2
            ranked = sorted(self.population, key=lambda n: self._fitness(n), reverse=True)
            new_pop.extend(ranked[:2])

            while len(new_pop) < self.pop_size:
                p1 = self._tournament_select()
                p2 = self._tournament_select()

                if self.rng.random() < self.crossover_rate:
                    c1, c2 = self._crossover(p1, p2)
                else:
                    c1, c2 = p1.copy(), p2.copy()

                if self.rng.random() < self.mutation_rate:
                    c1 = self._mutate(c1)
                if self.rng.random() < self.mutation_rate:
                    c2 = self._mutate(c2)

                new_pop.append(c1)
                if len(new_pop) < self.pop_size:
                    new_pop.append(c2)

            self.population = new_pop[: self.pop_size]

            best_f = self._fitness(ranked[0])
            self.fitness_history.append(best_f)
            if verbose and (gen % 20 == 0 or gen == self.n_generations - 1):
                print(f"Network evolution gen {gen:4d} | Best fitness: {best_f:.4f}")

        ranked = sorted(self.population, key=lambda n: self._fitness(n), reverse=True)
        return ranked


# ─────────────────────────────────────────────────────────────────────────────
# Convenience Runner
# ─────────────────────────────────────────────────────────────────────────────


def run_optimization(
    objectives: Optional[List[str]] = None,
    param_ranges: Optional[Sequence[Tuple[float, float]]] = None,
    n_generations: int = 200,
    pop_size: int = 100,
    seed: Optional[int] = None,
    verbose: bool = True,
) -> ParetoFront:
    """Run multi-objective optimisation of PCM design parameters.

    Args:
        objectives: List of objective names from OBJECTIVE_REGISTRY.  Defaults
            to all five objectives.
        param_ranges: Per-parameter (lower, upper) bounds.  Defaults to the
            standard DesignParameters bounds.
        n_generations: Number of generations.
        pop_size: Population size.
        seed: Random seed.
        verbose: Print progress.

    Returns:
        ParetoFront with the optimal trade-off solutions.
    """
    if objectives is None:
        objectives = list(OBJECTIVE_REGISTRY.keys())
    objective_fns = [OBJECTIVE_REGISTRY[name] for name in objectives]

    if param_ranges is None:
        param_ranges = get_param_bounds()

    optimizer = ParetoOptimizer(
        objectives=objective_fns,
        objective_names=objectives,
        param_bounds=param_ranges,
        pop_size=pop_size,
        n_generations=n_generations,
        seed=seed,
    )
    return optimizer.optimize(verbose=verbose)


# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────


def plot_pareto_front(
    pareto_front: ParetoFront,
    objectives: Optional[List[str]] = None,
    save_path: Optional[str] = None,
) -> None:
    """Plot the Pareto front as a pairplot of all objective combinations.

    Args:
        pareto_front: Results from run_optimization.
        objectives: Subset of objective names to plot.  Defaults to all.
        save_path: If provided, save figure to this path.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plot.")
        return

    if objectives is None:
        objectives = pareto_front.objective_names
    n = len(objectives)
    obj_matrix = pareto_front.objective_matrix
    name_indices = {
        name: i for i, name in enumerate(pareto_front.objective_names)
    }
    selected = [name_indices[name] for name in objectives]
    data = obj_matrix[:, selected]

    fig, axes = plt.subplots(n, n, figsize=(3 * n, 3 * n))
    if n == 1:
        axes = np.array([[axes]])

    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            if i == j:
                ax.hist(data[:, i], bins=20, color="steelblue", edgecolor="white")
                ax.set_xlabel(objectives[i])
                ax.set_ylabel("Count")
            else:
                ax.scatter(data[:, j], data[:, i], s=15, alpha=0.7, c="crimson")
                ax.set_xlabel(objectives[j])
                ax.set_ylabel(objectives[i])

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Pareto front plot saved to {save_path}")
    else:
        plt.show()
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Usage
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Quick demonstration
    print("=" * 60)
    print("PTR-94 Pareto Optimizer — Quick Demo")
    print("=" * 60)

    result = run_optimization(
        objectives=["maximize_atp_yield", "maximize_thermodynamic_efficiency", "minimize_entropy_production"],
        n_generations=50,
        pop_size=40,
        seed=42,
    )
    print(result)
    print("\nPareto-optimal designs:")
    for ind in result.individuals[:5]:
        atp = ObjectiveFunction.compute_atp_yield(ind.design)
        print(f"  {ind.design} → ATP yield = {atp:.2f}")
