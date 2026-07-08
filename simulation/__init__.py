#!/usr/bin/env python3
"""
PTR-94 Simulation Package.

Computational modelling, multi-objective optimisation, sensitivity analysis,
and uncertainty quantification for the Perfect Coupling Module (PCM) of the
PTR-94 pathway.
"""

from .thermodynamics import (
    ThermodynamicState,
    compute_delta_g,
    compute_equilibrium_constant,
    standard_conditions_correction,
)
from .energy_balance import (
    EnergyBalance,
    compute_atp_yield,
    compute_proton_motive_force,
    compute_redox_potential,
)
from .reaction_network import (
    ReactionNetwork,
    ReactionRule,
    ReactionNetworkEvolution,
)
from .kinetics import (
    KineticModel,
    michaelis_menten_rate,
    hill_rate,
    compute_flux,
)
from .pareto_optimizer import (
    DesignParameters,
    ObjectiveFunction,
    ParetoOptimizer,
    ParetoFront,
    Individual,
    run_optimization,
    plot_pareto_front,
    get_param_bounds,
    random_design,
    OBJECTIVE_REGISTRY,
)
from .sensitivity_analysis import (
    LocalSensitivity,
    GlobalSensitivity,
    SensitivityReport,
    sensitivity_of_atp_yield,
    sensitivity_of_efficiency,
    sensitivity_of_entropy,
    parameter_ranking,
    plot_sensitivity,
    build_sensitivity_report,
)
from .monte_carlo import (
    ParameterDistribution,
    MonteCarloEngine,
    MonteCarloResult,
    UncertaintyQuantification,
    UncertaintyResult,
    monte_carlo_atp_yield,
    monte_carlo_efficiency,
    monte_carlo_proton_flux,
    confidence_intervals,
    plot_uncertainty,
)

__all__ = [
    "DesignParameters",
    "ObjectiveFunction",
    "ParetoOptimizer",
    "ParetoFront",
    "Individual",
    "ReactionRule",
    "ReactionNetwork",
    "ReactionNetworkEvolution",
    "run_optimization",
    "plot_pareto_front",
    "get_param_bounds",
    "random_design",
    "OBJECTIVE_REGISTRY",
    "LocalSensitivity",
    "GlobalSensitivity",
    "SensitivityReport",
    "sensitivity_of_atp_yield",
    "sensitivity_of_efficiency",
    "sensitivity_of_entropy",
    "parameter_ranking",
    "plot_sensitivity",
    "build_sensitivity_report",
    "ParameterDistribution",
    "MonteCarloEngine",
    "MonteCarloResult",
    "UncertaintyQuantification",
    "UncertaintyResult",
    "monte_carlo_atp_yield",
    "monte_carlo_efficiency",
    "monte_carlo_proton_flux",
    "confidence_intervals",
    "plot_uncertainty",
]
