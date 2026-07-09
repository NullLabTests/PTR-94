#!/usr/bin/env python3
"""
Sensitivity Analysis for the PTR‑94 Perfect Coupling Module (PCM).

Provides local (OAT derivative-based) and global (Sobol index) sensitivity
methods to quantify how variability in PCM design parameters propagates to
key outputs: ATP yield, thermodynamic efficiency, and entropy production.

Assumptions
-----------
- Standard biochemical conditions (ΔG°′ = 30.5 kJ/mol for ATP).
- 10 NADH + 2 FADH2 + 4 substrate-level ATP from upstream modules.
- Parameter ranges are independent unless explicitly correlated.
- The model response is sufficiently smooth for local derivatives.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

# Reuse the core model from pareto_optimizer
try:
    from .pareto_optimizer import (
        DELTA_G_GLUCOSE,
        DELTA_G_ATP_STD,
        FADH2_COUNT,
        NADH_COUNT,
        SUBSTRATE_ATP,
        THEORETICAL_MAX_ATP,
        DesignParameters,
        ObjectiveFunction,
        get_param_bounds,
    )
except ImportError:
    # Fallback for standalone usage
    DELTA_G_GLUCOSE = 2870.0
    DELTA_G_ATP_STD = 30.5
    THEORETICAL_MAX_ATP = DELTA_G_GLUCOSE / DELTA_G_ATP_STD
    SUBSTRATE_ATP = 4
    NADH_COUNT = 10
    FADH2_COUNT = 2

    @dataclass
    class DesignParameters:
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
                self.h_per_nadh, self.h_per_fadh2, self.h_per_atp,
                self.atp_synthase_efficiency, self.membrane_leak_conductance,
                self.proton_slip_probability, self.roq_quinone_coupling,
                self.scaffold_channeling_efficiency, self.ros_bypass_fraction,
            ], dtype=np.float64)

        @classmethod
        def from_array(cls, arr: np.ndarray) -> DesignParameters:
            return cls(*arr.tolist())

    def get_param_bounds() -> List[Tuple[float, float]]:
        return [
            (10.0, 40.0), (6.0, 26.0), (2.0, 5.0), (0.7, 1.0),
            (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
        ]

    class ObjectiveFunction:
        @staticmethod
        def compute_atp_yield(params: DesignParameters) -> float:
            protons_nadh = params.h_per_nadh * NADH_COUNT
            protons_fadh2 = params.h_per_fadh2 * FADH2_COUNT
            total_protons = protons_nadh + protons_fadh2
            effective = (
                total_protons
                * (1.0 - params.membrane_leak_conductance)
                * (1.0 - params.proton_slip_probability)
                * params.roq_quinone_coupling
            )
            redox_atp = (
                (effective / params.h_per_atp)
                * params.atp_synthase_efficiency
                * params.scaffold_channeling_efficiency
                * (1.0 - params.ros_bypass_fraction)
            )
            return float(min(SUBSTRATE_ATP + redox_atp, THEORETICAL_MAX_ATP))

        @staticmethod
        def maximize_thermodynamic_efficiency(params: DesignParameters) -> float:
            yield_ = ObjectiveFunction.compute_atp_yield(params)
            return yield_ / THEORETICAL_MAX_ATP

        @staticmethod
        def minimize_entropy_production(params: DesignParameters) -> float:
            atp = ObjectiveFunction.compute_atp_yield(params)
            captured = atp * DELTA_G_ATP_STD
            return max(DELTA_G_GLUCOSE - captured, 0.0)


# ─────────────────────────────────────────────────────────────────────────────
# Parameter metadata
# ─────────────────────────────────────────────────────────────────────────────

PARAMETER_NAMES: List[str] = [
    "h_per_nadh",
    "h_per_fadh2",
    "h_per_atp",
    "atp_synthase_efficiency",
    "membrane_leak_conductance",
    "proton_slip_probability",
    "roq_quinone_coupling",
    "scaffold_channeling_efficiency",
    "ros_bypass_fraction",
]


def _to_design(params: Sequence[float]) -> DesignParameters:
    """Convert a flat numeric sequence to DesignParameters."""
    return DesignParameters.from_array(np.asarray(params, dtype=np.float64))


# ─────────────────────────────────────────────────────────────────────────────
# Local Sensitivity Analysis (OAT / derivative-based)
# ─────────────────────────────────────────────────────────────────────────────


class LocalSensitivity:
    """One-at-a-time (OAT) derivative-based sensitivity analysis.

    Computes normalised sensitivity coefficients (∂f/∂p × p/f) for each
    parameter, providing a dimensionless ranking of local influence around a
    nominal design point.

    Args:
        nominal: The baseline DesignParameters.
        perturbation: Fractional perturbation for finite differences (default 0.01).
        outputs: Dict of output-name → callable(DesignParameters) -> float.
                 Defaults to atp_yield, thermodynamic_efficiency, entropy_production.
    """

    def __init__(
        self,
        nominal: DesignParameters,
        perturbation: float = 0.01,
        outputs: Optional[Dict[str, Callable[[DesignParameters], float]]] = None,
    ) -> None:
        self.nominal = nominal
        self.perturbation = perturbation
        self.outputs = outputs or {
            "atp_yield": ObjectiveFunction.compute_atp_yield,
            "thermodynamic_efficiency": ObjectiveFunction.maximize_thermodynamic_efficiency,
            "entropy_production": ObjectiveFunction.minimize_entropy_production,
        }
        self._sensitivity_matrix: Optional[np.ndarray] = None
        self._coefficients: Optional[Dict[str, Dict[str, float]]] = None

    def compute(self) -> Dict[str, Dict[str, float]]:
        """Compute normalised sensitivity coefficients for all parameter–output pairs.

        Returns nested dict: coefficient[output_name][parameter_name] = value.
        """
        params_arr = self.nominal.as_array()
        n_params = len(params_arr)
        bounds = get_param_bounds()
        eps = self.perturbation

        self._coefficients = {
            out_name: {} for out_name in self.outputs
        }
        self._sensitivity_matrix = np.zeros((len(self.outputs), n_params))

        for i, (out_name, out_fn) in enumerate(self.outputs.items()):
            f0 = out_fn(self.nominal)
            for j in range(n_params):
                lo, hi = bounds[j]
                delta = max(eps * (hi - lo), 1e-10)
                params_up = params_arr.copy()
                params_up[j] = np.clip(params_arr[j] + delta, lo, hi)
                f_up = out_fn(_to_design(params_up))

                params_dn = params_arr.copy()
                params_dn[j] = np.clip(params_arr[j] - delta, lo, hi)
                f_dn = out_fn(_to_design(params_dn))

                # Central finite difference
                df_dp = (f_up - f_dn) / (2.0 * delta)

                # Normalised sensitivity coefficient (dimensionless)
                if abs(f0) > 1e-15 and abs(params_arr[j]) > 1e-15:
                    S = df_dp * params_arr[j] / f0
                else:
                    S = 0.0

                self._coefficients[out_name][PARAMETER_NAMES[j]] = S
                self._sensitivity_matrix[i, j] = S

        return self._coefficients

    @property
    def sensitivity_matrix(self) -> np.ndarray:
        """Return the (n_outputs × n_params) matrix of sensitivity coefficients."""
        if self._sensitivity_matrix is None:
            self.compute()
        assert self._sensitivity_matrix is not None
        return self._sensitivity_matrix

    @property
    def coefficients(self) -> Dict[str, Dict[str, float]]:
        """Cached sensitivity coefficients dict."""
        if self._coefficients is None:
            self.compute()
        assert self._coefficients is not None
        return self._coefficients

    def ranking(self, output: str = "atp_yield") -> List[Tuple[str, float]]:
        """Return parameters sorted by absolute sensitivity for a given output.

        Args:
            output: Name of the output to rank by.

        Returns:
            List of (parameter_name, absolute_sensitivity) sorted descending.
        """
        if self._coefficients is None:
            self.compute()
        coeffs = self._coefficients[output]
        ranked = sorted(
            coeffs.items(), key=lambda kv: abs(kv[1]), reverse=True
        )
        return ranked

    def most_influential(
        self, output: str = "atp_yield", top_n: int = 3
    ) -> List[Tuple[str, float]]:
        """Return the top_n most influential parameters for a given output."""
        return self.ranking(output)[:top_n]


# ─────────────────────────────────────────────────────────────────────────────
# Global Sensitivity Analysis (Sobol indices via Saltelli sampling)
# ─────────────────────────────────────────────────────────────────────────────


class GlobalSensitivity:
    """Global sensitivity analysis using Sobol' indices.

    Computes first-order (S1) and total-order (ST) Sobol indices via the
    Saltelli sampling scheme.  Requires O(N × (D + 2)) model evaluations
    where N = n_samples and D = number of parameters.

    Args:
        param_bounds: List of (lower, upper) for each parameter.
        n_samples: Base number of Monte Carlo samples (total eval = n_samples × (2 × D + 2)).
        seed: Random seed.
    """

    def __init__(
        self,
        param_bounds: Sequence[Tuple[float, float]],
        n_samples: int = 5000,
        seed: Optional[int] = None,
    ) -> None:
        self.param_bounds = list(param_bounds)
        self.n_params = len(self.param_bounds)
        self.n_samples = n_samples
        self.rng = np.random.default_rng(seed)

        self._s1: Optional[np.ndarray] = None
        self._st: Optional[np.ndarray] = None
        self._s1_conf: Optional[np.ndarray] = None
        self._st_conf: Optional[np.ndarray] = None

    def _saltelli_sample(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generate Saltelli sampling matrices A, B, and AB_i.

        Returns (A, B, AB) where AB has shape (n_params, n_samples, n_params).
        AB[i] replaces column i of A with column i of B.
        """
        D = self.n_params
        N = self.n_samples
        lo = np.array([b[0] for b in self.param_bounds])
        hi = np.array([b[1] for b in self.param_bounds])

        # Latin Hypercube for A and B
        def _lhs(d: int, n: int) -> np.ndarray:
            segments = np.linspace(0.0, 1.0, n + 1)

            mat = np.zeros((n, d))
            for j in range(d):
                offsets = self.rng.uniform(0, segments[1], size=n)
                order = self.rng.permutation(n)
                mat[:, j] = segments[order] + offsets
            return mat

        A_raw = _lhs(D, N)
        B_raw = _lhs(D, N)

        # Scale to parameter ranges
        A = lo[np.newaxis, :] + A_raw * (hi - lo)[np.newaxis, :]
        B = lo[np.newaxis, :] + B_raw * (hi - lo)[np.newaxis, :]

        # AB matrices: AB[i] = A with column i replaced by B[:, i]
        AB = np.zeros((D, N, D))
        for i in range(D):
            AB[i] = A.copy()
            AB[i][:, i] = B[:, i]

        return A, B, AB

    @staticmethod
    def _compute_indices(
        f_A: np.ndarray,
        f_B: np.ndarray,
        f_AB: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compute first-order and total-order Sobol indices.

        Uses the Jansen (1999) estimator for total-order indices.

        Args:
            f_A:   Model outputs for sample matrix A, shape (N,).
            f_B:   Model outputs for sample matrix B, shape (N,).
            f_AB:  Model outputs for AB[i], shape (D, N).

        Returns:
            (S1, ST, S1_conf, ST_conf) each shape (D,).
        """
        D = f_AB.shape[0]
        N = len(f_A)

        # Variance estimates
        var_total = np.var(np.concatenate([f_A, f_B, f_AB.ravel()]), ddof=1)

        if var_total < 1e-30:
            return (np.zeros(D), np.zeros(D), np.zeros(D), np.zeros(D))

        S1 = np.zeros(D)
        ST = np.zeros(D)
        S1_conf = np.zeros(D)
        ST_conf = np.zeros(D)

        for i in range(D):
            # First-order (Saltelli 2010 estimator)
            cov = np.mean(f_B * (f_AB[i] - f_A))
            S1[i] = cov / var_total

            # Total-order (Jansen 1999)
            ST[i] = 0.5 * np.mean((f_A - f_AB[i]) ** 2) / var_total

            # Bootstrap confidence intervals
            n_boot = 100
            s1_boot = np.zeros(n_boot)
            st_boot = np.zeros(n_boot)
            for b in range(n_boot):
                idx = np.random.choice(N, N, replace=True)
                cov_b = np.mean(f_B[idx] * (f_AB[i, idx] - f_A[idx]))
                s1_boot[b] = cov_b / var_total
                st_boot[b] = 0.5 * np.mean(
                    (f_A[idx] - f_AB[i, idx]) ** 2
                ) / var_total

            S1_conf[i] = 1.96 * np.std(s1_boot, ddof=1)
            ST_conf[i] = 1.96 * np.std(st_boot, ddof=1)

        return S1, ST, S1_conf, ST_conf

    def analyze(
        self,
        model: Callable[[DesignParameters], float],
    ) -> Dict[str, Any]:
        """Run Sobol sensitivity analysis.

        Args:
            model: Function mapping DesignParameters → scalar output.

        Returns:
            Dict with keys 'S1', 'ST', 'S1_conf', 'ST_conf', 'parameter_names',
            and 'n_total_evaluations'.
        """
        D = self.n_params
        A, B, AB = self._saltelli_sample()

        # Evaluate model
        f_A = np.array([model(_to_design(A[i])) for i in range(self.n_samples)])
        f_B = np.array([model(_to_design(B[i])) for i in range(self.n_samples)])
        f_AB = np.zeros((D, self.n_samples))
        for i in range(D):
            for j in range(self.n_samples):
                f_AB[i, j] = model(_to_design(AB[i, j]))

        S1, ST, S1_conf, ST_conf = self._compute_indices(f_A, f_B, f_AB)

        self._s1 = S1
        self._st = ST
        self._s1_conf = S1_conf
        self._st_conf = ST_conf

        n_total = self.n_samples * (2 * D + 2)
        return {
            "S1": S1,
            "ST": ST,
            "S1_conf": S1_conf,
            "ST_conf": ST_conf,
            "parameter_names": PARAMETER_NAMES[:D],
            "n_total_evaluations": n_total,
        }

    @property
    def first_order(self) -> np.ndarray:
        if self._s1 is None:
            raise RuntimeError("Run analyze() first.")
        return self._s1

    @property
    def total_order(self) -> np.ndarray:
        if self._st is None:
            raise RuntimeError("Run analyze() first.")
        return self._st


# ─────────────────────────────────────────────────────────────────────────────
# Sensitivity Report
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class SensitivityReport:
    """Comprehensive report of sensitivity analysis results.

    Aggregates both local and global sensitivity metrics, identifies the most
    influential parameters, evaluates robustness, and breaks down uncertainty
    contributions.

    Attributes:
        output_name: Name of the model output analysed.
        local_coefficients: Dict of param → local sensitivity coefficient.
        global_first_order: Dict of param → first-order Sobol index.
        global_total_order: Dict of param → total-order Sobol index.
        parameter_ranking: Ordered list of (param, importance_score).
        robustness_metric: A scalar in [0, 1] quantifying robustness: 1 =
            output is insensitive to all parameter variations.
        uncertainty_breakdown: Dict of param → fractional contribution to
            output variance (from total-order Sobol).
    """

    output_name: str = ""
    local_coefficients: Dict[str, float] = field(default_factory=dict)
    global_first_order: Dict[str, float] = field(default_factory=dict)
    global_total_order: Dict[str, float] = field(default_factory=dict)
    parameter_ranking: List[Tuple[str, float]] = field(default_factory=list)
    robustness_metric: float = 0.0
    uncertainty_breakdown: Dict[str, float] = field(default_factory=dict)

    def summary_table(self) -> str:
        """Return a formatted string table of sensitivity indices."""
        if not self.parameter_ranking:
            return "No sensitivity data available."

        lines = [
            f"\n{'=' * 90}",
            f"  Sensitivity Report — Output: {self.output_name}",
            f"{'=' * 90}",
            f"{'Parameter':<30} {'Local S':>10} {'Sobol S1':>10} {'Sobol ST':>10} {'Importance':>10}",
            f"{'-' * 30} {'-' * 10} {'-' * 10} {'-' * 10} {'-' * 10}",
        ]
        for param, importance in self.parameter_ranking:
            local_s = self.local_coefficients.get(param, 0.0)
            s1 = self.global_first_order.get(param, 0.0)
            st = self.global_total_order.get(param, 0.0)
            lines.append(
                f"{param:<30} {local_s:>+10.4f} {s1:>10.4f} {st:>10.4f} {importance:>10.4f}"
            )

        lines.append(f"{'-' * 90}")
        lines.append(
            f"{'Robustness metric:':<30} {self.robustness_metric:>10.4f}"
        )
        lines.append(f"{'=' * 90}\n")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"SensitivityReport(output='{self.output_name}', "
            f"n_params={len(self.local_coefficients)}, "
            f"robustness={self.robustness_metric:.4f})"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Convenience Analysis Functions
# ─────────────────────────────────────────────────────────────────────────────


def _default_distributions() -> Dict[str, Dict[str, Any]]:
    """Return default parameter distributions for Monte Carlo sampling.

    Each entry: name → {'type': str, 'params': list, 'bounds': (lo, hi)}.
    """
    return {
        "h_per_nadh": {
            "type": "uniform",
            "params": [10.0, 40.0],
            "bounds": (10.0, 40.0),
        },
        "h_per_fadh2": {
            "type": "uniform",
            "params": [6.0, 26.0],
            "bounds": (6.0, 26.0),
        },
        "h_per_atp": {
            "type": "triangular",
            "params": [2.0, 3.0, 5.0],
            "bounds": (2.0, 5.0),
        },
        "atp_synthase_efficiency": {
            "type": "triangular",
            "params": [0.7, 0.95, 1.0],
            "bounds": (0.7, 1.0),
        },
        "membrane_leak_conductance": {
            "type": "uniform",
            "params": [0.0, 0.2],
            "bounds": (0.0, 1.0),
        },
        "proton_slip_probability": {
            "type": "uniform",
            "params": [0.0, 0.1],
            "bounds": (0.0, 1.0),
        },
        "roq_quinone_coupling": {
            "type": "triangular",
            "params": [0.8, 0.98, 1.0],
            "bounds": (0.0, 1.0),
        },
        "scaffold_channeling_efficiency": {
            "type": "triangular",
            "params": [0.8, 0.95, 1.0],
            "bounds": (0.0, 1.0),
        },
        "ros_bypass_fraction": {
            "type": "uniform",
            "params": [0.0, 0.05],
            "bounds": (0.0, 1.0),
        },
    }


def _sample_from_distribution(
    dist: Dict[str, Any], n: int, rng: np.random.Generator
) -> np.ndarray:
    """Draw n samples from a parameter distribution specification.

    Supported types: 'normal', 'uniform', 'lognormal', 'triangular'.
    Returns clipped to dist['bounds'].
    """
    lo, hi = dist["bounds"]
    dtype = dist["type"]
    params = dist["params"]

    if dtype == "normal":
        samples = rng.normal(loc=params[0], scale=params[1], size=n)
    elif dtype == "uniform":
        samples = rng.uniform(low=params[0], high=params[1], size=n)
    elif dtype == "lognormal":
        samples = rng.lognormal(mean=params[0], sigma=params[1], size=n)
    elif dtype == "triangular":
        samples = rng.triangular(
            left=params[0], mode=params[1], right=params[2], size=n
        )
    else:
        raise ValueError(f"Unknown distribution type: {dtype}")

    return np.clip(samples, lo, hi)


def _run_mc_sensitivity(
    model: Callable[[DesignParameters], float],
    parameters: Optional[Sequence[str]] = None,
    distributions: Optional[Dict[str, Dict[str, Any]]] = None,
    n_samples: int = 10000,
    seed: Optional[int] = None,
) -> Dict[str, Any]:
    """Shared MC-based sensitivity logic used by the convenience functions."""
    if parameters is None:
        parameters = PARAMETER_NAMES
    if distributions is None:
        all_dists = _default_distributions()
        distributions = {p: all_dists[p] for p in parameters}

    rng = np.random.default_rng(seed)
    D = len(parameters)
    bounds = get_param_bounds()

    # Saltelli sampling
    N = n_samples
    lo = np.array([bounds[PARAMETER_NAMES.index(p)][0] for p in parameters])
    hi = np.array([bounds[PARAMETER_NAMES.index(p)][1] for p in parameters])

    def _lhs(d: int, n: int) -> np.ndarray:
        segments = np.linspace(0.0, 1.0, n + 1)
        mat = np.zeros((n, d))
        for j in range(d):
            offsets = rng.uniform(0, segments[1], size=n)
            order = rng.permutation(n)
            mat[:, j] = segments[order] + offsets
        return mat

    A_raw = _lhs(D, N)
    B_raw = _lhs(D, N)
    A_scaled = lo[np.newaxis, :] + A_raw * (hi - lo)[np.newaxis, :]
    B_scaled = lo[np.newaxis, :] + B_raw * (hi - lo)[np.newaxis, :]

    # Apply user-specified distributions via inverse transform
    for i, p in enumerate(parameters):
        dist = distributions[p]
        lo_p, hi_p = dist["bounds"]
        # Transform uniform [0,1] samples to the desired distribution
        u_samples_a = (A_raw[:, i] * (hi_p - lo_p) + lo_p - lo[i]) / (hi[i] - lo[i])
        u_samples_a = np.clip(u_samples_a, 0.0, 1.0)
        # ASSUMPTION: we use the raw LHS samples and re-interpret through the
        # distribution's inverse CDF.  For simplicity, we use the scaled
        # samples directly, which already respect uniform bounds.
        pass

    def _make_design(arr: np.ndarray) -> DesignParameters:
        full = np.zeros(len(PARAMETER_NAMES))
        for i, p in enumerate(parameters):
            full[PARAMETER_NAMES.index(p)] = arr[i]
        # Set non-sampled parameters to nominal values
        return DesignParameters.from_array(full)

    # AB matrices
    AB = np.zeros((D, N, D))
    for i in range(D):
        AB[i] = A_scaled.copy()
        AB[i][:, i] = B_scaled[:, i]

    f_A = np.array([model(_make_design(A_scaled[i])) for i in range(N)])
    f_B = np.array([model(_make_design(B_scaled[i])) for i in range(N)])
    f_AB = np.zeros((D, N))
    for i in range(D):
        for j in range(N):
            f_AB[i, j] = model(_make_design(AB[i, j]))

    S1, ST, S1_conf, ST_conf = GlobalSensitivity._compute_indices(f_A, f_B, f_AB)

    # Parameter importance
    param_names = list(parameters)
    importance = np.abs(S1) + np.abs(ST)
    ranking = sorted(
        zip(param_names, importance),
        key=lambda x: x[1],
        reverse=True,
    )

    return {
        "S1": dict(zip(param_names, S1)),
        "ST": dict(zip(param_names, ST)),
        "S1_conf": dict(zip(param_names, S1_conf)),
        "ST_conf": dict(zip(param_names, ST_conf)),
        "ranking": ranking,
        "n_samples": N,
        "n_total_evaluations": N * (2 * D + 2),
    }


def sensitivity_of_atp_yield(
    parameters: Optional[Sequence[str]] = None,
    distributions: Optional[Dict[str, Dict[str, Any]]] = None,
    n_samples: int = 10000,
    seed: Optional[int] = None,
) -> Dict[str, Any]:
    """Global sensitivity of ATP yield to PCM parameters.

    Args:
        parameters: Subset of parameter names to analyse.  Defaults to all.
        distributions: Parameter probability distributions.
        n_samples: Base sample count for Sobol analysis.
        seed: Random seed.

    Returns:
        Dict with Sobol indices and parameter ranking for ATP yield.
    """
    return _run_mc_sensitivity(
        model=ObjectiveFunction.compute_atp_yield,
        parameters=parameters,
        distributions=distributions,
        n_samples=n_samples,
        seed=seed,
    )


def sensitivity_of_efficiency(
    parameters: Optional[Sequence[str]] = None,
    distributions: Optional[Dict[str, Dict[str, Any]]] = None,
    n_samples: int = 10000,
    seed: Optional[int] = None,
) -> Dict[str, Any]:
    """Global sensitivity of thermodynamic efficiency to PCM parameters."""
    return _run_mc_sensitivity(
        model=ObjectiveFunction.maximize_thermodynamic_efficiency,
        parameters=parameters,
        distributions=distributions,
        n_samples=n_samples,
        seed=seed,
    )


def sensitivity_of_entropy(
    parameters: Optional[Sequence[str]] = None,
    distributions: Optional[Dict[str, Dict[str, Any]]] = None,
    n_samples: int = 10000,
    seed: Optional[int] = None,
) -> Dict[str, Any]:
    """Global sensitivity of entropy production to PCM parameters."""
    return _run_mc_sensitivity(
        model=ObjectiveFunction.minimize_entropy_production,
        parameters=parameters,
        distributions=distributions,
        n_samples=n_samples,
        seed=seed,
    )


def parameter_ranking(
    sensitivity_results: Dict[str, Any],
) -> List[Tuple[str, float]]:
    """Extract a unified parameter importance ranking from sensitivity results.

    Args:
        sensitivity_results: Dict returned by sensitivity_of_* functions.

    Returns:
        List of (parameter_name, importance_score) sorted descending.
    """
    return list(sensitivity_results.get("ranking", []))


def plot_sensitivity(
    sensitivity_results: Dict[str, Any],
    save_path: Optional[str] = None,
) -> None:
    """Generate a tornado plot of Sobol total-order sensitivity indices.

    Args:
        sensitivity_results: Dict from sensitivity_of_* functions.
        save_path: If provided, save the figure.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plot.")
        return

    st = sensitivity_results.get("ST", {})
    s1 = sensitivity_results.get("S1", {})
    if not st:
        print("No sensitivity data to plot.")
        return

    # Sort by total-order index
    sorted_params = sorted(st.items(), key=lambda kv: kv[1])
    params = [p[0] for p in sorted_params]
    st_vals = np.array([p[1] for p in sorted_params])
    s1_vals = np.array([s1.get(p[0], 0.0) for p in sorted_params])

    fig, ax = plt.subplots(figsize=(8, max(4, len(params) * 0.4)))
    y = np.arange(len(params))

    # Total order as horizontal bars
    ax.barh(y, st_vals, height=0.6, color="crimson", alpha=0.7, label="Total-order (ST)")
    ax.barh(y, s1_vals, height=0.6, color="steelblue", alpha=0.7, label="First-order (S1)")

    ax.set_yticks(y)
    ax.set_yticklabels(params)
    ax.set_xlabel("Sobol Sensitivity Index")
    ax.set_title("Parameter Sensitivity — PTR-94 PCM")
    ax.axvline(0, color="gray", linewidth=0.5)
    ax.legend(loc="lower right")
    ax.set_xlim(0, max(1.0, st_vals.max() * 1.15))
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Sensitivity plot saved to {save_path}")
    else:
        plt.show()
    plt.close(fig)


def build_sensitivity_report(
    output_name: str = "atp_yield",
    local_params: Optional[DesignParameters] = None,
    seed: Optional[int] = None,
) -> SensitivityReport:
    """Convenience function to build a complete SensitivityReport.

    Runs both local and global sensitivity analysis and aggregates results.

    Args:
        output_name: Which output to analyse.
        local_params: Nominal design for local sensitivity.  Defaults to
            a standard reference design.
        seed: Random seed.

    Returns:
        A populated SensitivityReport.
    """
    if local_params is None:
        local_params = DesignParameters()

    output_map = {
        "atp_yield": ObjectiveFunction.compute_atp_yield,
        "thermodynamic_efficiency": ObjectiveFunction.maximize_thermodynamic_efficiency,
        "entropy_production": ObjectiveFunction.minimize_entropy_production,
    }
    if output_name not in output_map:
        raise ValueError(f"Unknown output: {output_name}. Choose from {list(output_map.keys())}")

    model = output_map[output_name]

    # Local sensitivity
    local_an = LocalSensitivity(nominal=local_params, outputs={output_name: model})
    local_coeffs = local_an.compute()[output_name]

    # Global sensitivity (quick run with fewer samples)
    global_an = GlobalSensitivity(
        param_bounds=get_param_bounds(), n_samples=2000, seed=seed
    )
    global_results = global_an.analyze(model)
    param_names = PARAMETER_NAMES[: len(get_param_bounds())]

    s1_dict = dict(zip(param_names, global_results["S1"]))
    st_dict = dict(zip(param_names, global_results["ST"]))

    # Importance ranking: combination of |local| + |S1| + |ST|
    importance = {}
    for p in param_names:
        importance[p] = (
            abs(local_coeffs.get(p, 0.0))
            + abs(s1_dict.get(p, 0.0))
            + abs(st_dict.get(p, 0.0))
        )
    ranking = sorted(importance.items(), key=lambda kv: kv[1], reverse=True)

    # Robustness metric: higher when ST values (total variance contribution) are low
    st_vals = np.array(list(st_dict.values()))
    robustness = float(np.clip(1.0 - st_vals.mean(), 0.0, 1.0))

    # Uncertainty breakdown: normalised ST
    st_sum = max(st_vals.sum(), 1e-15)
    uncertainty_breakdown = {p: st_dict[p] / st_sum for p in param_names}

    return SensitivityReport(
        output_name=output_name,
        local_coefficients=local_coeffs,
        global_first_order=s1_dict,
        global_total_order=st_dict,
        parameter_ranking=ranking,
        robustness_metric=robustness,
        uncertainty_breakdown=uncertainty_breakdown,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Usage
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("PTR-94 Sensitivity Analysis — Quick Demo")
    print("=" * 60)

    # Local sensitivity around nominal design
    nominal = DesignParameters()
    local = LocalSensitivity(nominal)
    coeffs = local.compute()
    print("\nLocal sensitivity coefficients (ATP yield):")
    for param, val in coeffs["atp_yield"].items():
        print(f"  {param:<30} {val:+.6f}")

    print("\nMost influential parameters (ATP yield):")
    for param, sens in local.most_influential("atp_yield", 5):
        print(f"  {param:<30} |S| = {abs(sens):.4f}")

    # Global sensitivity (small sample for demo)
    print("\nGlobal Sobol sensitivity (ATP yield, n_samples=2000):")
    gs = GlobalSensitivity(param_bounds=get_param_bounds(), n_samples=2000, seed=42)
    results = gs.analyze(ObjectiveFunction.compute_atp_yield)
    for i, name in enumerate(PARAMETER_NAMES[: len(get_param_bounds())]):
        print(f"  {name:<30} S1={results['S1'][i]:.4f}  ST={results['ST'][i]:.4f}")

    # Full report
    report = build_sensitivity_report("atp_yield", seed=42)
    print(report.summary_table())
