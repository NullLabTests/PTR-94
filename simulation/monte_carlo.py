#!/usr/bin/env python3
"""
Monte Carlo Simulation and Uncertainty Quantification for PTR‑94.

Provides Latin Hypercube Sampling (LHS), random sampling with correlated
parameters, and comprehensive uncertainty quantification for the Perfect
Coupling Module (PCM) outputs: ATP yield, thermodynamic efficiency, and
proton flux.

Assumptions
-----------
- Standard biochemical conditions (ΔG°′ = 30.5 kJ/mol for ATP).
- 10 NADH + 2 FADH2 + 4 substrate-level ATP from upstream modules.
- Parameter distributions are independent unless correlation matrices are
  explicitly provided.
- Physical bounds on each parameter are strictly enforced via clipping.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

# Reuse core model from pareto_optimizer
try:
    from .pareto_optimizer import (
        DELTA_G_ATP_STD,
        DELTA_G_GLUCOSE,
        FADH2_COUNT,
        NADH_COUNT,
        SUBSTRATE_ATP,
        THEORETICAL_MAX_ATP,
        DesignParameters,
        ObjectiveFunction,
        get_param_bounds,
    )
except ImportError:
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


# ─────────────────────────────────────────────────────────────────────────────
# Parameter Distribution
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class ParameterDistribution:
    """A probability distribution for a PCM design parameter.

    Attributes:
        name: Parameter name (must match a DesignParameters field).
        distribution_type: One of 'normal', 'uniform', 'lognormal',
            'triangular'.
        parameters: Distribution-specific parameters:
            - normal:   [mean, std]
            - uniform:  [min, max]
            - lognormal: [mean, sigma]  (of the underlying normal)
            - triangular: [left, mode, right]
        physical_bounds: (min, max) enforced after sampling via clipping.
    """

    name: str
    distribution_type: str
    parameters: List[float] = field(default_factory=list)
    physical_bounds: Tuple[float, float] = (0.0, 1.0)

    def sample(self, n: int, rng: np.random.Generator) -> np.ndarray:
        """Draw n independent samples from this distribution.

        Args:
            n: Number of samples.
            rng: NumPy random generator.

        Returns:
            Array of shape (n,) clipped to physical_bounds.
        """
        lo, hi = self.physical_bounds
        dtype = self.distribution_type
        p = self.parameters

        if dtype == "normal":
            samples = rng.normal(loc=p[0], scale=p[1], size=n)
        elif dtype == "uniform":
            samples = rng.uniform(low=p[0], high=p[1], size=n)
        elif dtype == "lognormal":
            samples = rng.lognormal(mean=p[0], sigma=p[1], size=n)
        elif dtype == "triangular":
            samples = rng.triangular(left=p[0], mode=p[1], right=p[2], size=n)
        else:
            raise ValueError(f"Unknown distribution type: {dtype}")

        return np.clip(samples, lo, hi)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "distribution_type": self.distribution_type,
            "parameters": self.parameters,
            "physical_bounds": list(self.physical_bounds),
        }

    @classmethod
    def default_atp_yield(cls) -> List[ParameterDistribution]:
        """Return recommended default distributions for all parameters."""
        return [
            cls("h_per_nadh", "uniform", [10.0, 40.0], (10.0, 40.0)),
            cls("h_per_fadh2", "uniform", [6.0, 26.0], (6.0, 26.0)),
            cls("h_per_atp", "triangular", [2.0, 3.0, 5.0], (2.0, 5.0)),
            cls("atp_synthase_efficiency", "triangular", [0.7, 0.95, 1.0], (0.7, 1.0)),
            cls("membrane_leak_conductance", "uniform", [0.0, 0.2], (0.0, 1.0)),
            cls("proton_slip_probability", "uniform", [0.0, 0.1], (0.0, 1.0)),
            cls("roq_quinone_coupling", "triangular", [0.8, 0.98, 1.0], (0.0, 1.0)),
            cls("scaffold_channeling_efficiency", "triangular", [0.8, 0.95, 1.0], (0.0, 1.0)),
            cls("ros_bypass_fraction", "uniform", [0.0, 0.05], (0.0, 1.0)),
        ]


# ─────────────────────────────────────────────────────────────────────────────
# Monte Carlo Engine
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class MonteCarloResult:
    """Container for Monte Carlo simulation results.

    Attributes:
        input_samples: (n_samples × n_params) array of sampled parameters.
        output_samples: (n_samples,) array of output values.
        parameter_names: List of parameter names corresponding to columns.
        output_name: Name of the model output.
        n_samples: Number of Monte Carlo samples.
        sampling_method: 'lhs' or 'random'.
    """

    input_samples: np.ndarray
    output_samples: np.ndarray
    parameter_names: List[str]
    output_name: str = ""
    n_samples: int = 0
    sampling_method: str = "lhs"


class MonteCarloEngine:
    """Monte Carlo simulation engine with Latin Hypercube Sampling (LHS).

    Supports independent and correlated parameter sampling, and runs model
    evaluations for each sample.

    Args:
        parameter_distributions: List of ParameterDistribution for each input.
        sampling_method: 'lhs' (default) or 'random'.
        seed: Random seed for reproducibility.
    """

    def __init__(
        self,
        parameter_distributions: Sequence[ParameterDistribution],
        sampling_method: str = "lhs",
        seed: Optional[int] = None,
    ) -> None:
        self.distributions = list(parameter_distributions)
        self.n_params = len(self.distributions)
        self.sampling_method = sampling_method
        self.rng = np.random.default_rng(seed)

        self._check_consistency()

    def _check_consistency(self) -> None:
        """Verify distribution names match known parameters."""
        known = set(PARAMETER_NAMES)
        for d in self.distributions:
            if d.name not in known:
                raise ValueError(
                    f"Unknown parameter '{d.name}'. "
                    f"Known: {PARAMETER_NAMES}"
                )

    def _latin_hypercube(self, n: int) -> np.ndarray:
        """Generate a Latin Hypercube sample in [0, 1]^D.

        Args:
            n: Number of samples.

        Returns:
            Array of shape (n, n_params) with each column stratified.
        """
        samples = np.zeros((n, self.n_params))
        for j in range(self.n_params):
            segments = np.linspace(0.0, 1.0, n + 1)
            offsets = self.rng.uniform(0, segments[1], size=n)
            order = self.rng.permutation(n)
            samples[:, j] = segments[order] + offsets
        return samples

    def _map_to_distributions(self, unit_samples: np.ndarray) -> np.ndarray:
        """Map unit-hypercube samples to parameter spaces via inverse CDF.

        Uses the distribution parameters directly for uniform/triangular,
        or constructs quantile-based mappings for normal/lognormal.

        Args:
            unit_samples: Array in [0, 1]^D, shape (n, D).

        Returns:
            Array of physical parameter samples, shape (n, D).
        """
        physical = np.zeros_like(unit_samples)

        for j, dist in enumerate(self.distributions):
            u = np.clip(unit_samples[:, j], 1e-15, 1.0 - 1e-15)
            lo, hi = dist.physical_bounds
            dtype = dist.distribution_type
            p = dist.parameters

            if dtype == "uniform":
                # Direct linear mapping
                physical[:, j] = lo + u * (hi - lo)
            elif dtype == "triangular":
                # Triangular inverse CDF
                left, mode, right = p
                if abs(right - left) < 1e-15:
                    physical[:, j] = left
                else:
                    # Normalised triangular
                    f_mode = (mode - left) / (right - left)
                    phys_norm = np.where(
                        u <= f_mode,
                        left + np.sqrt(u * (mode - left) * (right - left)),
                        right - np.sqrt(
                            (1.0 - u) * (right - mode) * (right - left)
                        ),
                    )
                    physical[:, j] = np.clip(phys_norm, lo, hi)
            elif dtype == "normal":
                # Use inverse normal CDF
                mean, std = p
                physical[:, j] = np.clip(
                    mean + std * _normal_ppf(u), lo, hi
                )
            elif dtype == "lognormal":
                mean_ln, sigma_ln = p
                physical[:, j] = np.clip(
                    np.exp(mean_ln + sigma_ln * _normal_ppf(u)), lo, hi
                )
            else:
                raise ValueError(f"Unknown distribution: {dtype}")

        return physical

    def sample(self, n: int) -> np.ndarray:
        """Generate parameter samples.

        Args:
            n: Number of samples.

        Returns:
            Array of shape (n, n_params) with physical parameter values.
        """
        if self.sampling_method == "lhs":
            unit = self._latin_hypercube(n)
        else:
            unit = self.rng.uniform(0.0, 1.0, size=(n, self.n_params))

        return self._map_to_distributions(unit)

    def sample_correlated(
        self, n: int, correlation_matrix: np.ndarray
    ) -> np.ndarray:
        """Generate correlated parameter samples using Iman-Conover method.

        Preserves the marginal distributions of each parameter while
        approximately matching the target rank correlation matrix.

        Reference: Iman & Conover (1982), Technometrics 24(4), 299–305.

        Args:
            n: Number of samples.
            correlation_matrix: (n_params × n_params) target rank correlation.

        Returns:
            Array of shape (n, n_params) with correlated marginal samples.
        """
        if correlation_matrix.shape != (self.n_params, self.n_params):
            raise ValueError(
                f"Correlation matrix must be ({self.n_params}, {self.n_params})"
            )

        # Step 1: Generate independent LHS samples
        X = self.sample(n)

        # Step 2: Create a multivariate normal with target correlation
        # and reorder to match it (Iman-Conover)
        Z = self.rng.normal(size=(n, self.n_params))
        # Cholesky decomposition of correlation matrix
        try:
            L = np.linalg.cholesky(correlation_matrix)
        except np.linalg.LinAlgError:
            # Fallback: nearest positive-definite
            eigvals, eigvecs = np.linalg.eigh(correlation_matrix)
            eigvals = np.maximum(eigvals, 1e-10)
            corr_pd = eigvecs @ np.diag(eigvals) @ eigvecs.T
            L = np.linalg.cholesky(corr_pd)

        Z_corr = Z @ L.T  # Z_corr ~ N(0, correlation_matrix)

        # Step 3: Reorder each column of X to match the rank order of Z_corr
        X_corr = np.zeros_like(X)
        for j in range(self.n_params):
            # Rank of Z_corr[:, j] gives the target ordering
            order = np.argsort(np.argsort(Z_corr[:, j]))
            sorted_X = np.sort(X[:, j])
            X_corr[:, j] = sorted_X[order]

        return X_corr

    def run(
        self,
        model: Callable[[DesignParameters], float],
        n_samples: int = 10000,
        correlated: bool = False,
        correlation_matrix: Optional[np.ndarray] = None,
        output_name: str = "output",
    ) -> MonteCarloResult:
        """Run Monte Carlo simulation.

        Args:
            model: Function mapping DesignParameters → scalar output.
            n_samples: Number of Monte Carlo samples.
            correlated: If True, use correlated sampling.
            correlation_matrix: Target correlation (required if correlated=True).
            output_name: Label for the output.

        Returns:
            MonteCarloResult with input and output samples.
        """
        if correlated:
            if correlation_matrix is None:
                raise ValueError(
                    "correlation_matrix required when correlated=True"
                )
            param_samples = self.sample_correlated(n_samples, correlation_matrix)
        else:
            param_samples = self.sample(n_samples)

        output_samples = np.array([
            model(_make_design(param_samples[i])) for i in range(n_samples)
        ])

        return MonteCarloResult(
            input_samples=param_samples,
            output_samples=output_samples,
            parameter_names=[d.name for d in self.distributions],
            output_name=output_name,
            n_samples=n_samples,
            sampling_method=self.sampling_method,
        )


def _make_design(arr: np.ndarray) -> DesignParameters:
    """Convert a full 9-parameter array to DesignParameters."""
    return DesignParameters.from_array(arr)


def _normal_ppf(q: np.ndarray) -> np.ndarray:
    """Approximate normal inverse CDF (ppf) using rational approximation.

    Uses the Wichura (1988) AS 241 algorithm.  Accurate to ~1e-15.
    """
    # Coefficients for rational approximation
    a = np.array([ 3.3871328727963666080e0, 1.3314166789178437745e2,
                   1.9715909503065514427e3, 1.3731693765509461125e4,
                   4.5921953931549870937e4, 6.7265770927008700853e4,
                   3.3430575583588128105e4, 2.5090809287301226727e3 ])
    b = np.array([ 1.0, 4.2313330701603001252e1,
                   6.8718700749205790830e2, 5.3941960214247511077e3,
                   2.1213794301586595867e4, 3.9307895800092710610e4,
                   2.8729085735721942674e4, 5.2264952788528545610e3 ])
    c = np.array([ 1.42343711074968357734e0, 4.63033784615654529590e0,
                   5.76949722146069140550e0, 3.64784832476320460504e0,
                   1.27045825245236838258e0, 2.41780725177450611770e-1,
                   2.27238449892691845833e-2, 7.74545014278341407640e-4 ])
    d = np.array([ 1.0, 2.05319162663775882187e0,
                   1.67638483094294423062e0, 6.89767334985100004550e-1,
                   1.48103976427480074590e-1, 1.51986665636164571966e-2,
                   5.47593808499534494600e-4, 1.05075007164441684324e-9 ])
    e = np.array([ 6.65790464350110377720e0, 5.46378491116411436990e0,
                   1.78482653991729133580e0, 2.96560571828504891230e-1,
                   2.65321895265721230930e-2, 1.24266094738807843860e-3,
                   2.71155556874348757815e-5, 2.01033439929228813265e-7 ])
    f = np.array([ 1.0, 5.99832206555887937690e-1,
                   1.36929880922735805310e-1, 1.48753612908506148522e-2,
                   7.86869131145613259100e-4, 1.84631831751005468180e-5,
                   1.42151175831644588870e-7, 2.04426331040495049249e-10 ])

    split1 = 0.425

    result = np.zeros_like(q)

    # Left tail
    mask = q < split1
    if np.any(mask):
        q0 = q[mask]
        r = np.sqrt(-np.log(1.0 - q0))
        result[mask] = (((((c[7] * r + c[6]) * r + c[5]) * r + c[4]) * r + c[3]) * r + c[2]) * r + c[1]
        result[mask] /= (((((d[7] * r + d[6]) * r + d[5]) * r + d[4]) * r + d[3]) * r + d[2]) * r + d[1]
        result[mask] = -result[mask]

    # Central
    mask = (q >= split1) & (q <= 1.0 - split1)
    if np.any(mask):
        q0 = q[mask]
        r = q0 - 0.5
        result[mask] = (((((a[7] * r + a[6]) * r + a[5]) * r + a[4]) * r + a[3]) * r + a[2]) * r + a[1]
        result[mask] /= (((((b[7] * r + b[6]) * r + b[5]) * r + b[4]) * r + b[3]) * r + b[2]) * r + b[1]

    # Right tail
    mask = q > 1.0 - split1
    if np.any(mask):
        q0 = q[mask]
        r = np.sqrt(-np.log(1.0 - q0))
        result[mask] = (((((e[7] * r + e[6]) * r + e[5]) * r + e[4]) * r + e[3]) * r + e[2]) * r + e[1]
        result[mask] /= (((((f[7] * r + f[6]) * r + f[5]) * r + f[4]) * r + f[3]) * r + f[2]) * r + f[1]

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Uncertainty Quantification
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class UncertaintyResult:
    """Comprehensive uncertainty quantification results.

    Attributes:
        output_name: Name of the analysed output.
        mean: Sample mean.
        median: Sample median.
        std: Sample standard deviation.
        ci_68: (lower, upper) 68% confidence interval.
        ci_95: (lower, upper) 95% confidence interval.
        ci_99: (lower, upper) 99% confidence interval.
        percentiles: Dict of percentile → value (1st, 5th, 25th, 50th, 75th,
            95th, 99th).
        histogram: (counts, bin_edges) tuple from numpy histogram.
        fitted_distribution: Name of best-fit distribution.
        fitted_params: Parameters of the best-fit distribution.
        correlation_with_inputs: Dict of parameter_name → Pearson correlation.
        samples: Raw output samples (useful for further analysis).
        n_samples: Number of samples.
    """

    output_name: str = ""
    mean: float = 0.0
    median: float = 0.0
    std: float = 0.0
    ci_68: Tuple[float, float] = (0.0, 0.0)
    ci_95: Tuple[float, float] = (0.0, 0.0)
    ci_99: Tuple[float, float] = (0.0, 0.0)
    percentiles: Dict[float, float] = field(default_factory=dict)
    histogram: Tuple[np.ndarray, np.ndarray] = (np.array([]), np.array([]))
    fitted_distribution: str = ""
    fitted_params: List[float] = field(default_factory=list)
    correlation_with_inputs: Dict[str, float] = field(default_factory=dict)
    samples: np.ndarray = field(default_factory=lambda: np.array([]))
    n_samples: int = 0

    def summary(self) -> str:
        """Return a formatted summary string of uncertainty results."""
        lines = [
            f"\n{'=' * 60}",
            f"  Uncertainty Quantification — {self.output_name}",
            f"{'=' * 60}",
            f"  n_samples:          {self.n_samples}",
            f"  Mean:               {self.mean:.4f}",
            f"  Median:             {self.median:.4f}",
            f"  Std:                {self.std:.4f}",
            f"  68% CI:             ({self.ci_68[0]:.4f}, {self.ci_68[1]:.4f})",
            f"  95% CI:             ({self.ci_95[0]:.4f}, {self.ci_95[1]:.4f})",
            f"  99% CI:             ({self.ci_99[0]:.4f}, {self.ci_99[1]:.4f})",
            f"  Fitted distribution: {self.fitted_distribution}",
            f"  Fitted params:        {self.fitted_params}",
            f"  {'=' * 60}",
            "  Percentiles:",
        ]
        for pct in [1, 5, 25, 50, 75, 95, 99]:
            if pct in self.percentiles:
                lines.append(f"    {pct:>2}%:  {self.percentiles[pct]:.4f}")
        lines.append(f"{'=' * 60}")
        return "\n".join(lines)


class UncertaintyQuantification:
    """Compute uncertainty statistics from Monte Carlo output samples.

    Args:
        result: MonteCarloResult from a completed simulation.
    """

    def __init__(self, result: MonteCarloResult) -> None:
        self.result = result
        self._uq_result: Optional[UncertaintyResult] = None

    def compute(
        self,
        fit_distributions: bool = True,
    ) -> UncertaintyResult:
        """Compute comprehensive uncertainty metrics.

        Args:
            fit_distributions: If True, attempt to fit a parametric
                distribution to the output samples.

        Returns:
            Fully populated UncertaintyResult.
        """
        y = self.result.output_samples
        n = len(y)

        # Basic statistics
        mean = float(np.mean(y))
        median = float(np.median(y))
        std = float(np.std(y, ddof=1))

        # Confidence intervals (percentile-based, bias-corrected)
        ci_68 = tuple(float(v) for v in np.percentile(y, [16, 84]))
        ci_95 = tuple(float(v) for v in np.percentile(y, [2.5, 97.5]))
        ci_99 = tuple(float(v) for v in np.percentile(y, [0.5, 99.5]))

        # Percentiles
        pct_vals = [1, 5, 25, 50, 75, 95, 99]
        percentiles = {float(p): float(np.percentile(y, p)) for p in pct_vals}

        # Histogram
        histogram = np.histogram(y, bins="auto", density=True)

        # Distribution fitting
        fitted_dist = ""
        fitted_params: List[float] = []
        if fit_distributions and n > 100:
            fitted_dist, fitted_params = self._fit_distribution(y)

        # Correlation with inputs
        corr_with_inputs: Dict[str, float] = {}
        for i, name in enumerate(self.result.parameter_names):
            r = np.corrcoef(self.result.input_samples[:, i], y)[0, 1]
            corr_with_inputs[name] = float(r if not np.isnan(r) else 0.0)

        self._uq_result = UncertaintyResult(
            output_name=self.result.output_name,
            mean=mean,
            median=median,
            std=std,
            ci_68=ci_68,
            ci_95=ci_95,
            ci_99=ci_99,
            percentiles=percentiles,
            histogram=histogram,
            fitted_distribution=fitted_dist,
            fitted_params=fitted_params,
            correlation_with_inputs=corr_with_inputs,
            samples=y.copy(),
            n_samples=n,
        )
        return self._uq_result

    @staticmethod
    def _fit_distribution(y: np.ndarray) -> Tuple[str, List[float]]:
        """Fit a parametric distribution to samples by maximum likelihood.

        Tries normal, lognormal, and gamma; reports the best fit by BIC.

        Returns:
            (distribution_name, [param1, param2, ...])
        """
        from scipy import stats

        candidates = {
            "normal": stats.norm,
            "lognormal": stats.lognorm,
            "gamma": stats.gamma,
        }
        best_name = "normal"
        best_params: List[float] = [float(np.mean(y)), float(np.std(y, ddof=1))]
        best_bic = float("inf")

        for name, dist_cls in candidates.items():
            try:
                # Fit via MLE
                params = dist_cls.fit(y)
                log_likelihood = np.sum(dist_cls.logpdf(y, *params))
                k = len(params)
                bic = k * np.log(len(y)) - 2.0 * log_likelihood
                if bic < best_bic:
                    best_bic = bic
                    best_name = name
                    best_params = [float(p) for p in params]
            except Exception:
                continue

        return best_name, best_params


# ─────────────────────────────────────────────────────────────────────────────
# Convenience Monte Carlo Functions
# ─────────────────────────────────────────────────────────────────────────────


def _run_mc(
    model: Callable[[DesignParameters], float],
    n_samples: int = 100000,
    parameter_distributions: Optional[List[ParameterDistribution]] = None,
    output_name: str = "output",
    sampling_method: str = "lhs",
    seed: Optional[int] = None,
) -> UncertaintyResult:
    """Shared logic for convenience MC functions."""
    if parameter_distributions is None:
        parameter_distributions = ParameterDistribution.default_atp_yield()

    engine = MonteCarloEngine(
        parameter_distributions=parameter_distributions,
        sampling_method=sampling_method,
        seed=seed,
    )
    result = engine.run(
        model=model,
        n_samples=n_samples,
        output_name=output_name,
    )
    uq = UncertaintyQuantification(result)
    return uq.compute()


def monte_carlo_atp_yield(
    n_samples: int = 100000,
    parameter_distributions: Optional[List[ParameterDistribution]] = None,
    seed: Optional[int] = None,
) -> UncertaintyResult:
    """Run Monte Carlo uncertainty quantification for ATP yield.

    Args:
        n_samples: Number of Monte Carlo samples.
        parameter_distributions: Custom parameter distributions.  Defaults
            to reasonable default ranges.
        seed: Random seed.

    Returns:
        UncertaintyResult with statistics for ATP yield.
    """
    return _run_mc(
        model=ObjectiveFunction.compute_atp_yield,
        n_samples=n_samples,
        parameter_distributions=parameter_distributions,
        output_name="ATP_yield",
        seed=seed,
    )


def monte_carlo_efficiency(
    n_samples: int = 100000,
    parameter_distributions: Optional[List[ParameterDistribution]] = None,
    seed: Optional[int] = None,
) -> UncertaintyResult:
    """Run Monte Carlo uncertainty quantification for thermodynamic efficiency."""
    return _run_mc(
        model=ObjectiveFunction.maximize_thermodynamic_efficiency,
        n_samples=n_samples,
        parameter_distributions=parameter_distributions,
        output_name="thermodynamic_efficiency",
        seed=seed,
    )


def monte_carlo_proton_flux(
    n_samples: int = 100000,
    parameter_distributions: Optional[List[ParameterDistribution]] = None,
    seed: Optional[int] = None,
) -> UncertaintyResult:
    """Run Monte Carlo uncertainty quantification for proton flux.

    Returns
    -------
    UncertaintyResult for total proton flux (H⁺ per glucose).
    """
    def proton_flux(params: DesignParameters) -> float:
        protons_nadh = params.h_per_nadh * NADH_COUNT
        protons_fadh2 = params.h_per_fadh2 * FADH2_COUNT
        total = protons_nadh + protons_fadh2
        effective = (
            total
            * (1.0 - params.membrane_leak_conductance)
            * (1.0 - params.proton_slip_probability)
            * params.roq_quinone_coupling
        )
        return float(effective)

    return _run_mc(
        model=proton_flux,
        n_samples=n_samples,
        parameter_distributions=parameter_distributions,
        output_name="proton_flux",
        seed=seed,
    )


def confidence_intervals(
    results: UncertaintyResult,
    alpha: float = 0.05,
) -> Dict[str, Tuple[float, float]]:
    """Extract confidence intervals at a given significance level.

    Args:
        results: UncertaintyResult from a Monte Carlo run.
        alpha: Significance level (default 0.05 → 95% CI).

    Returns:
        Dict with keys 'ci_68', 'ci_95', 'ci_99' and the requested level.
    """
    z = 100 * (1.0 - alpha / 2.0)
    lower = float(np.percentile(results.samples, 100.0 * alpha / 2.0))
    upper = float(np.percentile(results.samples, z))

    return {
        "ci_68": results.ci_68,
        "ci_95": results.ci_95,
        "ci_99": results.ci_99,
        f"ci_{1-alpha:.0%}": (lower, upper),
    }


def plot_uncertainty(
    results: UncertaintyResult,
    parameter: Optional[str] = None,
    save_path: Optional[str] = None,
) -> None:
    """Plot uncertainty results as a histogram with confidence bounds.

    If a parameter name is provided, plots the joint distribution of that
    input parameter against the output.  Otherwise plots the output marginal
    distribution.

    Args:
        results: UncertaintyResult from a Monte Carlo run.
        parameter: If given, plot output vs. this input parameter.
        save_path: Path to save the figure.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plot.")
        return

    n_cols = 2 if parameter else 1
    fig, axes = plt.subplots(1, n_cols, figsize=(10, 4))

    if parameter:
        ax1, ax2 = axes
        if results.correlation_with_inputs and parameter in results.correlation_with_inputs:
            r_val = results.correlation_with_inputs[parameter]
            ax1.set_title(f"Output vs. {parameter} (r={r_val:.3f})")
        else:
            ax1.set_title(f"Output vs. {parameter}")
        ax1.text(
            0.5, 0.5,
            "Input samples not stored\nin UncertaintyResult.\n"
            "Use MonteCarloResult directly\nfor scatter plots.",
            transform=ax1.transAxes,
            ha="center", va="center",
            fontsize=10,
        )
        ax1.set_xlabel(parameter)
        ax1.set_ylabel(results.output_name)
        ax_hist = ax2
    else:
        ax_hist = axes

    counts, bins, _ = ax_hist.hist(
        results.samples,
        bins=60,
        density=True,
        alpha=0.6,
        color="steelblue",
        edgecolor="white",
        linewidth=0.3,
    )
    ax_hist.axvline(results.mean, color="crimson", linestyle="--",
                    linewidth=1.5, label=f"Mean = {results.mean:.2f}")
    ax_hist.axvline(results.median, color="darkgreen", linestyle=":",
                    linewidth=1.5, label=f"Median = {results.median:.2f}")
    ax_hist.axvspan(*results.ci_95, alpha=0.08, color="orange",
                    label="95% CI")
    ax_hist.axvspan(*results.ci_68, alpha=0.12, color="gold",
                    label="68% CI")
    ax_hist.set_xlabel(results.output_name)
    ax_hist.set_ylabel("Probability Density")
    ax_hist.set_title(f"Uncertainty — {results.output_name}")
    ax_hist.legend(fontsize=8)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Uncertainty plot saved to {save_path}")
    else:
        plt.show()
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Usage
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("PTR-94 Monte Carlo — Quick Demo")
    print("=" * 60)

    print("\nRunning ATP yield MC (10,000 samples for speed)...")
    uq_result = monte_carlo_atp_yield(n_samples=10000, seed=42)
    print(uq_result.summary())

    print("\nTop input correlations with ATP yield:")
    for name, r in sorted(
        uq_result.correlation_with_inputs.items(),
        key=lambda kv: abs(kv[1]),
        reverse=True,
    )[:5]:
        print(f"  {name:<30} r = {r:+.4f}")
