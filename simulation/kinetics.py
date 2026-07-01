#!/usr/bin/env python3
"""
PTR-94 Kinetics Module

Kinetic modelling of the PTR-94 pathway including Michaelis-Menten
enzyme kinetics, proton leak kinetics, ATP synthase rotary kinetics,
and pathway-level ODE integration for steady-state and time-course
simulations.

All concentrations in molar (M), times in seconds (s), rates in
M/s unless otherwise specified.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

try:
    from scipy.integrate import solve_ivp
    from scipy.optimize import fsolve
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# ---------------------------------------------------------------------------
# Kinetic parameters
# ---------------------------------------------------------------------------

@dataclass
class KineticParameter:
    """
    Enzyme kinetic parameters for a single reaction.

    Attributes
    ----------
    vmax : float
        Maximum forward reaction velocity (M/s).
    km : Dict[str, float]
        Michaelis constants for each substrate (M). Keyed by
        metabolite name.
    kcat : float
        Turnover number (1/s). Related to Vmax = kcat * [E]_total.
    enzyme_concentration : float
        Total enzyme concentration (M).
    ki : Dict[str, float]
        Inhibition constants (M). Keyed by inhibitor name.

    Notes
    -----
    ASSUMPTION: Michaelis-Menten formalism with optional competitive
    inhibition. Multi-substrate reactions use ordered sequential
    or ping-pong mechanisms as specified.
    """
    vmax: float = 1e-3
    km: Dict[str, float] = field(default_factory=dict)
    kcat: float = 1e2
    enzyme_concentration: float = 1e-6
    ki: Dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.vmax <= 0 and self.kcat > 0 and self.enzyme_concentration > 0:
            self.vmax = self.kcat * self.enzyme_concentration


# ---------------------------------------------------------------------------
# EnzymeKinetics
# ---------------------------------------------------------------------------

class EnzymeKinetics:
    """
    Collection of enzyme kinetic rate laws.

    All methods are static and return reaction rates in M/s.
    """

    @staticmethod
    def michaelis_menten(
        substrate_conc: float,
        vmax: float,
        km: float,
    ) -> float:
        """
        Classic Michaelis-Menten rate law.

        Parameters
        ----------
        substrate_conc : float
            Substrate concentration (M).
        vmax : float
            Maximum velocity (M/s).
        km : float
            Michaelis constant (M).

        Returns
        -------
        float
            Reaction rate in M/s.

        Notes
        -----
        v = Vmax * [S] / (Km + [S])

        ASSUMPTION: Irreversible reaction, single substrate,
        steady-state approximation. No product inhibition.
        """
        if substrate_conc < 0:
            return 0.0
        return vmax * substrate_conc / (km + substrate_conc)

    @staticmethod
    def competitive_inhibition(
        substrate_conc: float,
        inhibitor_conc: float,
        vmax: float,
        km: float,
        ki: float,
    ) -> float:
        """
        Michaelis-Menten with competitive inhibition.

        Parameters
        ----------
        substrate_conc : float
            Substrate concentration (M).
        inhibitor_conc : float
            Inhibitor concentration (M).
        vmax : float
            Maximum velocity (M/s).
        km : float
            Michaelis constant for substrate (M).
        ki : float
            Inhibition constant (M).

        Returns
        -------
        float
            Reaction rate in M/s.

        Notes
        -----
        v = Vmax * [S] / (Km * (1 + [I]/Ki) + [S])

        ASSUMPTION: Pure competitive inhibition. Inhibitor binds
        reversibly to the active site. No allosteric effects.
        """
        if substrate_conc < 0:
            return 0.0
        apparent_km = km * (1.0 + max(inhibitor_conc, 0.0) / max(ki, 1e-30))
        return vmax * substrate_conc / (apparent_km + substrate_conc)

    @staticmethod
    def allosteric_regulation(
        substrate_conc: float,
        vmax: float,
        k0_5: float,
        hill_coefficient: float,
    ) -> float:
        """
        Allosteric enzyme kinetics using the Hill equation.

        Parameters
        ----------
        substrate_conc : float
            Substrate concentration (M).
        vmax : float
            Maximum velocity (M/s).
        k0_5 : float
            Half-saturation constant (M). The [S] at v = Vmax/2.
        hill_coefficient : float
            Hill coefficient (n). n > 1 = positive cooperativity,
            n < 1 = negative cooperativity, n = 1 = hyperbolic.

        Returns
        -------
        float
            Reaction rate in M/s.

        Notes
        -----
        v = Vmax * [S]^n / (K0.5^n + [S]^n)

        ASSUMPTION: Homotropic allostery with a single type of
        binding site. The Hill coefficient is an empirical measure
        of cooperativity.
        """
        if substrate_conc <= 0:
            return 0.0
        s_n = substrate_conc ** hill_coefficient
        k_n = k0_5 ** hill_coefficient
        return vmax * s_n / (k_n + s_n)

    @staticmethod
    def reversible_michaelis_menten(
        substrate_conc: float,
        product_conc: float,
        vmax_forward: float,
        vmax_reverse: float,
        km_substrate: float,
        km_product: float,
    ) -> float:
        """
        Reversible Michaelis-Menten kinetics (Haldane relationship).

        Parameters
        ----------
        substrate_conc : float
            Substrate concentration (M).
        product_conc : float
            Product concentration (M).
        vmax_forward : float
            Maximum forward velocity (M/s).
        vmax_reverse : float
            Maximum reverse velocity (M/s).
        km_substrate : float
            Michaelis constant for substrate (M).
        km_product : float
            Michaelis constant for product (M).

        Returns
        -------
        float
            Net reaction rate in M/s. Positive = forward.

        Notes
        -----
        v = (Vmax_f * [S]/Km_s - Vmax_r * [P]/Km_p)
            / (1 + [S]/Km_s + [P]/Km_p)

        At equilibrium, Keq = (Vmax_f * Km_p) / (Vmax_r * Km_s).

        ASSUMPTION: Reversible Michaelis-Menten with competitive
        binding of substrate and product. No other intermediates.
        """
        if substrate_conc < 0 or product_conc < 0:
            return 0.0
        km_s = max(km_substrate, 1e-30)
        km_p = max(km_product, 1e-30)
        numerator = (
            vmax_forward * substrate_conc / km_s
            - vmax_reverse * product_conc / km_p
        )
        denominator = 1.0 + substrate_conc / km_s + product_conc / km_p
        return numerator / denominator

    @staticmethod
    def multi_substrate_ordered(
        substrate_a_conc: float,
        substrate_b_conc: float,
        vmax: float,
        km_a: float,
        km_b: float,
        ki_a: float,
    ) -> float:
        """
        Ordered sequential bi-substrate kinetics.

        Parameters
        ----------
        substrate_a_conc : float
            First (leading) substrate concentration (M).
        substrate_b_conc : float
            Second (following) substrate concentration (M).
        vmax : float
            Maximum velocity (M/s).
        km_a : float
            Michaelis constant for substrate A (M).
        km_b : float
            Michaelis constant for substrate B (M).
        ki_a : float
            Inhibition constant for substrate A (dissociation
            constant from the enzyme, M).

        Returns
        -------
        float
            Reaction rate in M/s.

        Notes
        -----
        v = Vmax * [A] * [B]
            / (Ki_A * Km_B + Km_A * [B] + [A] * [B])

        Ordered binding: A binds first, then B binds. Product
        release is also ordered.

        ASSUMPTION: Ordered sequential mechanism (e.g., many
        NAD-dependent dehydrogenases). No abortive complexes.
        """
        if substrate_a_conc <= 0 or substrate_b_conc <= 0:
            return 0.0
        denom = (
            ki_a * km_b
            + km_a * substrate_b_conc
            + substrate_a_conc * substrate_b_conc
        )
        if denom <= 0:
            return 0.0
        return vmax * substrate_a_conc * substrate_b_conc / denom

    @staticmethod
    def ping_pong_bi_bi(
        substrate_a_conc: float,
        substrate_b_conc: float,
        vmax: float,
        km_a: float,
        km_b: float,
    ) -> float:
        """
        Ping-pong (double displacement) bi-substrate kinetics.

        Parameters
        ----------
        substrate_a_conc : float
            First substrate concentration (M).
        substrate_b_conc : float
            Second substrate concentration (M).
        vmax : float
            Maximum velocity (M/s).
        km_a : float
            Michaelis constant for substrate A (M).
        km_b : float
            Michaelis constant for substrate B (M).

        Returns
        -------
        float
            Reaction rate in M/s.

        Notes
        -----
        v = Vmax * [A] * [B] / (Km_B * [A] + Km_A * [B] + [A] * [B])

        In a ping-pong mechanism, A binds and a product is released
        (modified enzyme intermediate), then B binds and the second
        product is released.

        ASSUMPTION: Classical ping-pong mechanism (e.g.,
        aminotransferases). Product inhibition not modelled.
        """
        if substrate_a_conc <= 0 or substrate_b_conc <= 0:
            return 0.0
        denom = (
            km_b * substrate_a_conc
            + km_a * substrate_b_conc
            + substrate_a_conc * substrate_b_conc
        )
        if denom <= 0:
            return 0.0
        return vmax * substrate_a_conc * substrate_b_conc / denom


# ---------------------------------------------------------------------------
# ProtonLeakKinetics
# ---------------------------------------------------------------------------

class ProtonLeakKinetics:
    """
    Models proton leak across the membrane as a function of the
    proton motive force.

    The leak represents the primary energy dissipation mechanism in
    chemiosmotic systems. In the PCM, leaks are minimised by
    engineered membrane composition.
    """

    def __init__(
        self,
        l0: float = 1e-8,       # basal conductance (M/s per mV)
        alpha: float = 0.08,     # exponential factor (1/mV)
        q10: float = 1.5,        # temperature coefficient
        temp_ref: float = 310.15,  # reference temperature (K)
    ) -> None:
        """
        Parameters
        ----------
        l0 : float
            Basal proton conductance (M/s per unit driving force).
        alpha : float
            Exponential factor characterising the voltage dependence
            of the leak (1/mV).
        q10 : float
            Factor by which leak rate increases per 10 K rise.
        temp_ref : float
            Reference temperature for Q10 scaling.

        ASSUMPTION: Leak is primarily driven by the proton motive
        force. The exponential factor alpha captures the non-ohmic
        behaviour observed in mitochondrial membranes (Nicholls &
        Ferguson, 2013). Basal conductance for synthetic membranes
        is assumed lower than natural membranes.
        """
        self.l0 = l0
        self.alpha = alpha
        self.q10 = q10
        self.temp_ref = temp_ref

    def leak_rate(
        self,
        pmf_kj_per_mol: float,
        temperature: float = 310.15,
    ) -> float:
        """
        Compute the proton leak rate.

        Parameters
        ----------
        pmf_kj_per_mol : float
            Proton motive force in kJ/mol. Larger PMF drives larger
            leak.
        temperature : float
            Current temperature in Kelvin.

        Returns
        -------
        float
            Proton leak rate in M/s (equivalent H+ flux).

        Notes
        -----
        J_leak = L0 * exp(alpha * PMF) * Q10_factor

        The non-ohmic (exponential) behaviour is well-documented in
        mitochondrial membranes. The Q10 factor accounts for
        temperature-dependent membrane fluidity.

        ASSUMPTION: Leak is a passive, non-saturable process. At
        very high PMF (>30 kJ/mol), leak becomes the dominant
        dissipative mechanism.
        """
        pmf_positive = max(pmf_kj_per_mol, 0.0)
        tc = (temperature - self.temp_ref) / 10.0
        q10_factor = self.q10 ** tc
        return self.l0 * math.exp(self.alpha * pmf_positive) * q10_factor

    def temperature_dependence(
        self,
        temp_range: np.ndarray,
        pmf_kj_per_mol: float = 20.0,
    ) -> np.ndarray:
        """
        Compute leak rate across a temperature range.

        Parameters
        ----------
        temp_range : np.ndarray
            Array of temperatures (K).
        pmf_kj_per_mol : float
            Fixed proton motive force (kJ/mol).

        Returns
        -------
        np.ndarray
            Leak rates at each temperature (M/s).
        """
        return np.array([
            self.leak_rate(pmf_kj_per_mol, float(t))
            for t in np.atleast_1d(temp_range)
        ])

    def effective_conductance(self, pmf_kj_per_mol: float) -> float:
        """
        Compute the effective proton conductance at a given PMF.

        Parameters
        ----------
        pmf_kj_per_mol : float
            Proton motive force (kJ/mol).

        Returns
        -------
        float
            Apparent conductance (M/s per kJ/mol).

        Notes
        -----
        C_H = J_leak / PMF

        ASSUMPTION: Linear approximation around a given operating
        point. Useful for small-signal analysis.
        """
        j = self.leak_rate(pmf_kj_per_mol)
        if abs(pmf_kj_per_mol) < 1e-10:
            return self.l0
        return j / pmf_kj_per_mol


# ---------------------------------------------------------------------------
# ATPSynthaseKinetics
# ---------------------------------------------------------------------------

class ATPSynthaseKinetics:
    """
    Kinetic model of ATP synthase driven by proton motive force.

    Models the rotary mechanism with H+/ATP stoichiometry dependence,
    slip probability, and load dependence (ATP/ADP ratio feedback).
    """

    def __init__(
        self,
        h_per_atp: float = 3.0,
        k_rotation: float = 100.0,    # 1/s, max rotation rate
        km_atp_site: float = 1e-5,    # M, affinity for ADP at catalytic site
        km_pi_site: float = 1e-3,     # M, affinity for Pi
        slip_probability: float = 0.001,  # fraction of rotations that slip
        atp_adp_ratio_threshold: float = 50.0,  # ratio above which slip increases
    ) -> None:
        """
        Parameters
        ----------
        h_per_atp : float
            Number of protons required per synthesised ATP.
        k_rotation : float
            Maximum rotation rate of the central stalk (1/s).
        km_atp_site : float
            Michaelis constant for ADP at the catalytic site (M).
        km_pi_site : float
            Michaelis constant for Pi at the catalytic site (M).
        slip_probability : float
            Fraction of 120-degree rotations that fail to produce
            ATP (dimensionless, 0-1).
        atp_adp_ratio_threshold : float
            [ATP]/[ADP] ratio above which slip probability
            increases due to product inhibition.

        ASSUMPTIONS:
        - The rotary mechanism proceeds in 120-degree steps, each
          consuming H+.
        - Slip is a stochastic decoupling event where H+ transit
          does not produce ATP.
        - The Km values account for physiological pH and [Mg2+].
        """
        self.h_per_atp = h_per_atp
        self.k_rotation = k_rotation
        self.km_atp_site = km_atp_site
        self.km_pi_site = km_pi_site
        self.slip_probability = slip_probability
        self.atp_adp_ratio_threshold = atp_adp_ratio_threshold

    def rotation_rate(
        self,
        pmf_kj_per_mol: float,
        adp_conc: float,
        pi_conc: float,
        atp_adp_ratio: float,
    ) -> float:
        """
        Compute the ATP synthesis rate based on proton-driven rotation.

        Parameters
        ----------
        pmf_kj_per_mol : float
            Proton motive force (kJ/mol).
        adp_conc : float
            ADP concentration (M).
        pi_conc : float
            Phosphate concentration (M).
        atp_adp_ratio : float
            [ATP] / [ADP] ratio (product inhibition).

        Returns
        -------
        float
            ATP synthesis rate in M/s (per enzyme), or negative
            if hydrolysis dominates.

        Notes
        -----
        The driving force for rotation:

            Delta_G_drive = PMF * h_per_atp - Delta_G_ATP_hydrolysis

        where Delta_G_ATP_hydrolysis is the actual free energy of
        ATP hydrolysis under current conditions.

        Rotation rate is modulated by substrate saturation:

            f_substrate = [ADP] / (Km_ADP + [ADP])
                         * [Pi] / (Km_Pi + [Pi])

        and by slip:

            effective_rate = k_rotation * f_substrate * (1 - slip)

        Slip increases at high ATP/ADP ratios:

            slip = slip_0 + (1 - slip_0) * (ratio - threshold)
                   / (ratio - threshold + K_slip)

        ASSUMPTION: The PMF-rotation relationship is approximately
        linear above a threshold of ~3 kJ/mol per H+. Below this,
        no rotation occurs.
        """
        if pmf_kj_per_mol <= 1e-6:
            return 0.0

        # Substrate saturation
        f_adp = adp_conc / (self.km_atp_site + adp_conc)
        f_pi = pi_conc / (self.km_pi_site + pi_conc)
        substrate_factor = f_adp * f_pi

        # Load-dependent slip
        ratio = max(atp_adp_ratio, 0.0)
        if ratio > self.atp_adp_ratio_threshold:
            excess = ratio - self.atp_adp_ratio_threshold
            k_slip = 100.0  # transition steepness
            slip = min(
                self.slip_probability
                + (1.0 - self.slip_probability)
                * excess / (excess + k_slip),
                1.0,
            )
        else:
            slip = self.slip_probability

        # Effective rotation rate
        rate = self.k_rotation * substrate_factor * (1.0 - slip)
        return max(rate, 0.0)

    def proton_flux(
        self,
        rotation_rate: float,
    ) -> float:
        """
        Compute the proton flux through ATP synthase.

        Parameters
        ----------
        rotation_rate : float
            ATP synthesis rate (M/s).

        Returns
        -------
        float
            Proton consumption rate (M/s). Positive = H+ flowing
            through the synthase from intermembrane space to matrix.
        """
        return rotation_rate * self.h_per_atp

    def power_output(
        self,
        rate: float,
        delta_g_atp_synthesis: float = 30.5,
    ) -> float:
        """
        Chemical power output of ATP synthase.

        Parameters
        ----------
        rate : float
            ATP synthesis rate (M/s).
        delta_g_atp_synthesis : float
            Gibbs free energy stored per ATP (kJ/mol).

        Returns
        -------
        float
            Power in kW/m^3 (kJ/(L*s) = MW/m^3 = 1e3 kW/m^3).
        """
        # Convert: M/s = mol/(L*s), delta_G in kJ/mol
        # Power = rate * delta_G (kJ/(L*s) = kW/L = 1000 kW/m^3)
        return rate * delta_g_atp_synthesis * 1000.0


# ---------------------------------------------------------------------------
# PathwayKinetics
# ---------------------------------------------------------------------------

class PathwayKinetics:
    """
    Kinetic model of the complete PTR-94 pathway.

    Supports steady-state flux calculation and ODE-based time-course
    simulation of metabolite concentrations.
    """

    def __init__(self) -> None:
        # Named enzyme kinetics functions
        self.enzymes: Dict[str, Callable[..., float]] = {}
        # Metabolite concentrations (M) indexed by name
        self.metabolites: Dict[str, float] = {}
        # Flux through each reaction (M/s)
        self.fluxes: Dict[str, float] = {}
        # Parameters for each enzyme
        self.params: Dict[str, KineticParameter] = {}

        # Proton leak and ATP synthase models
        self.proton_leak = ProtonLeakKinetics()
        self.atp_synthase = ATPSynthaseKinetics()

    def set_metabolite(self, name: str, concentration: float) -> None:
        """Set the concentration of a metabolite."""
        self.metabolites[name] = max(concentration, 0.0)

    def get_metabolite(self, name: str, default: float = 0.0) -> float:
        """Get the concentration of a metabolite."""
        return self.metabolites.get(name, default)

    def add_enzyme(
        self,
        name: str,
        kinetics_func: Callable[..., float],
        params: KineticParameter,
    ) -> None:
        """
        Register an enzyme with its kinetic parameters.

        Parameters
        ----------
        name : str
            Enzyme/reaction name.
        kinetics_func : callable
            Rate function (from EnzymeKinetics or custom).
        params : KineticParameter
            Kinetic parameters.
        """
        self.enzymes[name] = kinetics_func
        self.params[name] = params

    def compute_flux_michaelis_menten(
        self,
        enzyme_name: str,
        substrate_name: str,
    ) -> float:
        """
        Compute flux for a simple Michaelis-Menten reaction.

        Parameters
        ----------
        enzyme_name : str
            Registered enzyme name.
        substrate_name : str
            Name of the substrate metabolite.

        Returns
        -------
        float
            Flux in M/s.
        """
        p = self.params[enzyme_name]
        s = self.get_metabolite(substrate_name)
        return EnzymeKinetics.michaelis_menten(
            s, p.vmax, p.km.get(substrate_name, 1e-3)
        )

    def compute_steady_state_flux(
        self,
        input_metabolite: str = "Glucose",
        input_rate: float = 1e-6,  # M/s
        tolerance: float = 1e-12,
        max_iterations: int = 100,
    ) -> Dict[str, float]:
        """
        Compute steady-state fluxes through the pathway.

        Uses an iterative approach: given an input rate, find the
        flux distribution that balances all metabolite pools.

        Parameters
        ----------
        input_metabolite : str
            Entry point metabolite.
        input_rate : float
            Supply rate of input metabolite (M/s).
        tolerance : float
            Convergence tolerance for metabolite residuals.
        max_iterations : int
            Maximum iteration count.

        Returns
        -------
        dict
            Reaction name -> steady-state flux (M/s).

        Notes
        -----
        ASSUMPTION: Simple substrate-limited kinetics. No regulation,
        no cofactor limitations. Useful for approximate yield
        calculations.

        The steady-state condition is: S * v = 0 (mass balance).
        We approximate by matching consumption and production rates.
        """
        if not HAS_SCIPY:
            return self._steady_state_iterative(
                input_metabolite, input_rate, tolerance, max_iterations
            )
        return self._steady_state_via_fsolve(
            input_metabolite, input_rate
        )

    def _steady_state_iterative(
        self,
        input_metabolite: str,
        input_rate: float,
        tolerance: float,
        max_iterations: int,
    ) -> Dict[str, float]:
        """Iterative steady-state solver (no scipy fallback)."""
        self.fluxes = {
            name: input_rate for name in self.enzymes
        }
        for iteration in range(max_iterations):
            max_delta = 0.0
            for name in self.enzymes:
                old_flux = self.fluxes.get(name, 0.0)
                if name in self.params:
                    p = self.params[name]
                    # Try to find consumed metabolites
                    consumed = [
                        m for m in p.km if self.get_metabolite(m, 0) > 0
                    ]
                    if consumed:
                        s = self.get_metabolite(consumed[0])
                        new_flux = EnzymeKinetics.michaelis_menten(
                            s, p.vmax, p.km.get(consumed[0], 1e-3)
                        )
                        self.fluxes[name] = new_flux
                        delta = abs(new_flux - old_flux)
                        if delta > max_delta:
                            max_delta = delta
            if max_delta < tolerance:
                break
        return self.fluxes

    def _steady_state_via_fsolve(
        self,
        input_metabolite: str,
        input_rate: float,
    ) -> Dict[str, float]:
        """Steady-state solver using scipy.fsolve."""
        enzyme_names = list(self.enzymes.keys())
        n_enzymes = len(enzyme_names)

        def equations(vars_flat: np.ndarray) -> np.ndarray:
            fluxes = dict(zip(enzyme_names, vars_flat))
            residuals = np.zeros(n_enzymes)
            for i, name in enumerate(enzyme_names):
                p = self.params.get(name)
                if p is None:
                    residuals[i] = fluxes[name]
                    continue
                # Find primary substrate
                substrates = [m for m in p.km if self.get_metabolite(m, 0) > 0]
                if not substrates:
                    residuals[i] = fluxes[name]
                    continue
                s = self.get_metabolite(substrates[0])
                expected = EnzymeKinetics.michaelis_menten(
                    s, p.vmax, p.km[substrates[0]]
                )
                residuals[i] = fluxes[name] - expected
            return residuals

        initial_guess = np.full(n_enzymes, input_rate)
        try:
            solution = fsolve(equations, initial_guess, xtol=1e-8, maxfev=1000)
            self.fluxes = dict(zip(enzyme_names, solution))
        except Exception:
            self.fluxes = dict(zip(enzyme_names, initial_guess))
        return self.fluxes

    def time_course(
        self,
        initial_conditions: Dict[str, float],
        t_span: Tuple[float, float] = (0.0, 100.0),
        t_eval: Optional[np.ndarray] = None,
        method: str = "RK45",
    ) -> Dict[str, np.ndarray]:
        """
        Simulate metabolite concentration time course using ODE
        integration.

        Parameters
        ----------
        initial_conditions : dict
            Metabolite name -> initial concentration (M).
        t_span : tuple
            (t_start, t_end) in seconds.
        t_eval : np.ndarray, optional
            Time points for output. Default: 100 evenly spaced points.
        method : str
            ODE solver method for scipy.integrate.solve_ivp.
            Default: 'RK45'.

        Returns
        -------
        dict
            Keys: 't' (time array), and metabolite names (arrays).
            Each value is a 1D numpy array.

        Notes
        -----
        Integrates the ODE system:

            d[M_i]/dt = sum_j S_ij * v_j(t, [M])

        where S_ij is the stoichiometric matrix and v_j are
        reaction rates computed from enzyme kinetics.

        Uses scipy.integrate.solve_ivp with dense output.

        ASSUMPTION: Well-mixed compartments. No spatial gradients.
        Enzyme concentrations are constant (no gene expression
        dynamics). Temperature is constant at 310.15 K.
        """
        if not HAS_SCIPY:
            return self._time_course_euler(initial_conditions, t_span)

        # Set initial conditions
        for name, conc in initial_conditions.items():
            self.set_metabolite(name, conc)

        met_names = sorted(initial_conditions.keys())
        n_mets = len(met_names)

        if t_eval is None:
            t_eval = np.linspace(t_span[0], t_span[1], 100)

        def ode_system(t: float, y: np.ndarray) -> np.ndarray:
            # Update metabolite pool
            for j, name in enumerate(met_names):
                self.set_metabolite(name, max(float(y[j]), 0.0))

            # Compute all fluxes
            dydt = np.zeros(n_mets)
            for enzyme_name in self.enzymes:
                try:
                    flux = self.compute_flux_michaelis_menten(
                        enzyme_name,
                        list(self.params[enzyme_name].km.keys())[0]
                        if enzyme_name in self.params
                        and self.params[enzyme_name].km else ""
                    )
                except Exception:
                    flux = 0.0
                # Apply stoichiometry (crude: assume each enzyme
                # consumes one substrate and produces one product)
                p = self.params.get(enzyme_name)
                if p is None:
                    continue
                for substrate in p.km:
                    if substrate in met_names:
                        si = met_names.index(substrate)
                        dydt[si] -= flux
                        # Find product (any metabolite not a substrate)
                        for prod in met_names:
                            if prod not in p.km:
                                pi = met_names.index(prod)
                                dydt[pi] += flux
                                break
                        break
            return dydt

        try:
            result = solve_ivp(
                ode_system,
                t_span,
                np.array([initial_conditions[m] for m in met_names]),
                method=method,
                t_eval=t_eval,
                rtol=1e-6,
                atol=1e-12,
            )
            output: Dict[str, np.ndarray] = {"t": result.t}
            for i, name in enumerate(met_names):
                output[name] = result.y[i]
            return output
        except Exception:
            return self._time_course_euler(initial_conditions, t_span)

    def _time_course_euler(
        self,
        initial_conditions: Dict[str, float],
        t_span: Tuple[float, float],
        n_steps: int = 1000,
    ) -> Dict[str, np.ndarray]:
        """Simple forward Euler integration as fallback."""
        met_names = sorted(initial_conditions.keys())
        t = np.linspace(t_span[0], t_span[1], n_steps)
        dt = t[1] - t[0]
        y = np.zeros((len(met_names), n_steps))
        for i, name in enumerate(met_names):
            y[i, 0] = initial_conditions[name]

        for step in range(1, n_steps):
            for i, name in enumerate(met_names):
                self.set_metabolite(name, max(float(y[i, step - 1]), 0.0))
            for i, name in enumerate(met_names):
                flux = 0.0
                for ename in self.enzymes:
                    p = self.params.get(ename)
                    if p is None:
                        continue
                    for subs in p.km:
                        if subs == name:
                            try:
                                f = EnzymeKinetics.michaelis_menten(
                                    self.get_metabolite(subs, 0),
                                    p.vmax, p.km[subs]
                                )
                                flux -= f
                            except Exception:
                                pass
                            break
                        elif name not in p.km:
                            # Potential product
                            pass
                y[i, step] = max(y[i, step - 1] + dt * flux, 0.0)

        output: Dict[str, np.ndarray] = {"t": t}
        for i, name in enumerate(met_names):
            output[name] = y[i]
        return output

    def atp_production_rate(self) -> float:
        """
        Compute the total ATP production rate.

        Scans all enzyme fluxes for reactions that produce ATP.

        Returns
        -------
        float
            Net ATP production rate in M/s.
        """
        rate = 0.0
        # Use ATP synthase model for a more accurate estimate
        adp = self.get_metabolite("ADP", 1e-3)
        pi = self.get_metabolite("Pi", 5e-3)
        atp = self.get_metabolite("ATP", 3e-3)
        atp_adp_ratio = atp / max(adp, 1e-30)

        # Use typical PMF estimate
        pmf = self.get_metabolite("PMF", 20.0)

        rate = self.atp_synthase.rotation_rate(
            pmf_kj_per_mol=pmf,
            adp_conc=adp,
            pi_conc=pi,
            atp_adp_ratio=atp_adp_ratio,
        )
        return rate

    def nadh_turnover_rate(self) -> float:
        """
        Compute the NADH consumption (oxidation) rate.

        Returns
        -------
        float
            NADH consumption rate in M/s.
        """
        nadh = self.get_metabolite("NADH", 0.0)
        # ASSUMPTION: NADH is consumed by Complex I in the PCM
        p = self.params.get("PCM_CI")
        if p is not None:
            return EnzymeKinetics.michaelis_menten(
                nadh, p.vmax, p.km.get("NADH", 1e-4)
            )
        return 0.0

    def yield_atp_per_glucose(self) -> float:
        """
        Estimate ATP yield per glucose from the kinetic model.

        Returns
        -------
        float
            Estimated ATP molecules per glucose.

        Notes
        -----
        ASSUMPTION: The ratio of ATP production rate to glucose
        consumption rate gives a steady-state yield estimate.
        This is a simplified approach; real yield requires
        full stoichiometric accounting.
        """
        glucose_flux = self.get_metabolite("Glucose_flux", 1e-6)
        atp_flux = self.atp_production_rate()
        if glucose_flux > 0:
            return atp_flux / glucose_flux
        return 0.0
