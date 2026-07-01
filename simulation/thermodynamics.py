#!/usr/bin/env python3
"""
PTR-94 Thermodynamics Module

Rigorous thermodynamic calculations for the Perfect Coupling Module and
theoretical maximum ATP yield from glucose oxidation.

All energies in kJ/mol unless otherwise specified. Standard biochemical
conditions: 298.15 K, 1 atm, pH 7.0, [Mg2+] = 1 mM, ionic strength 0.25 M.
Physiological temperature taken as 310.15 K (37 degrees C) where applicable.

References:
  - Alberty, R.A. (2003) Thermodynamics of Biochemical Reactions
  - Nicholls, D.G. & Ferguson, S.J. (2013) Bioenergetics, 4th Ed.
  - Berg, J.M. et al. (2015) Biochemistry, 8th Ed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Tuple
import numpy as np


# ---------------------------------------------------------------------------
# Fundamental constants
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ThermodynamicConstants:
    """Grouped thermodynamic and biophysical constants."""
    R: float = 8.314462618       # J/(mol*K) — Gas constant
    F: float = 96485.335         # C/mol — Faraday constant
    T_std: float = 298.15        # K — Standard temperature (25 C)
    T_phys: float = 310.15       # K — Physiological temperature (37 C)
    pH_std: float = 7.0          # — Standard biochemical pH


CONST = ThermodynamicConstants()

# Standard Gibbs free energy of glucose oxidation (complete to CO2 + H2O)
#   ASSUMPTION: Standard biochemical conditions (pH 7, 298.15 K, 1 M solutes,
#   excluding water). The value -2870 kJ/mol is the consensus from Berg et al.
#   and Alberty. Some sources quote -2840 to -2880 kJ/mol.
STANDARD_GIBBS_GLUCOSE: float = 2870.0  # kJ/mol (absolute value, exergonic)

# Standard Gibbs free energy of ATP synthesis from ADP + Pi
#   ASSUMPTION: Standard biochemical conditions. The widely accepted value
#   is +30.5 kJ/mol. Under cellular conditions this ranges from +45 to +65.
STANDARD_GIBBS_ATP: float = 30.5  # kJ/mol

# Typical physiological Gibbs free energy of ATP hydrolysis
#   ASSUMPTION: In vivo [ATP]/[ADP] ~ 10, [Pi] ~ 1-5 mM, giving
#   Delta_G ~ -55 to -60 kJ/mol. We use 55 kJ/mol as typical.
PHYSIOLOGICAL_GIBBS_ATP: float = 55.0  # kJ/mol (absolute value, exergonic)

FARADAY_CONSTANT: float = CONST.F  # C/mol — exposed for convenience
GAS_CONSTANT: float = CONST.R      # J/(mol*K) — exposed for convenience


# ---------------------------------------------------------------------------
# Core thermodynamic computations
# ---------------------------------------------------------------------------

def compute_theoretical_max_atp(
    delta_g_glucose: float = STANDARD_GIBBS_GLUCOSE,
    delta_g_atp: float = STANDARD_GIBBS_ATP,
) -> float:
    """
    Compute the theoretical maximum number of ATP molecules that can be
    synthesised from the complete oxidation of one glucose molecule.

    Parameters
    ----------
    delta_g_glucose : float
        Gibbs free energy released by glucose oxidation (kJ/mol).
        Default: STANDARD_GIBBS_GLUCOSE = 2870.0
    delta_g_atp : float
        Gibbs free energy required to synthesise one ATP (kJ/mol).
        Default: STANDARD_GIBBS_ATP = 30.5

    Returns
    -------
    float
        Maximum ATP yield (dimensionless). A fractional value indicates
        the thermodynamic ceiling; the integer floor is the practical max.

    Notes
    -----
    The theoretical ceiling is simply the ratio of available energy to
    the cost per ATP: n_max = Delta_G_glucose / Delta_G_ATP.

    Under standard conditions: 2870 / 30.5 = 94.1
    Under physiological conditions: 2870 / 55.0 = 52.2
    """
    return delta_g_glucose / delta_g_atp


def compute_coupling_efficiency(
    atp_yield: float,
    delta_g_glucose: float = STANDARD_GIBBS_GLUCOSE,
    delta_g_atp: float = STANDARD_GIBBS_ATP,
) -> float:
    """
    Compute the thermodynamic coupling efficiency of ATP synthesis.

    Parameters
    ----------
    atp_yield : float
        Number of ATP molecules produced per glucose.
    delta_g_glucose : float
        Gibbs free energy released by glucose oxidation (kJ/mol).
    delta_g_atp : float
        Gibbs free energy required per ATP (kJ/mol).

    Returns
    -------
    float
        Efficiency as a percentage (0-100).

    Notes
    -----
    Efficiency = (ATP_yield * Delta_G_ATP / Delta_G_glucose) * 100
    """
    energy_captured = atp_yield * delta_g_atp
    return (energy_captured / delta_g_glucose) * 100.0


def nernst_potential(
    delta_psi: float,
    ion_charge: int,
    temp: float = CONST.T_phys,
) -> float:
    """
    Compute the ion concentration ratio at equilibrium for a given
    membrane potential using the Nernst equation.

    Parameters
    ----------
    delta_psi : float
        Membrane potential in volts (V). Positive means inside negative
        relative to outside (typical for mitochondria).
    ion_charge : int
        Charge number of the ion (e.g., +1 for H+, -1 for Cl-).
    temp : float
        Temperature in Kelvin (default: 310.15 K).

    Returns
    -------
    float
        Equilibrium concentration ratio [C_out] / [C_in].

    Notes
    -----
    Derived from: Delta_psi = (R*T / z*F) * ln([C_out] / [C_in]).

    Re-arranged: [C_out] / [C_in] = exp(Delta_psi * z * F / (R * T))

    ASSUMPTION: The ion is at equilibrium across the membrane; no
    active transport contributions.
    """
    rt = CONST.R * temp
    exponent = delta_psi * ion_charge * CONST.F / rt
    return math.exp(exponent)


def proton_motive_force(
    delta_psi: float,
    delta_ph: float,
    temp: float = CONST.T_phys,
) -> float:
    """
    Calculate the proton motive force across a membrane.

    Parameters
    ----------
    delta_psi : float
        Membrane potential in volts. Positive = inside negative.
    delta_ph : float
        Transmembrane pH gradient (pH_out - pH_in). Positive means the
        outside is more acidic (typical for mitochondria).
    temp : float
        Temperature in Kelvin (default: 310.15 K).

    Returns
    -------
    float
        Proton motive force in kJ/mol. Positive PMF means protons are
        driven from the intermembrane space to the matrix.

    Notes
    -----
    PMF (in J/mol) = F * (Delta_psi - (2.303 * R * T / F) * Delta_pH)

    The factor 2.303 * R * T / F converts Delta_pH to voltage:
      - At 310.15 K: ~0.06154 V / pH unit
      - At 298.15 K: ~0.05916 V / pH unit

    ASSUMPTION: Ideal solution, linear response.
    """
    thermal_voltage = 2.303 * CONST.R * temp / CONST.F
    pmf_volts = delta_psi - thermal_voltage * delta_ph
    return CONST.F * pmf_volts / 1000.0


def gibbs_from_ph_and_potential(
    delta_psi: float,
    delta_ph: float,
    h_stoichiometry: int,
    temp: float = CONST.T_phys,
) -> float:
    """
    Compute the Gibbs free energy available from a proton gradient
    for a given H+ stoichiometry.

    Parameters
    ----------
    delta_psi : float
        Membrane potential in volts. Positive = inside negative.
    delta_ph : float
        Transmembrane pH gradient (pH_out - pH_in).
    h_stoichiometry : int
        Number of protons translocated per process (e.g., per NADH
        oxidised, or per ATP synthesised).
    temp : float
        Temperature in Kelvin (default: 310.15 K).

    Returns
    -------
    float
        Gibbs free energy in kJ/mol.

    Notes
    -----
    Delta_G = n * F * (Delta_psi - (2.303*R*T/F) * Delta_pH) / 1000

    ASSUMPTION: The gradient is the sole energy intermediate (delocalised
    chemiosmotic coupling). Direct coupling mechanisms may differ.
    """
    pmf = proton_motive_force(delta_psi, delta_ph, temp)
    return float(h_stoichiometry) * pmf


def entropy_production(
    delta_g_total: float,
    delta_g_captured: float,
    temp: float = CONST.T_phys,
) -> float:
    """
    Calculate the entropy production (irreversibility) of an energy
    transduction process.

    Parameters
    ----------
    delta_g_total : float
        Total Gibbs free energy released by the process (absolute value,
        kJ/mol). Must be > delta_g_captured.
    delta_g_captured : float
        Gibbs free energy captured in useful work/products (kJ/mol).
    temp : float
        Temperature in Kelvin (default: 310.15 K).

    Returns
    -------
    float
        Entropy production in J/(K*mol).

    Notes
    -----
    S_prod = (Delta_G_total - Delta_G_captured) / T

    ASSUMPTION: All non-captured energy is dissipated as heat.
    """
    if delta_g_captured > delta_g_total:
        raise ValueError(
            f"Captured energy ({delta_g_captured} kJ/mol) cannot exceed "
            f"total energy ({delta_g_total} kJ/mol)."
        )
    delta_g_dissipated = delta_g_total - delta_g_captured
    return delta_g_dissipated * 1000.0 / temp


def chemical_potential(
    concentration: float,
    standard_potential: float = 0.0,
    temp: float = CONST.T_phys,
) -> float:
    """
    Calculate the chemical potential of a solute at a given concentration.

    Parameters
    ----------
    concentration : float
        Concentration in molar (mol/L). Must be > 0.
    standard_potential : float
        Standard-state chemical potential (kJ/mol). Default 0.0.
    temp : float
        Temperature in Kelvin (default: 310.15 K).

    Returns
    -------
    float
        Chemical potential in kJ/mol.

    Notes
    -----
    mu = mu_0 + R * T * ln(c / c_0) where c_0 = 1 M (standard state).

    ASSUMPTION: Ideal dilute solution; activity coefficient = 1.
    """
    if concentration <= 0.0:
        raise ValueError(f"Concentration must be positive, got {concentration}")
    rt_term = CONST.R * temp / 1000.0
    return standard_potential + rt_term * math.log(concentration)


def redox_potential(
    e0_prime: float,
    n_electrons: int,
    oxidized: float,
    reduced: float,
    temp: float = CONST.T_phys,
) -> float:
    """
    Calculate the actual reduction potential under non-standard conditions
    using the Nernst equation.

    Parameters
    ----------
    e0_prime : float
        Standard reduction potential at pH 7 (V).
    n_electrons : int
        Number of electrons transferred.
    oxidized : float
        Concentration of the oxidised species (M). Must be > 0.
    reduced : float
        Concentration of the reduced species (M). Must be > 0.
    temp : float
        Temperature in Kelvin (default: 310.15 K).

    Returns
    -------
    float
        Actual reduction potential in volts.

    Notes
    -----
    E = E0' + (R * T / (n * F)) * ln([ox] / [red])

    ASSUMPTION: Activity coefficients are unity; pH effects are absorbed
    into E0' (which is already pH 7-corrected).
    """
    if oxidized <= 0.0:
        raise ValueError(f"Oxidised concentration must be > 0, got {oxidized}")
    if reduced <= 0.0:
        raise ValueError(f"Reduced concentration must be > 0, got {reduced}")
    rt_over_nf = CONST.R * temp / (float(n_electrons) * CONST.F)
    return e0_prime + rt_over_nf * math.log(oxidized / reduced)


def electrochemical_gradient(
    concentration_in: float,
    concentration_out: float,
    membrane_potential: float,
    ion_charge: int,
    temp: float = CONST.T_phys,
) -> float:
    """
    Compute the total electrochemical potential difference for an ion
    across a membrane.

    Parameters
    ----------
    concentration_in : float
        Intracellular (or matrix) concentration (M). Must be > 0.
    concentration_out : float
        Extracellular (or intermembrane space) concentration (M). Must be > 0.
    membrane_potential : float
        Membrane potential in volts (inside - outside convention;
        positive means inside negative).
    ion_charge : int
        Charge number of the ion.
    temp : float
        Temperature in Kelvin (default: 310.15 K).

    Returns
    -------
    float
        Electrochemical potential difference in kJ/mol (inside - outside).

    Notes
    -----
    Delta_mu = z * F * Delta_psi + R * T * ln(c_in / c_out)

    ASSUMPTION: Membrane is impermeable to the ion except through
    specific channels/transporters.
    """
    if concentration_in <= 0.0:
        raise ValueError(f"Inside concentration must be > 0, got {concentration_in}")
    if concentration_out <= 0.0:
        raise ValueError(f"Outside concentration must be > 0, got {concentration_out}")
    electrical = ion_charge * CONST.F * membrane_potential / 1000.0
    chemical = CONST.R * temp / 1000.0 * math.log(concentration_in / concentration_out)
    return electrical + chemical


def atp_hydrolysis_free_energy(
    atp_adp_ratio: float,
    pi_concentration: float,
    temp: float = CONST.T_phys,
) -> float:
    """
    Calculate the actual Gibbs free energy of ATP hydrolysis under
    cellular (non-standard) conditions.

    Parameters
    ----------
    atp_adp_ratio : float
        Cellular [ATP] / [ADP] ratio. Typical range: 5-100.
    pi_concentration : float
        Inorganic phosphate concentration in M. Typical: 0.001-0.01 M.
    temp : float
        Temperature in Kelvin (default: 310.15 K).

    Returns
    -------
    float
        Gibbs free energy change in kJ/mol. Negative means ATP hydrolysis
        is exergonic.

    Notes
    -----
    Delta_G = Delta_G0' + R * T * ln([Pi] / (ATP/ADP_ratio))

    Under typical conditions ([ATP]/[ADP] = 10, [Pi] = 1 mM):
        Delta_G = -30.5 + RT*ln(0.001/10) = -54.0 kJ/mol

    ASSUMPTION: Delta_G0' = -30.5 kJ/mol (standard biochemical).
    """
    if atp_adp_ratio <= 0.0:
        raise ValueError(f"ATP/ADP ratio must be > 0, got {atp_adp_ratio}")
    if pi_concentration <= 0.0:
        raise ValueError(f"Pi concentration must be > 0, got {pi_concentration}")
    rt = CONST.R * temp / 1000.0
    delta_g_std = -STANDARD_GIBBS_ATP
    mass_action_ratio = pi_concentration / atp_adp_ratio
    return delta_g_std + rt * math.log(mass_action_ratio)


def max_work_from_redox(
    electron_donor_potential: float,
    electron_acceptor_potential: float,
    n_electrons: int,
) -> float:
    """
    Calculate the maximum useful work obtainable from a redox reaction.

    Parameters
    ----------
    electron_donor_potential : float
        Reduction potential of the electron donor (E0', V).
    electron_acceptor_potential : float
        Reduction potential of the electron acceptor (E0', V).
    n_electrons : int
        Number of electrons transferred.

    Returns
    -------
    float
        Maximum work available in kJ/mol. Positive value indicates
        energy released.

    Notes
    -----
    Delta_G = -n * F * Delta_E, where Delta_E = E_acceptor - E_donor.

    For NADH -> O2: Delta_G = -2 * 96485 * 1.14 = -220 kJ/mol
    For FADH2 -> O2: Delta_G = -2 * 96485 * 0.84 = -162 kJ/mol

    ASSUMPTION: Standard reduction potentials at pH 7 (E0').
    """
    delta_e = electron_acceptor_potential - electron_donor_potential
    delta_g_joules = -float(n_electrons) * CONST.F * delta_e
    return delta_g_joules / 1000.0


def thermodynamic_ceiling(
    available_energy: float,
    atp_cost: float = STANDARD_GIBBS_ATP,
) -> int:
    """
    Compute the integer floor of the maximum ATP yield given an
    available free energy budget.

    Parameters
    ----------
    available_energy : float
        Free energy available for ATP synthesis (kJ/mol).
    atp_cost : float
        Free energy required per ATP synthesised (kJ/mol).
        Default: STANDARD_GIBBS_ATP = 30.5

    Returns
    -------
    int
        Maximum integer number of ATP molecules. This is the floor
        of the thermodynamic ratio.

    Notes
    -----
    n_max = floor(available_energy / atp_cost)

    For PTR-94: floor(2870 / 30.5) = floor(94.1) = 94

    ASSUMPTION: No energy is required for ATP export, transport, or
    other coupled processes.
    """
    if available_energy <= 0.0:
        return 0
    return int(available_energy // atp_cost)


# ---------------------------------------------------------------------------
# Convenience constants reference
# ---------------------------------------------------------------------------

STANDARD_REDOX_POTENTIALS: dict[str, Tuple[float, int]] = {
    "NAD_NADH": (-0.320, 2),
    "FAD_FADH2": (-0.220, 2),
    "ubiquinone_ubiquinol": (0.040, 2),
    "O2_H2O": (0.820, 4),
    "cytochrome_c_ox_red": (0.250, 1),
    "fumarate_succinate": (0.030, 2),
    "pyruvate_lactate": (-0.190, 2),
    "glutathione_GSSG_GSH": (-0.240, 2),
    "ferredoxin_ox_red": (-0.430, 1),
    "plastoquinone_ox_red": (0.110, 2),
}
"""
ASSUMPTION: These are textbook values from Berg et al. (2015) and
Nicholls & Ferguson (2013). Values vary slightly between sources.
"""

NATURAL_COMPARISON: dict[str, dict[str, float]] = {
    "eukaryote": {
        "atp_per_glucose": 30.5,
        "redox_atp": 28.0,
        "h_per_nadh": 10.0,
        "h_per_atp": 3.7,
        "efficiency_pct": 32.4,
    },
    "prokaryote": {
        "atp_per_glucose": 37.0,
        "redox_atp": 34.0,
        "h_per_nadh": 10.0,
        "h_per_atp": 3.3,
        "efficiency_pct": 39.4,
    },
    "ptr94": {
        "atp_per_glucose": 94.0,
        "redox_atp": 90.0,
        "h_per_nadh": 30.0,
        "h_per_atp": 3.0,
        "efficiency_pct": 99.9,
    },
}
"""
ASSUMPTION: Eukaryotic yield (30-32) accounts for transport costs;
prokaryotic yield (36-38) has no mitochondrial shuttles.
"""
