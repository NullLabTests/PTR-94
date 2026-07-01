#!/usr/bin/env python3
"""
PTR-94 Reaction Network Module

Defines and analyses metabolic reaction networks with full stoichiometric
tracking. Supports mass balance, charge balance, stoichiometric matrix
construction, nullspace analysis, and pathway enumeration.

Designed for both natural respiratory pathways and the engineered PCM.
"""

from __future__ import annotations

import math
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

try:
    from scipy.linalg import null_space
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# ---------------------------------------------------------------------------
# Core data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Reaction:
    """
    A biochemical reaction with thermodynamic and enzymatic annotation.

    Attributes
    ----------
    name : str
        Human-readable reaction name.
    reactants : Dict[str, float]
        Reactant metabolite names -> stoichiometric coefficient.
    products : Dict[str, float]
        Product metabolite names -> stoichiometric coefficient.
    delta_g_prime : float
        Standard transformed Gibbs free energy (kJ/mol).
        Negative = exergonic (favours forward direction).
    enzyme : str
        Enzyme catalysing the reaction (or 'spontaneous').
    reversibility : bool
        True if the reaction is considered reversible under cellular
        conditions.
    """
    name: str
    reactants: Dict[str, float]
    products: Dict[str, float]
    delta_g_prime: float = 0.0
    enzyme: str = "unknown"
    reversibility: bool = False

    @property
    def all_metabolites(self) -> Set[str]:
        """Set of all metabolite names in this reaction."""
        return set(self.reactants) | set(self.products)

    @property
    def net_stoichiometry(self) -> Dict[str, float]:
        """Net stoichiometric coefficients (products - reactants)."""
        net: Dict[str, float] = {}
        for m, c in self.reactants.items():
            net[m] = net.get(m, 0.0) - c
        for m, c in self.products.items():
            net[m] = net.get(m, 0.0) + c
        return net

    def __repr__(self) -> str:
        left = " + ".join(
            f"{'' if c == 1 else f'{c:.0f} '}{m}"
            for m, c in self.reactants.items()
        )
        arrow = " <-> " if self.reversibility else " -> "
        right = " + ".join(
            f"{'' if c == 1 else f'{c:.0f} '}{m}"
            for m, c in self.products.items()
        )
        return f"{left}{arrow}{right}  [{self.name}, {self.enzyme}]"


@dataclass(frozen=True)
class Metabolite:
    """
    A metabolite with formula, typical concentration, and compartment.

    Attributes
    ----------
    name : str
        Canonical name (matches Reaction reactant/product keys).
    formula : str
        Elemental formula (e.g., 'C6H12O6' for glucose).
    concentration_default : float
        Typical cellular concentration in molar (M).
    compartment : str
        Cellular compartment ('cytosol', 'matrix', 'ims', 'mitochondrial').
    charge : int
        Net charge at pH 7 (default 0).
    """
    name: str
    formula: str = ""
    concentration_default: float = 1e-3
    compartment: str = "cytosol"
    charge: int = 0


# ---------------------------------------------------------------------------
# Elemental composition helpers
# ---------------------------------------------------------------------------

ELEMENT_ATOMIC_MASSES: Dict[str, float] = {
    "C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007,
    "P": 30.974, "S": 32.065, "Mg": 24.305, "Fe": 55.845,
}


def parse_formula(formula: str) -> Dict[str, float]:
    """
    Parse a chemical formula string into element -> count mapping.

    Parameters
    ----------
    formula : str
        e.g., 'C6H12O6', 'C2H3O2Na'.

    Returns
    -------
    Dict[str, float]
        Element symbol -> number of atoms.

    ASSUMPTION: Supports simple formula notation. Multi-character
    element symbols (e.g., 'Mg') are recognised. Unknown tokens
    are returned as-is with count 1.
    """
    if not formula:
        return {}
    mapping: Dict[str, float] = {}
    i = 0
    while i < len(formula):
        if formula[i].isupper():
            j = i + 1
            while j < len(formula) and formula[j].islower():
                j += 1
            element = formula[i:j]
            num_start = j
            while num_start < len(formula) and formula[num_start].isdigit():
                num_start += 1
            count_str = formula[j:num_start]
            count = float(count_str) if count_str else 1.0
            mapping[element] = mapping.get(element, 0.0) + count
            i = num_start
        else:
            j = i + 1
            while j < len(formula) and not formula[j].isupper() and not formula[j].isdigit():
                j += 1
            if j == i + 1:
                mapping[formula[i]] = mapping.get(formula[i], 0.0) + 1.0
                i += 1
            else:
                element = formula[i:j]
                mapping[element] = mapping.get(element, 0.0) + 1.0
                i = j
    return mapping


def compute_atom_balances(
    reactions: List[Reaction],
    metabolite_formulas: Dict[str, str],
    elements: Optional[List[str]] = None,
) -> Dict[str, List[float]]:
    """
    Compute the atom balance for each element across all reactions.

    Parameters
    ----------
    reactions : list of Reaction
        Reactions to check.
    metabolite_formulas : dict
        Metabolite name -> formula string.
    elements : list of str, optional
        Elements to check. Default: ['C', 'H', 'O', 'N', 'P', 'S'].

    Returns
    -------
    dict
        element -> list of net atom counts (one per reaction).
        Near-zero means the reaction is balanced for that element.
    """
    if elements is None:
        elements = ["C", "H", "O", "N", "P", "S"]
    balances: Dict[str, List[float]] = {e: [] for e in elements}
    parsed: Dict[str, Dict[str, float]] = {}
    for name, formula in metabolite_formulas.items():
        parsed[name] = parse_formula(formula)
    for rxn in reactions:
        net = rxn.net_stoichiometry
        for element in elements:
            net_atoms = 0.0
            for met, coeff in net.items():
                if met in parsed:
                    net_atoms += coeff * parsed[met].get(element, 0.0)
            balances[element].append(net_atoms)
    return balances


def carbon_balance(reactions: List[Reaction],
                   metabolite_formulas: Dict[str, str]) -> float:
    """Compute the overall carbon balance for a set of reactions."""
    atom_bal = compute_atom_balances(reactions, metabolite_formulas, elements=["C"])
    return sum(abs(v) for v in atom_bal["C"])


def charge_balance(reactions: List[Reaction],
                   metabolite_charges: Dict[str, int]) -> Dict[str, float]:
    """Compute the charge balance for each reaction."""
    balances: Dict[str, float] = {}
    for rxn in reactions:
        net = rxn.net_stoichiometry
        net_charge = 0.0
        for met, coeff in net.items():
            net_charge += coeff * metabolite_charges.get(met, 0)
        balances[rxn.name] = net_charge
    return balances


# ---------------------------------------------------------------------------
# ReactionNetwork class
# ---------------------------------------------------------------------------

@dataclass
class ReactionNetwork:
    """
    A metabolic network defined by a set of reactions and metabolites.

    Attributes
    ----------
    name : str
        Network name for identification.
    reactions : list of Reaction
        All reactions in the network.
    metabolites : Dict[str, Metabolite]
        Metabolite name -> Metabolite object.
    external_metabolites : Set[str]
        Metabolites considered external (exchangeable with environment).
    """
    name: str = "unnamed_network"
    reactions: List[Reaction] = field(default_factory=list)
    metabolites: Dict[str, Metabolite] = field(default_factory=dict)
    external_metabolites: Set[str] = field(default_factory=set)

    def _get_internal_metabolites(self) -> List[str]:
        return [
            m for m in self.metabolite_names
            if m not in self.external_metabolites
        ]

    @property
    def metabolite_names(self) -> List[str]:
        names: Set[str] = set()
        for rxn in self.reactions:
            names.update(rxn.all_metabolites)
        return sorted(names)

    def add_reaction(self, reaction: Reaction) -> None:
        self.reactions.append(reaction)

    def remove_reaction(self, name: str) -> None:
        self.reactions = [r for r in self.reactions if r.name != name]

    def add_metabolite(self, metabolite: Metabolite) -> None:
        self.metabolites[metabolite.name] = metabolite

    def verify_mass_balance(
        self,
        metabolite_formulas: Optional[Dict[str, str]] = None,
    ) -> Dict[str, bool]:
        """
        Check mass balance for each reaction (C, H, O, N, P, S).

        Parameters
        ----------
        metabolite_formulas : dict, optional
            Metabolite name -> formula. If None, uses self.metabolites.

        Returns
        -------
        dict
            Reaction name -> True if all elements are balanced.
        """
        if metabolite_formulas is None:
            formulas = {name: m.formula for name, m in self.metabolites.items()}
        else:
            formulas = metabolite_formulas
        elements = ["C", "H", "O", "N", "P", "S"]
        balances = compute_atom_balances(self.reactions, formulas, elements)
        results: Dict[str, bool] = {}
        for i, rxn in enumerate(self.reactions):
            balanced = all(abs(balances[el][i]) < 0.5 for el in elements)
            results[rxn.name] = balanced
        return results

    def verify_charge_balance(
        self,
        metabolite_charges: Optional[Dict[str, int]] = None,
    ) -> Dict[str, bool]:
        """
        Check charge balance for each reaction.

        Parameters
        ----------
        metabolite_charges : dict, optional
            Metabolite name -> charge. If None, uses self.metabolites.

        Returns
        -------
        dict
            Reaction name -> True if charge is balanced.
        """
        if metabolite_charges is None:
            charges = {name: m.charge for name, m in self.metabolites.items()}
        else:
            charges = metabolite_charges
        imbalances = charge_balance(self.reactions, charges)
        return {name: abs(imb) < 0.5 for name, imb in imbalances.items()}

    def stoichiometric_matrix(
        self,
        metabolite_order: Optional[List[str]] = None,
        reaction_order: Optional[List[str]] = None,
    ) -> np.ndarray:
        """
        Construct the stoichiometric matrix S.

        S[i, j] = net stoichiometric coefficient of metabolite i in
        reaction j. Rows = metabolites, Columns = reactions.

        Parameters
        ----------
        metabolite_order : list of str, optional
            Ordering of metabolite rows. Default: sorted order.
        reaction_order : list of str, optional
            Ordering of reaction columns. Default: insertion order.

        Returns
        -------
        np.ndarray
            Stoichiometric matrix, shape (n_metabolites, n_reactions).
        """
        if metabolite_order is None:
            metabolite_order = self.metabolite_names
        if reaction_order is None:
            reaction_order = [r.name for r in self.reactions]

        met_idx = {m: i for i, m in enumerate(metabolite_order)}
        rxn_idx = {r: i for i, r in enumerate(reaction_order)}

        n_mets = len(metabolite_order)
        n_rxns = len(reaction_order)
        S = np.zeros((n_mets, n_rxns), dtype=np.float64)

        rxn_map = {r.name: r for r in self.reactions}
        for rname, rcol in rxn_idx.items():
            if rname not in rxn_map:
                continue
            rxn = rxn_map[rname]
            net = rxn.net_stoichiometry
            for mname, coeff in net.items():
                if mname in met_idx:
                    S[met_idx[mname], rcol] = coeff

        return S

    def nullspace_analysis(
        self,
        tol: float = 1e-10,
    ) -> Dict[str, Any]:
        """
        Compute the nullspace of the stoichiometric matrix to find
        conserved moieties (e.g., total NAD+NADH, total CoA).

        Parameters
        ----------
        tol : float
            Singular value threshold for nullspace computation.

        Returns
        -------
        dict
            Keys: 'nullspace' (ndarray), 'dimension' (int),
                  'conserved_moieties' (list of lists),
                  'use_scipy' (bool).

        The nullspace of S^T (left nullspace) gives conserved moiety
        vectors v such that v^T * S = 0.
        """
        S = self.stoichiometric_matrix()
        S_st = S.T

        if HAS_SCIPY:
            ns = null_space(S_st, rcond=tol)
        else:
            u, s, vh = np.linalg.svd(S_st, full_matrices=True)
            ns = vh[s <= tol].T
            if ns.size == 0:
                ns = np.zeros((S_st.shape[1], 0))

        dim = ns.shape[1]
        moieties: List[List[str]] = []
        met_names = self.metabolite_names
        for col in range(dim):
            vec = ns[:, col]
            participating = [
                met_names[i] for i in range(len(met_names))
                if abs(vec[i]) > tol
            ]
            moieties.append(participating)

        return {
            "nullspace": ns,
            "dimension": dim,
            "conserved_moieties": moieties,
            "use_scipy": HAS_SCIPY,
        }

    def flux_balance_analysis_skeleton(
        self,
        biomass_reaction: Optional[str] = None,
        atp_drain: bool = True,
    ) -> Dict[str, Any]:
        """
        Set up the constraint-based skeleton for flux balance analysis.

        Parameters
        ----------
        biomass_reaction : str, optional
            Name of a biomass reaction if present.
        atp_drain : bool
            Whether to include an ATP drain (ATPM) reaction.

        Returns
        -------
        dict
            Keys: 'S' (stoichiometric matrix), 'reactions',
                  'metabolites', 'objective_idx',
                  'bounds' (list of (lb, ub) tuples).

        ASSUMPTION: Steady state (S*v = 0). All reactions have
        default bounds [0, 1000] for irreversible, [-1000, 1000]
        for reversible. Units are mmol/(gDW*h).
        """
        met_order = self.metabolite_names
        rxn_order = [r.name for r in self.reactions]

        extra_rxns: List[str] = []
        if atp_drain:
            drain_name = "ATP_drain"
            if drain_name not in rxn_order:
                extra_rxns.append(drain_name)
                rxn_order = rxn_order + [drain_name]

        S = self.stoichiometric_matrix(
            metabolite_order=met_order,
            reaction_order=[r for r in rxn_order if r != "ATP_drain"],
        )

        if atp_drain and "ATP_drain" in rxn_order:
            drain_col = np.zeros((S.shape[0], 1))
            for met_name, col_name in [("ATP", -1.0), ("ADP", 1.0), ("Pi", 1.0)]:
                if met_name in met_order:
                    drain_col[met_order.index(met_name), 0] = col_name
            S = np.hstack([S, drain_col])

        bounds: List[Tuple[float, float]] = []
        for rxn in self.reactions:
            if rxn.reversibility:
                bounds.append((-1000.0, 1000.0))
            else:
                bounds.append((0.0, 1000.0))
        if atp_drain:
            bounds.append((0.0, 1000.0))

        n_rxns_total = len(rxn_order)
        objective_idx = rxn_order.index("ATP_drain") if "ATP_drain" in rxn_order else -1
        objective = np.zeros(n_rxns_total)
        if objective_idx >= 0:
            objective[objective_idx] = 1.0

        return {
            "S": S,
            "reactions": rxn_order,
            "metabolites": met_order,
            "objective_idx": objective_idx,
            "objective": objective,
            "bounds": bounds,
        }

    def pathway_enumerate(
        self,
        source_metabolite: str,
        target_metabolite: str,
        max_length: int = 5,
    ) -> List[List[str]]:
        """
        Enumerate simple linear pathways from source to target.

        Performs BFS over the reaction graph.

        Parameters
        ----------
        source_metabolite : str
            Starting metabolite name.
        target_metabolite : str
            Target metabolite name.
        max_length : int
            Maximum number of reactions in the pathway.

        Returns
        -------
        list of list of str
            Each pathway is a list of reaction names.

        ASSUMPTION: Reversible reactions are treated as bidirectional.
        """
        adjacency: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
        for rxn in self.reactions:
            for prod in rxn.products:
                adjacency[prod].append((rxn.name, "out"))
            for reac in rxn.reactants:
                adjacency[reac].append((rxn.name, "in"))
            if rxn.reversibility:
                for reac in rxn.reactants:
                    adjacency[reac].append((rxn.name, "out"))
                for prod in rxn.products:
                    adjacency[prod].append((rxn.name, "in"))

        pathways: List[List[str]] = []
        queue: deque = deque()
        queue.append((source_metabolite, []))

        while queue:
            current_met, path = queue.popleft()
            if len(path) > max_length:
                continue
            if current_met == target_metabolite and path:
                pathways.append(path)
                continue
            for rxn_name, direction in adjacency.get(current_met, []):
                if rxn_name in path:
                    continue
                rxn = next((r for r in self.reactions if r.name == rxn_name), None)
                if rxn is None:
                    continue
                next_mets = rxn.products if direction == "out" else rxn.reactants
                for next_met in next_mets:
                    if next_met != current_met:
                        queue.append((next_met, path + [rxn_name]))

        return pathways

    def atp_yield_from_network(
        self,
        glucose_consumed: float = 1.0,
    ) -> float:
        """
        Calculate the net ATP yield from the reaction network.

        Parameters
        ----------
        glucose_consumed : float
            Number of glucose molecules consumed (default: 1).

        Returns
        -------
        float
            Net ATP produced per glucose consumed.
        """
        net_atp = 0.0
        glucose_net = 0.0
        for rxn in self.reactions:
            net = rxn.net_stoichiometry
            atp_coeff = net.get("ATP", 0.0)
            glucose_coeff = abs(net.get("glucose", 0.0) or net.get("Glucose", 0.0))
            net_atp += atp_coeff
            glucose_net += glucose_coeff

        if glucose_net > 0:
            return net_atp * glucose_consumed / glucose_net
        return net_atp * glucose_consumed

    def carbon_balance_network(self) -> float:
        """Compute overall carbon balance for the entire network."""
        formulas = {name: m.formula for name, m in self.metabolites.items()}
        atom_bal = compute_atom_balances(self.reactions, formulas, elements=["C"])
        return sum(abs(v) for v in atom_bal["C"])

    def summary(self) -> Dict[str, Any]:
        """Return a summary dictionary for the network."""
        return {
            "name": self.name,
            "n_reactions": len(self.reactions),
            "n_metabolites": len(self.metabolite_names),
            "internal_metabolites": self._get_internal_metabolites(),
            "external_metabolites": list(self.external_metabolites),
            "carbon_balance_total": self.carbon_balance_network(),
            "estimated_atp_yield": self.atp_yield_from_network(),
        }


# ---------------------------------------------------------------------------
# Specific network constructors (factory functions)
# ---------------------------------------------------------------------------

def _make_glycolysis_reactions() -> List[Reaction]:
    """Return the standard glycolysis reaction set."""
    return [
        Reaction("HK", {"Glucose": 1, "ATP": 1}, {"G6P": 1, "ADP": 1},
                 delta_g_prime=-16.7, enzyme="hexokinase", reversibility=False),
        Reaction("PGI", {"G6P": 1}, {"F6P": 1},
                 delta_g_prime=1.7, enzyme="phosphoglucose isomerase", reversibility=True),
        Reaction("PFK", {"F6P": 1, "ATP": 1}, {"F16BP": 1, "ADP": 1},
                 delta_g_prime=-14.2, enzyme="phosphofructokinase-1", reversibility=False),
        Reaction("ALDO", {"F16BP": 1}, {"DHAP": 1, "G3P": 1},
                 delta_g_prime=23.8, enzyme="aldolase", reversibility=True),
        Reaction("TIM", {"DHAP": 1}, {"G3P": 1},
                 delta_g_prime=7.5, enzyme="triose phosphate isomerase", reversibility=True),
        Reaction("GAPDH", {"G3P": 1, "NAD": 1, "Pi": 1}, {"BPG": 1, "NADH": 1, "H": 1},
                 delta_g_prime=6.3, enzyme="glyceraldehyde-3-phosphate dehydrogenase",
                 reversibility=True),
        Reaction("PGK", {"BPG": 1, "ADP": 1}, {"3PG": 1, "ATP": 1},
                 delta_g_prime=-18.8, enzyme="phosphoglycerate kinase", reversibility=True),
        Reaction("PGM", {"3PG": 1}, {"2PG": 1},
                 delta_g_prime=4.4, enzyme="phosphoglycerate mutase", reversibility=True),
        Reaction("ENO", {"2PG": 1}, {"PEP": 1, "H2O": 1},
                 delta_g_prime=1.8, enzyme="enolase", reversibility=True),
        Reaction("PYK", {"PEP": 1, "ADP": 1}, {"Pyr": 1, "ATP": 1},
                 delta_g_prime=-31.7, enzyme="pyruvate kinase", reversibility=False),
    ]


def glycolysis_network() -> ReactionNetwork:
    """Construct the glycolysis network."""
    net = ReactionNetwork(name="glycolysis")
    for rxn in _make_glycolysis_reactions():
        net.add_reaction(rxn)

    for name, formula, conc, charge in [
        ("Glucose", "C6H12O6", 5e-3, 0),
        ("G6P", "C6H11O9P", 1e-4, -2),
        ("F6P", "C6H11O9P", 5e-5, -2),
        ("F16BP", "C6H10O12P2", 3e-5, -4),
        ("DHAP", "C3H5O6P", 2e-5, -2),
        ("G3P", "C3H5O6P", 4e-5, -2),
        ("BPG", "C3H4O10P2", 1e-6, -4),
        ("3PG", "C3H4O7P", 1e-4, -3),
        ("2PG", "C3H4O7P", 1e-5, -3),
        ("PEP", "C3H2O6P", 2e-5, -3),
        ("Pyr", "C3H3O3", 5e-4, -1),
        ("ATP", "C10H12N5O13P3", 3e-3, -4),
        ("ADP", "C10H11N5O10P2", 1e-3, -3),
        ("Pi", "HO4P", 5e-3, -2),
        ("NAD", "C21H26N7O14P2", 5e-4, -1),
        ("NADH", "C21H25N7O14P2", 1e-4, -2),
        ("H2O", "H2O", 55.5, 0),
        ("H", "H", 1e-7, 1),
    ]:
        net.add_metabolite(Metabolite(name=name, formula=formula,
                                      concentration_default=conc,
                                      compartment="cytosol", charge=charge))
    net.external_metabolites = {"Glucose", "Pyr", "H2O", "H"}
    return net


def _make_tca_reactions() -> List[Reaction]:
    """Return PDH + TCA cycle reactions (per turn)."""
    return [
        Reaction("PDH", {"Pyr": 1, "CoA": 1, "NAD": 1},
                 {"AcCoA": 1, "NADH": 1, "CO2": 1, "H": 1},
                 delta_g_prime=-33.4, enzyme="pyruvate dehydrogenase",
                 reversibility=False),
        Reaction("CS", {"AcCoA": 1, "OAA": 1, "H2O": 1},
                 {"Cit": 1, "CoA": 1, "H": 1},
                 delta_g_prime=-32.2, enzyme="citrate synthase",
                 reversibility=False),
        Reaction("ACO", {"Cit": 1}, {"IsoCit": 1},
                 delta_g_prime=6.7, enzyme="aconitase", reversibility=True),
        Reaction("IDH", {"IsoCit": 1, "NAD": 1},
                 {"AKG": 1, "NADH": 1, "CO2": 1, "H": 1},
                 delta_g_prime=-21.0, enzyme="isocitrate dehydrogenase",
                 reversibility=False),
        Reaction("AKGDH", {"AKG": 1, "CoA": 1, "NAD": 1},
                 {"SucCoA": 1, "NADH": 1, "CO2": 1, "H": 1},
                 delta_g_prime=-33.5, enzyme="alpha-ketoglutarate dehydrogenase",
                 reversibility=False),
        Reaction("SCS", {"SucCoA": 1, "ADP": 1, "Pi": 1},
                 {"Succ": 1, "ATP": 1, "CoA": 1},
                 delta_g_prime=-3.5, enzyme="succinyl-CoA synthetase",
                 reversibility=True),
        Reaction("SDH", {"Succ": 1, "FAD": 1}, {"Fum": 1, "FADH2": 1},
                 delta_g_prime=0.0, enzyme="succinate dehydrogenase",
                 reversibility=True),
        Reaction("FUM", {"Fum": 1, "H2O": 1}, {"Mal": 1},
                 delta_g_prime=-3.8, enzyme="fumarase", reversibility=True),
        Reaction("MDH", {"Mal": 1, "NAD": 1}, {"OAA": 1, "NADH": 1, "H": 1},
                 delta_g_prime=29.7, enzyme="malate dehydrogenase",
                 reversibility=True),
    ]


def tca_network() -> ReactionNetwork:
    """Construct the PDH + TCA cycle network (per turn)."""
    net = ReactionNetwork(name="PDH_TCA")
    for rxn in _make_tca_reactions():
        net.add_reaction(rxn)

    for name, formula, conc, charge in [
        ("Pyr", "C3H3O3", 5e-4, -1),
        ("AcCoA", "C23H34N7O17P3S", 1e-4, -4),
        ("CoA", "C21H32N7O16P3S", 1e-4, -4),
        ("Cit", "C6H5O7", 3e-4, -3),
        ("IsoCit", "C6H5O7", 1e-4, -3),
        ("AKG", "C5H3O5", 2e-4, -2),
        ("SucCoA", "C25H34N7O19P3S", 1e-5, -4),
        ("Succ", "C4H4O4", 2e-4, -2),
        ("Fum", "C4H2O4", 1e-4, -2),
        ("Mal", "C4H3O5", 2e-4, -2),
        ("OAA", "C4H2O5", 1e-5, -2),
        ("FAD", "C27H29N9O15P2", 1e-4, -2),
        ("FADH2", "C27H31N9O15P2", 5e-5, -2),
        ("CO2", "CO2", 1e-3, 0),
        ("NAD", "C21H26N7O14P2", 5e-4, -1),
        ("NADH", "C21H25N7O14P2", 1e-4, -2),
        ("ATP", "C10H12N5O13P3", 3e-3, -4),
        ("ADP", "C10H11N5O10P2", 1e-3, -3),
        ("Pi", "HO4P", 5e-3, -2),
        ("H2O", "H2O", 55.5, 0),
        ("H", "H", 1e-7, 1),
    ]:
        net.add_metabolite(Metabolite(name=name, formula=formula,
                                      concentration_default=conc,
                                      compartment="matrix", charge=charge))
    net.external_metabolites = {"Pyr", "CO2", "H2O", "H", "OAA"}
    return net


def _make_pcm_reactions() -> List[Reaction]:
    """
    Return the Perfect Coupling Module reactions.

    ASSUMPTION: The PCM couples NADH/FADH2 oxidation to ATP synthesis
    with near-perfect efficiency. The stoichiometry reflects 30 H+ per
    NADH and 3 H+ per ATP.
    """
    return [
        Reaction("PCM_CI",
                 {"NADH": 1, "Q": 1, "H_in": 5},
                 {"NAD": 1, "QH2": 1, "H_out": 5},
                 delta_g_prime=-50.0, enzyme="PCM complex I (30 H+/NADH)",
                 reversibility=False),
        Reaction("PCM_CII",
                 {"FADH2": 1, "Q": 1},
                 {"FAD": 1, "QH2": 1},
                 delta_g_prime=-30.0, enzyme="PCM complex II",
                 reversibility=False),
        Reaction("PCM_CIII",
                 {"QH2": 1, "CytC_ox": 2, "H_in": 2},
                 {"Q": 1, "CytC_red": 2, "H_out": 2},
                 delta_g_prime=-25.0, enzyme="PCM complex III",
                 reversibility=False),
        Reaction("PCM_CIV",
                 {"CytC_red": 2, "O2": 0.5, "H_in": 4},
                 {"CytC_ox": 2, "H2O": 1, "H_out": 4},
                 delta_g_prime=-60.0, enzyme="PCM complex IV",
                 reversibility=False),
        Reaction("PCM_ATPase",
                 {"ADP": 1, "Pi": 1, "H_out": 3},
                 {"ATP": 1, "H2O": 1, "H_in": 3},
                 delta_g_prime=30.5, enzyme="PCM ATP synthase (H+/ATP=3)",
                 reversibility=False),
        Reaction("PCM_leak",
                 {"H_out": 1}, {"H_in": 1},
                 delta_g_prime=0.0, enzyme="proton leak (minimised)",
                 reversibility=True),
    ]


def full_natural_network() -> ReactionNetwork:
    """Construct the full natural respiratory network."""
    net = ReactionNetwork(name="full_natural")
    for rxn in _make_glycolysis_reactions():
        net.add_reaction(rxn)
    for rxn in _make_tca_reactions():
        net.add_reaction(rxn)

    # Natural OXPHOS (10 H+/NADH, 6 H+/FADH2, H+/ATP = 3.7)
    net.add_reaction(Reaction("NADH_OX",
                              {"NADH": 1, "O2": 0.5, "H_in": 1},
                              {"NAD": 1, "H2O": 1, "H_out": 4},
                              delta_g_prime=-60.0,
                              enzyme="natural NADH dehydrogenase",
                              reversibility=False))
    net.add_reaction(Reaction("FADH2_OX",
                              {"FADH2": 1, "O2": 0.5, "H_in": 1},
                              {"FAD": 1, "H2O": 1, "H_out": 6},
                              delta_g_prime=-40.0,
                              enzyme="natural FADH2 oxidase",
                              reversibility=False))
    net.add_reaction(Reaction("NAT_ATPase",
                              {"ADP": 1, "Pi": 1, "H_out": 3.7},
                              {"ATP": 1, "H2O": 1, "H_in": 3.7},
                              delta_g_prime=30.5,
                              enzyme="natural ATP synthase",
                              reversibility=False))

    for name, formula, charge in [
        ("Glucose", "C6H12O6", 0), ("Pyr", "C3H3O3", -1),
        ("ATP", "C10H12N5O13P3", -4), ("ADP", "C10H11N5O10P2", -3),
        ("Pi", "HO4P", -2), ("NAD", "C21H26N7O14P2", -1),
        ("NADH", "C21H25N7O14P2", -2), ("FAD", "C27H29N9O15P2", -2),
        ("FADH2", "C27H31N9O15P2", -2), ("H2O", "H2O", 0),
        ("O2", "O2", 0), ("CO2", "CO2", 0), ("H_in", "H", 1), ("H_out", "H", 1),
    ]:
        net.add_metabolite(Metabolite(name=name, formula=formula,
                                      concentration_default=1e-3,
                                      compartment="mitochondrial", charge=charge))
    net.external_metabolites = {"Glucose", "CO2", "O2", "H2O", "H_in", "H_out"}
    return net


def full_ptr94_network() -> ReactionNetwork:
    """Construct the complete PTR-94 network: glycolysis + PDH + TCA + PCM."""
    net = ReactionNetwork(name="PTR94_full")
    for rxn in _make_glycolysis_reactions():
        net.add_reaction(rxn)
    for rxn in _make_tca_reactions():
        net.add_reaction(rxn)
    for rxn in _make_pcm_reactions():
        net.add_reaction(rxn)

    for name, formula, charge in [
        ("Glucose", "C6H12O6", 0), ("G6P", "C6H11O9P", -2),
        ("F6P", "C6H11O9P", -2), ("F16BP", "C6H10O12P2", -4),
        ("DHAP", "C3H5O6P", -2), ("G3P", "C3H5O6P", -2),
        ("BPG", "C3H4O10P2", -4), ("3PG", "C3H4O7P", -3),
        ("2PG", "C3H4O7P", -3), ("PEP", "C3H2O6P", -3),
        ("Pyr", "C3H3O3", -1), ("AcCoA", "C23H34N7O17P3S", -4),
        ("CoA", "C21H32N7O16P3S", -4), ("Cit", "C6H5O7", -3),
        ("IsoCit", "C6H5O7", -3), ("AKG", "C5H3O5", -2),
        ("SucCoA", "C25H34N7O19P3S", -4), ("Succ", "C4H4O4", -2),
        ("Fum", "C4H2O4", -2), ("Mal", "C4H3O5", -2), ("OAA", "C4H2O5", -2),
        ("ATP", "C10H12N5O13P3", -4), ("ADP", "C10H11N5O10P2", -3),
        ("Pi", "HO4P", -2), ("NAD", "C21H26N7O14P2", -1),
        ("NADH", "C21H25N7O14P2", -2), ("FAD", "C27H29N9O15P2", -2),
        ("FADH2", "C27H31N9O15P2", -2), ("Q", "C6H4O2", 0),
        ("QH2", "C6H6O2", 0), ("CytC_ox", "Fe", 1), ("CytC_red", "Fe", 1),
        ("O2", "O2", 0), ("CO2", "CO2", 0), ("H2O", "H2O", 0),
        ("H_in", "H", 1), ("H_out", "H", 1),
    ]:
        net.add_metabolite(Metabolite(name=name, formula=formula,
                                      concentration_default=1e-3,
                                      compartment="mitochondrial", charge=charge))
    net.external_metabolites = {"Glucose", "O2", "CO2", "H2O", "H_in", "H_out"}
    return net
