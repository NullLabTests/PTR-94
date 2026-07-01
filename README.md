# PTR-94: Perfect Thermodynamic Respiration

**Achieving the Theoretical Maximum of 94 ATP per Glucose — A Computational Research Platform**

[![Status: Conceptual](https://img.shields.io/badge/status-conceptual-orange?style=flat-square)](https://github.com/NullLabTests/PTR-94)
[![Bioenergetics](https://img.shields.io/badge/focus-bioenergetics-blue?style=flat-square)](https://github.com/NullLabTests/PTR-94)
[![Synthetic Biology](https://img.shields.io/badge/focus-synthetic%20biology-green?style=flat-square)](https://github.com/NullLabTests/PTR-94)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](https://opensource.org/licenses/MIT)
[![CI](https://img.shields.io/badge/CI-pytest-brightgreen?style=flat-square)](.github/workflows/ci.yml)
[![PDF](https://img.shields.io/badge/PDF-manuscript-red?style=flat-square)](paper.pdf)
[![LaTeX](https://img.shields.io/badge/LaTeX-source-blue?style=flat-square)](docs/paper.tex)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue?style=flat-square)](https://www.python.org/)

> **A rigorously designed hypothetical metabolic pathway that captures the full thermodynamic free energy of glucose oxidation (~2870 kJ/mol) into 94 ATP molecules per glucose under standard biochemical conditions — the theoretical maximum.**
>
> This repository provides the complete conceptual design, thermodynamic analysis, stoichiometric verification, computational simulation framework, PCM engineering decomposition, and an arXiv-style manuscript — all as an open-source research platform.

---

## Why This Matters

ATP is the universal energy currency of life. A pathway that triples the ATP yield per glucose molecule would fundamentally reshape biotechnology:

| Domain | Impact |
|--------|--------|
| **Industrial bioproduction** | 2.5–3× yield for ATP-intensive products (amino acids, vitamins, antibiotics, bioplastics) |
| **Cell-free manufacturing** | Continuous-flow bioreactors operating near thermodynamic limits |
| **Synthetic minimal cells** | Energy module with headroom for genetic and metabolic loads |
| **Closed-loop life support** | Reduced resupply for space missions, submarines, isolated habitats |
| **Origins-of-life research** | Defines the hard upper bound on biological energy transduction |
| **Protein design challenge** | Grand challenge: build a membrane complex pumping 30 H⁺ per NADH |

---

## Project Status

**Status:** Conceptual / Hypothetical. No experimental validation yet.

All statements in this repository are explicitly classified:

- **Established biology** — Well-characterized natural processes with literature support
- **Engineering hypothesis** — Proposed design choices that can be tested
- **Computational prediction** — Results from the included simulation framework
- **Open research question** — Unknowns that require further investigation

---

## Repository Structure

```
PTR-94/
├── README.md                  # This file — project overview and documentation
├── LICENSE                    # MIT License
├── paper.pdf                  # Compiled arXiv-style manuscript (6 pages)
│
├── docs/                      # Extended documentation
│   ├── paper.tex              # LaTeX source for manuscript
│   ├── literature_review.md   # Comprehensive review (48 references, 6,000 words)
│   └── research_questions.md  # 15 testable hypotheses across 3 tiers
│
├── PCM/                       # Perfect Coupling Module decomposition
│   ├── 01_redox_capture.md           # Maximum work from NADH/FADH2
│   ├── 02_proton_coupling.md         # High H+/ATP stoichiometry
│   ├── 03_rotary_machine.md          # ATP synthase optimization
│   ├── 04_direct_phosphorylation.md  # Bypassing chemiosmosis
│   ├── 05_membrane_constraints.md    # Membrane integrity under high Δp
│   ├── 06_ros_management.md          # ROS prevention at high flux
│   ├── 07_energy_losses.md           # Complete dissipation accounting
│   └── 08_failure_modes.md           # System-level failure analysis
│
├── simulation/                # Computational research framework
│   ├── __init__.py            # Package exports
│   ├── thermodynamics.py      # Free energies, Nernst, PMF, entropy
│   ├── energy_balance.py      # Step-by-step energy accounting
│   ├── reaction_network.py    # Reaction definitions, mass balance, stoichiometric matrices
│   ├── kinetics.py            # Enzyme kinetics, proton leak, ATP synthase kinetics
│   ├── pareto_optimizer.py    # NSGA-II multi-objective optimization
│   ├── sensitivity_analysis.py # Sobol indices, local sensitivity
│   ├── monte_carlo.py         # LHS sampling, uncertainty quantification
│   ├── stoichiometry_verification.py # Core 94-ATP verification
│   ├── requirements.txt       # Python dependencies
│   ├── experiments/           # Reproducible simulation experiments
│   │   ├── leakage_sensitivity.py    # Proton leak effects on ATP yield
│   │   ├── atp_free_energy.py        # ΔG_ATP scan (30-60 kJ/mol)
│   │   ├── alternative_carriers.py   # Quinone redox potential effects
│   │   ├── membrane_potential.py     # Δψ scan (-100 to -300 mV)
│   │   ├── proton_slip.py            # ATP synthase slip probability
│   │   ├── enzyme_efficiency.py      # kcat/Km bottleneck scanning
│   │   └── temperature_effects.py    # Q10 temperature dependence
│   └── benchmarks/            # Comparison with natural systems
│       ├── natural_mitochondria.py   # Eukaryote: 31 ATP, 32.4% efficiency
│       ├── bacterial_respiration.py  # Prokaryote: 37 ATP, 39.4% efficiency
│       ├── fermentation.py           # Anaerobic: 2 ATP, 2.1% efficiency
│       └── artificial_pathways.py    # 6 synthetic systems comparison
│
├── tests/                     # Comprehensive test suite
│   ├── test_stoichiometry.py        # 33 tests: ATP yields, redox, carbon balance
│   ├── test_energy_conservation.py  # 18 tests: ΔG verification, efficiency
│   ├── test_mass_conservation.py    # 34 tests: C/H/O/N/P/charge balance
│   ├── test_thermodynamics.py       # 30 tests: Nernst, PMF, entropy, bounds
│   ├── test_kinetics.py             # 28 tests: MM, inhibition, proton leak
│   ├── test_optimization.py         # 24 tests: Pareto, bounds, evolution
│   └── test_sensitivity.py          # 25 tests: Sobol, MC, dimensionless coeffs
│
├── visualizations/            # Scientific figure generation
│   ├── sankey.py              # Energy flow Sankey diagrams
│   ├── energy_flow.py         # ATP yield bar charts
│   ├── reaction_graphs.py     # Network graphs + stoichiometric heatmaps
│   ├── pareto_fronts.py       # Multi-objective optimization fronts
│   ├── sensitivity_plots.py   # Tornado and Sobol plots
│   ├── stoichiometric_matrices.py  # Stoichiometric matrix + nullspace
│   └── efficiency_landscapes.py    # 2D/3D efficiency parameter scans
│
└── .github/workflows/         # CI configuration
    └── ci.yml                 # Automated testing on push/PR
```

---

## Quick Start

### Installation

```bash
git clone https://github.com/NullLabTests/PTR-94.git
cd PTR-94
pip install -r simulation/requirements.txt
```

Requires Python 3.10+ with numpy, scipy, matplotlib (optional), and pytest (for tests).

### Run Stoichiometry Verification

```bash
python simulation/stoichiometry_verification.py
```

Verifies that the pathway design sums to exactly 94 ATP with >99.9% thermodynamic efficiency.

### Run All Tests

```bash
python -m pytest tests/ -v
```

167+ tests covering stoichiometry, energy conservation, mass balance, thermodynamics, kinetics, optimization, and sensitivity.

### Run Experiments

```bash
python -c "
from simulation.experiments import *
LeakageSensitivityExperiment().run()
ATPFreeEnergyExperiment().run()
MembranePotentialExperiment().run()
ProtonSlipExperiment().run()
EnzymeEfficiencyExperiment().run()
TemperatureEffectsExperiment().run()
AlternativeCarrierExperiment().run()
"
```

### Run Benchmarks

```bash
python -c "
from simulation.benchmarks import *
NaturalMitochondriaBenchmark().run()
BacterialRespirationBenchmark().run()
FermentationBenchmark().run()
ArtificialPathwaysBenchmark().run()
"
```

### Generate Visualizations

```bash
python -c "
from visualizations.energy_flow import *
energy_flow_chart(save_path='/tmp/ptr94_energy_flow.png')
from visualizations.sankey import *
energy_sankey(save_path='/tmp/ptr94_sankey.png')
from visualizations.pareto_fronts import *
from simulation.pareto_optimizer import *
# See pareto_fronts.py docstring for full example
"
```

---

## Simulation Framework

The simulation package treats the PTR-94 pathway as a computational research platform for exploring the thermodynamic limits of ATP production.

### Core Modules

| Module | Description | Key Functions |
|--------|-------------|---------------|
| `thermodynamics.py` | Free-energy calculations, Nernst equation, PMF, entropy production | `compute_theoretical_max_atp()`, `proton_motive_force()`, `entropy_production()`, `redox_potential()` |
| `energy_balance.py` | Step-by-step energy accounting across all modules | `full_pathway_energy_balance()`, `compare_with_natural()` |
| `reaction_network.py` | Reaction definitions, mass/charge balance, stoichiometric matrices | `full_ptr94_network()`, `stoichiometric_matrix()`, `nullspace_analysis()` |
| `kinetics.py` | Enzyme kinetics, proton leak, ATP synthase, pathway ODEs | `EnzymeKinetics`, `ProtonLeakKinetics`, `ATPSynthaseKinetics`, `PathwayKinetics` |

### Analysis Modules

| Module | Description |
|--------|-------------|
| `pareto_optimizer.py` | NSGA-II genetic algorithm searching over 9 PCM design parameters; `ReactionNetworkEvolution` for novel coupling architectures |
| `sensitivity_analysis.py` | Local sensitivity (OAT) and global sensitivity (Sobol indices) identifying the most influential parameters |
| `monte_carlo.py` | Latin Hypercube Sampling with uncertainty quantification, confidence intervals, and distribution fitting |

### Key Numerical Targets

| Parameter | Value | Status |
|-----------|-------|--------|
| Glucose ΔG°′ (oxidation) | −2870 kJ/mol | Established biology |
| ATP ΔG°′ (synthesis, std.) | +30.5 kJ/mol | Established biology |
| Theoretical maximum ATP | 94.1 | Computational prediction |
| PTR-94 target yield | 94 ATP | Engineering hypothesis |
| PCM H⁺ per NADH | 30 | Engineering hypothesis |
| PCM H⁺ per FADH₂ | 20 | Engineering hypothesis |
| ATP synthase H⁺/ATP | 3.0 | Engineering hypothesis |
| Coupling efficiency target | >99.9% | Engineering hypothesis |
| Eukaryotic benchmark | 30-32 ATP | Established biology |
| Prokaryotic benchmark | 36-38 ATP | Established biology |

---

## Research Roadmap

### Phase 1: In Silico Exploration (current)
- Reaction-network evolution for high-yield coupling architectures
- Sensitivity analysis to identify limiting parameters
- Multi-objective optimization across yield, efficiency, and stability
- Uncertainty quantification via Monte Carlo simulation

### Phase 2: Cell-Free Prototyping (future)
- Reconstitute Modules 1+2 (natural enzymes, well-characterized)
- Artificial liposome/nanodisc assembly with candidate PCM components
- Microfluidic ATP assay for direct yield measurement
- Iterative design-build-test cycles

### Phase 3: Minimal Synthetic Cell (future)
- Port PTR-94 into JCVI-syn3.0 derivatives
- Engineered membranes for high Δp tolerance
- Growth-coupled ATP selection

### Phase 4: Industrial Scale-Up (future)
- High-yield bioproduction strain engineering
- Continuous cell-free bioreactor development
- Commercial evaluation

---

## Limitations, Assumptions, and Known Unknowns

### Explicit Assumptions

Standard biochemical conditions (ΔG°′ at pH 7, 298 K, 1 M concentrations) are used throughout unless stated otherwise. Key assumptions are marked with `# ASSUMPTION:` in all source code.

Major assumptions:
1. Perfect substrate channeling in the PCM eliminates transport overhead
2. Proton leakage can be reduced to near-zero with synthetic membranes
3. ATP synthase can achieve H⁺/ATP = 3 with zero slip
4. Extended proton-pumping stoichiometry (30 H⁺/NADH) is structurally feasible
5. No competing side reactions consume ATP or dissipate energy

### Known Unknowns

1. **Maximum sustainable H⁺/NADH ratio** — Nature uses ~10; is 30 structurally possible?
2. **Irreducible dissipation** — What entropy production is inevitable?
3. **Membrane stability limit** — Maximum Δp before dielectric breakdown
4. **Kinetic feasibility** — Can the pathway sustain sufficient flux?
5. **Regulation and control** — How is the high-yield state maintained against perturbations?
6. **Evolutionary accessibility** — Can such a pathway evolve incrementally?

---

## FAQ

**Q: Is PTR-94 physically possible?**  
A: The thermodynamics are sound — 94 ATP from 2870 kJ/mol at 30.5 kJ/mol per ATP is arithmetic. Whether the engineering is achievable is an open research question. See [research_questions.md](docs/research_questions.md) for testable hypotheses.

**Q: Why 94 and not 94.1?**  
A: The fractional ATP (0.1) represents energy that cannot be captured in integer ATP molecules. PTR-94 targets 94 whole ATP molecules, dissipating the remaining ~3 kJ/mol (~0.1 ATP equivalent) as unavoidable entropic cost.

**Q: How does this compare to natural systems?**  
A: Nature achieves 30-38 ATP (32-40% efficiency). The 3× gap is due to fixed proton stoichiometry, leakage, slip, and evolutionary constraints — not fundamental physics.

**Q: What are the biggest engineering challenges?**  
A: (1) Building a membrane complex that pumps 30 H⁺ per NADH, (2) eliminating proton leakage at high Δp, (3) reducing ATP synthase slip to zero. See the [PCM/](PCM/) directory for detailed breakdowns.

**Q: Has any of this been validated experimentally?**  
A: No. PTR-94 is a conceptual design and computational research platform. All modules, experiments, and benchmarks are in silico.

**Q: How do I contribute?**  
A: See Contributing below. All contributions are credited.

---

## How to Cite

```bibtex
@software{PTR94_2026,
  title        = {{PTR-94}: Perfect Thermodynamic Respiration -- Achieving the
                  Theoretical Maximum of 94 {ATP} per Glucose},
  author       = {{The PTR-94 Conceptual Design Collective}},
  year         = {2026},
  month        = jun,
  publisher    = {GitHub},
  url          = {https://github.com/NullLabTests/PTR-94}
}
```

---

## Contributing

Contributions are welcome in all forms:

- **Thermodynamic/kinetic modeling** — New experiments, refined parameters, additional constraints
- **Protein design** — Sequence proposals for PCM components
- **Literature** — Additions to the literature review, new references
- **Visualization** — Better figures, interactive dashboards
- **Documentation** — Clarifications, translations, tutorials
- **Code** — Bug fixes, optimizations, new features

Open an issue or pull request. All contributors will be credited.

---

## Bibliography

Detailed literature review available in [docs/literature_review.md](docs/literature_review.md) (48 references, 10 sections).

Key foundational works:
- Mitchell, P. (1961). Coupling of phosphorylation to electron and hydrogen transfer by a chemi-osmotic type of mechanism. *Nature*, 191, 144-148.
- Boyer, P. D. (1997). The ATP synthase — a splendid molecular machine. *Annual Review of Biochemistry*, 66, 717-749.
- Abrahams, J. P., Leslie, A. G. W., Lutter, R., & Walker, J. E. (1994). Structure at 2.8 Å resolution of F1-ATPase from bovine heart mitochondria. *Nature*, 370, 621-628.
- Noji, H., Yasuda, R., Yoshida, M., & Kinosita, K. (1997). Direct observation of the rotation of F1-ATPase. *Nature*, 386, 299-302.
- Sazanov, L. A. (2015). A giant molecular proton pump: structure and mechanism of respiratory complex I. *Nature Reviews Molecular Cell Biology*, 16, 375-388.
- Brand, M. D. (2000). Uncoupling to survive? The role of mitochondrial inefficiency in ageing. *Experimental Gerontology*, 35, 811-820.

---

**Repository maintained by the PTR-94 conceptual design collective.**  
Last updated: 2026-06-30

*"The theoretical maximum is not a limit to be accepted, but a target to be engineered."*
