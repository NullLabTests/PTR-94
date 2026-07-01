# PCM Document 4: Direct Redox-Driven Phosphorylation — Bypassing the Proton Gradient

## 1. Current Biological Knowledge

### Substrate-Level Phosphorylation

Nature achieves direct phosphorylation of ADP in two canonical pathways without using a proton gradient:

**Glycolysis** (phosphoglycerate kinase):

$$ 1,3\text{-BPG} + \text{ADP} \rightarrow 3\text{-PG} + \text{ATP} $$

$\Delta G^{\circ\prime} = -18.5$ kJ/mol — thermodynamically favorable due to the high-energy mixed anhydride bond in 1,3-BPG (Bridger & Harris, 1997, in *Textbook of Biochemistry with Clinical Correlations*, Wiley).

**TCA cycle** (succinyl-CoA synthetase):

$$ \text{Succinyl-CoA} + \text{GDP} + \text{P}_i \rightarrow \text{Succinate} + \text{GTP} + \text{CoA} $$

$\Delta G^{\circ\prime} = -3.4$ kJ/mol, utilizing the thioester energy of succinyl-CoA (Joyce & Matthews, 2015, *J Biol Chem* 290:11545-11555).

These represent **direct chemical coupling**: the energy of a high-energy intermediate is physically transferred to the phosphoanhydride bond of ATP/GTP via enzyme-bound phosphorylated intermediates.

### The Chemiosmotic Paradigm

Mitchell's chemiosmotic theory (Mitchell, 1961) proposed that the proton gradient is an obligatory intermediate between redox and phosphorylation. This has been overwhelmingly validated:

- Uncouplers (e.g., DNP, FCCP) abolish ATP synthesis while redox continues
- Artificially imposed Δp drives ATP synthesis in the dark
- All known respiratory ATP synthases are proton-driven

The success of chemiosmosis is so complete that direct redox-to-ATP coupling has received little serious attention since the 1960s.

### Points of Dissent

Some historical and contemporary observations raise questions:

- **Slater's chemical coupling hypothesis** (Slater, 1953, *Nature* 172:975-978) proposed direct electron-transfer-driven phosphorylation, similar to substrate-level phosphorylation. This was largely abandoned after Mitchell, but never definitively disproven.

- **Conformational coupling** (Boyer, 1965, in *Oxidases and Related Redox Systems*, Wiley) proposed that redox energy drives ATP synthase conformational changes directly. Boyer's binding-change mechanism is now accepted, but the source of conformational energy is always the proton gradient.

- **Guiling Qin's micro-rotary engine model**: Direct torque transmission from Complex I rotation to ATP synthase via connecting rotary shafts has been proposed in patent literature but lacks experimental support.

**Status:** Open research question — the possibility of direct coupling is not ruled out by evidence but is considered unlikely by the mainstream (see: Rich, 2003, *Biochem Soc Trans* 31:1095-1105).

## 2. Known Limitations

### Energetic Incompatibility

The redox energy of NADH (220 kJ/mol) and the energy of ATP synthesis (~50-55 kJ/mol) are mismatched in magnitude. Direct coupling would require a 4:1 (NADH-to-ATP) ratio, which is exactly what chemiosmosis achieves via the proton-counting mechanism. Direct coupling would need a mechanical or chemical gearing system.

### No Known Direct-Coupled Enzymes

No enzyme is known that directly transduces NAD(P)H redox energy into ATP synthesis without a chemiosmotic intermediate. While some electron-transferring proteins drive ATP synthesis "directly" (e.g., some ATPase subunits bind nucleotides in redox-sensitive domains), these never produce net ATP synthesis without Δp.

### Evolutionary Considerations

If direct coupling were feasible, evolution (which explored billions of years of sequence space) would likely have discovered it. The universal conservation of chemiosmotic mechanisms suggests it is the optimal — or only — scalable solution for redox-driven phosphorylation.

## 3. Unknowns

### Is Direct Coupling Thermodynamically Feasible at Scale?

**Status:** Open research question

The free energy of NADH oxidation is delivered as a series of small electron-transfer steps (~0.1-0.2 V each). Can this energy be captured and transduced into mechanical work (conformational change) at each step, and then summed to produce ATP? The requirement for temporal integration of multiple electron transfer events makes this challenging.

### Can Redox-Driven Conformational Changes Be Mechanically Coupled to ATP Synthesis?

**Status:** Open research question

Complex I undergoes long-range conformational changes during turnover (Sazanov, 2015). The ~30 Å movement of the NuoL domain could, in principle, drive a mechanical lever. Could a synthetic protein fuse the moving elements of Complex I to the γ-subunit of ATP synthase?

### What Would the Efficiency Be?

**Status:** Open research question

If direct coupling were achieved, would it inherently be more efficient than chemiosmosis? Direct coupling eliminates membrane proton leak, slip, and transport overhead, but introduces mechanical friction, conformational hysteresis, and elastic losses. A complete accounting:

- Chemiosmotic losses: ~15-25% (leak, slip, transport)
- Direct coupling potential losses: ~5-15% (friction, conformational hysteresis)
- Advantage to direct: ~10% potential efficiency gain

## 4. Engineering Targets

| Parameter | Chemiosmotic | Direct Coupling |
|-----------|-------------|-----------------|
| ATP/NADH | ~2.7 | ~4.0 (theoretical) |
| Proton gradient | Required | Not required |
| Membrane complexity | High (intact, pressurized) | Low (structural only) |
| Efficiency (theoretical) | ~65% | ~80% |
| Implementation timeline | Near-term (existing biology) | Far-term (synthetic) |

**Status:** Engineering hypothesis

## 5. Potential Implementation Ideas

### Nanomechanical Redox-to-ATP Coupler

Design a protein complex where:

1. A redox-sensitive domain (derived from Complex I or an engineered ferredoxin) undergoes a ~10 pN·nm conformational change upon electron transfer
2. Multiple (8-12) such domains are arranged radially around a central ATP synthase F₁ head
3. Each redox event rotates the central γ-subunit by 30-45°, accumulating to 120° per ATP

The device would be analogous to a multi-cylinder engine where each electron-transfer event provides a "power stroke."

Reference: Protein-based mechanical transducers exist in nature (e.g., PAS domains, von Willebrand factor A domains; Pervushin et al., 2007, *Biochemistry* 46:12618-12626).

### Hybrid Electrochemical Cell

Embed the ATP synthase F₁ head in a synthetic nanodisk with an integrated electron-transfer chain. Electrons from NADH flow through an engineered pathway to O₂, with the energy of each step captured by piezoelectric elements that physically push the γ-subunit.

- Piezoelectric peptide modules embedded between ETC complexes and the F₁ head
- Each electron-transfer step generates ~0.1-0.3 V piezoelectric potential
- Accumulated potential discharges to drive γ-subunit rotation

### Redox-Activated Rotary Switch

Design a synthetic single-molecule device where:

1. A [2Fe-2S] cluster in a polypeptide linker changes redox state
2. Reduction/oxidation alters the linker length by ~20% (analogous to protein conformational changes in ferredoxin-NADP⁺ reductase; Mora et al., 2017, *Annu Rev Biophys* 46:431-453)
3. Multiple such linkers in series act as a "ratchet" driving unidirectional rotation of the γ-subunit

**Status:** Engineering hypothesis

## 6. Open Questions

1. **Fundamental feasibility**: Can the free energy of electron transfer be directly converted to mechanical work with >50% efficiency in a protein-based system, or are chemiosmotic intermediates thermodynamically required?

2. **Scalability**: If direct coupling is possible at the single-molecule level, can it be scaled to the throughput required for cellular energy demands (10⁶ ATP per second per cell)?

3. **Efficiency comparison**: Under what conditions (Δp, temperature, ATP/ADP ratio) does direct coupling surpass chemiosmosis in overall free energy transduction efficiency?

4. **Evolutionary bypass**: Why did evolution universally adopt chemiosmosis if direct coupling is physically possible? Is this an "frozen accident," a kinetic optimization, or a thermodynamic necessity?

5. **Hybrid architectures**: Can a system that combines elements of direct coupling and chemiosmosis — using the proton gradient for energy storage but direct mechanical links for rapid response — exceed the efficiency of either alone?
