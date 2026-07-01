# PCM Document 1: Redox Energy Capture — Extracting Maximum Work from NADH and FADH2 Oxidation

## 1. Current Biological Knowledge

### The Natural Electron Transport Chain

The mitochondrial electron transport chain (ETC) couples the oxidation of NADH and FADH2 to proton translocation across the inner mitochondrial membrane. Complex I (NADH:ubiquinone oxidoreductase, ~1 MDa, 45 subunits in mammals) catalyzes:

$$ \text{NADH} + \text{Q} + 5\text{H}^+_{\text{matrix}} \rightarrow \text{NAD}^+ + \text{QH}_2 + 4\text{H}^+_{\text{IMS}} $$

The stoichiometry of 4 H⁺ pumped per NADH is well-established (Wikström & Hummer, 2012, *PNAS* 109:7449-7450; Jones et al., 2017, *Biochim Biophys Acta* 1858:557-564). Complex III (cytochrome *bc*1 complex) transfers electrons from QH₂ to cytochrome *c* via the Q-cycle, contributing 2 H⁺ per QH₂ oxidized (Trumpower, 1990, *J Biol Chem* 265:11409-11412). Complex IV (cytochrome *c* oxidase) pumps 1 H⁺ per electron (4 H⁺ per O₂ reduced) and uses 4 H⁺ for water formation (Wikström, 1977, *Nature* 266:271-273).

**Summary of natural stoichiometries:**
- NADH → 10 H⁺ pumped (4 CI + 4 CIII + 2 CIV) — textbook value
- FADH₂ → 6 H⁺ pumped (4 CIII + 2 CIV)
- ATP synthesis: ~3.7 H⁺ per ATP (Watt et al., 2010, *Nature* 467:1055-1059)
- Net: ~2.7 ATP per NADH, ~1.6 ATP per FADH₂

### The Mitchell Chemiosmotic Theory

Peter Mitchell's chemiosmotic hypothesis (Mitchell, 1961, *Nature* 191:144-148; Mitchell, 1966, *Biol Rev* 41:445-502) proposed that electron transport generates a proton electrochemical gradient (proton motive force, PMF) that drives ATP synthesis. The PMF has two components:

$$ \Delta p = \Delta\psi - \frac{2.303 RT}{F} \Delta\text{pH} $$

In energized mitochondria, $\Delta\psi \approx -150$ to $-180$ mV and $\Delta\text{pH} \approx 0.5$ units (alkaline matrix), yielding $\Delta p \approx 180$ mV (Nicholls & Ferguson, 2013, *Bioenergetics*, 4th ed., Academic Press).

### Thermodynamic Limit of NADH Oxidation

The complete oxidation of NADH by O₂:

$$ \text{NADH} + \text{H}^+ + \frac{1}{2}\text{O}_2 \rightarrow \text{NAD}^+ + \text{H}_2\text{O} $$

$\Delta G^{\circ\prime} = -220.0$ kJ/mol (Berg et al., 2015, *Biochemistry*, 8th ed., W.H. Freeman).

The energy available for proton pumping per electron pair:

$$ \Delta G_{\text{available}} = nF\Delta E_{\text{total}} = 2 \times 96.485 \times 1.14 = 220\ \text{kJ/mol} $$

where $\Delta E_{\text{total}} = E_{\text{O}_2/\text{H}_2\text{O}} - E_{\text{NAD}^+/\text{NADH}} = 0.82 - (-0.32) = 1.14$ V.

## 2. Known Limitations

### Suboptimal H⁺/e⁻ Ratios

Nature achieves ~2.5 H⁺/e⁻ for NADH oxidation (5 H⁺ per 2 e⁻ at CI, 2 H⁺ per e⁻ at CIII, 1 H⁺ per e⁻ at CIV). The theoretical maximum, limited by the membrane potential and the thermodynamic span, can be estimated:

$$ n_{\text{H}^+,\max} = \frac{\Delta G_{\text{redox}}}{F\Delta p + w} $$

where $w$ accounts for irreversible losses. At $\Delta p = 180$ mV and assuming $w \approx 10$ kJ/mol:

$$ n_{\text{H}^+,\max} \approx \frac{220}{96.485 \times 0.18 + 10} \approx 8.6 $$

This suggests 8-9 H⁺ per 2 e⁻ could be pumped — approximately double the natural yield.

**Status:** Engineering hypothesis

### Structural Constraints on Proton Pumping

Complex I operates through a long-range conformational change mechanism (Sazanov, 2015, *Nat Rev Mol Cell Biol* 16:375-388). The seven transmembrane helices of the NuoL/NuoM/NuoN subunits form a proton translocation channel. Each "quadruple" channel motif can translocate 1 H⁺ per turnover. Natural Complex I contains 3-4 such channels (Efremov & Sazanov, 2011, *Nature* 476:407-411).

To achieve 8-10 H⁺ per NADH, one would need 8-10 independent proton channels, requiring a substantial expansion of the membrane domain.

### Quinone Pool Limitations

Natural ubiquinone (Q₁₀) operates at near-diffusion-limited rates in the membrane. The Q-cycle at Complex III constrains electron bifurcation efficiency. Any increase in electron flux would require increased quinone pool concentrations or engineered quinone analogs with higher diffusion coefficients and more favorable midpoint potentials.

## 3. Unknowns

### Maximum Possible H⁺/e⁻ Ratio

**Status:** Open research question

What is the absolute thermodynamic ceiling for protons pumped per electron, given the proton motive force that must be maintained for ATP synthesis? If $\Delta p = 200$ mV is the maximum sustainable gradient (see Document 5), and each H⁺ requires $F\Delta p \approx 19.3$ kJ/mol, the maximum number of H⁺ from 220 kJ/mol is ~11, but this leaves no energy for kinetics or overcoming activation barriers.

### Protein-Scale Proton Channel Density

**Status:** Open research question

Can a single protein complex accommodate 8-10 independent proton channels without prohibitive structural instability? The membrane domain of Complex I would need to be 2-3× larger than its current ~60 Å width.

### Electron Tunneling Efficiency at High Flux

**Status:** Open research question

At the high electron fluxes required for PTR-94 (near-maximal respiration), internal electron transfer rates must exceed 10⁴ e⁻/s per complex. Does the protein matrix support such rates without significant charge recombination or side reactions?

## 4. Engineering Targets

| Parameter | Natural | Engineering Target |
|-----------|---------|-------------------|
| H⁺/NADH | 10 | 30 |
| H⁺/FADH₂ | 6 | 20 |
| H⁺/e⁻ ratio (CI) | 2 | 8-10 |
| Complex I efficiency | ~75% | >99% |
| Electron transfer rate | 10²-10³ e⁻/s | 10⁴ e⁻/s |
| Quinone midpoint potential | +30 mV | Tunable +30 to +100 mV |

**Status:** Engineering hypothesis

## 5. Potential Implementation Ideas

### Extended Complex I with Multichannel Membrane Domain

Design an engineered NADH dehydrogenase (dCI) with 8-10 modular proton-translocation subunits arranged in tandem. Each module would contain a canonical "E-channel" motif (conserved glutamate-lysine pair, as in *Thermus thermophilus* Nqo8; Baradaran et al., 2013, *Nature* 494:443-448). The modules would be connected by rigid α-helical linkers to ensure conformational coupling.

The total membrane domain would span ~240 Å (8 channels × 30 Å each), requiring a synthetic lipid environment (see Document 5).

### Engineered Quinone Analogs

Natural ubiquinone has limited membrane solubility and diffusion. Synthetic quinones with:
- Shortened isoprenoid tail (Q₂-Q₄) for faster diffusion
- Modified ring substituents to tune $E_m$ by ±50 mV
- Increased polarity for higher effective concentration in the membrane phase

Reference: Rich & Harper, 1990, *FEBS Lett* 269:139-144; Aubert et al., 2002, *Biochim Biophys Acta* 1556:92-102.

### Redox-Balanced Module Architecture

To handle both NADH and FADH₂:
- **Module A**: 10-channel Complex I for NADH (30 H⁺)
- **Module B**: 7-channel Complex I variant for FADH₂ (20 H⁺) — requiring an engineered FADH₂ entry site and modified FMN cofactor

Alternatively, use a single universal reductase that accepts both NADH and FADH₂ with different channel activation stoichiometries.

**Status:** Engineering hypothesis

## 6. Open Questions

1. **Stoichiometry ceiling**: What is the rigorous thermodynamic limit on H⁺/NADH given $\Delta p = -200$ mV, $T = 310$ K, and accounting for the full free energy landscape of electron transfer?

2. **Channel cooperativity**: Can 8-10 proton channels within a single complex operate independently, or do they exhibit cooperativity that reduces the effective stoichiometry?

3. **Kinetic bottleneck**: At 10⁴ e⁻/s throughput, does proton transfer (rather than electron transfer) become rate-limiting, given diffusion-limited H⁺ movement through protein channels (~10⁴-10⁵ s⁻¹)?

4. **Quinone pool capacity**: Can a synthetic membrane support the quinone concentrations required for 90 ATP/glucose throughput without phase separation or aggregation?

5. **Off-pathway electron transfer**: At near-maximal Δp, does reverse electron transfer become significant, reducing net H⁺ pumping efficiency?
