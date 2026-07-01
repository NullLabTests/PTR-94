# PCM Document 2: Proton Coupling Stoichiometry — Engineering High H⁺/ATP Ratios

## 1. Current Biological Knowledge

### The F₁F₀ ATP Synthase Rotary Machine

ATP synthase is a ~600 kDa molecular machine in eukaryotes (smaller in bacteria). It consists of two domains:

- **F₁ (catalytic head)**: α₃β₃γδε subunits. The three αβ pairs undergo binding-change mechanism (Boyer, 1997, *Annu Rev Biochem* 66:717-749). Each 120° rotation of the γ-subunit drives one ATP synthesis from ADP + Pi.

- **F₀ (membrane rotor)**: a₁b₂cₙ in bacteria; subunit stoichiometry varies (8-15 c-subunits depending on species). The c-ring rotates relative to the a-subunit, driven by proton flow through half-channels.

The H⁺/ATP ratio equals the number of c-subunits in the ring (cross et al., 1992, *Biochemistry* 31:10501-10506). Each 360° rotation produces 3 ATP, and each 360° requires n protons (one per c-subunit):

$$ \text{H}^+/\text{ATP} = \frac{n}{3} $$

### Known H⁺/ATP Ratios

| Organelle | c-ring stoichiometry | H⁺/ATP | Reference |
|-----------|---------------------|---------|-----------|
| *E. coli* | 10 | 3.33 | Meier et al., 2009, *Science* 325:887-890 |
| Mammalian mitochondria | 8 | 2.67 | Watt et al., 2010, *Nature* 467:1055-1059 |
| Yeast mitochondria | 10 | 3.33 | Stock et al., 1999, *Science* 286:1700-1705 |
| Spinach chloroplast | 14 | 4.67 | Seelert et al., 2000, *Nature* 405:418-419 |
| *Ilyobacter tartaricus* | 11 | 3.67 | Meier et al., 2005, *FEBS Lett* 579:5611-5615 |

The functional H⁺/ATP ratio in mitochondria is debated. Direct measurements suggest 3.7 H⁺/ATP under physiological conditions (Watt et al., 2010), while the c₈ ring structure predicts 8/3 = 2.67. The discrepancy may reflect H⁺ slip or additional H⁺ requirements for transport processes (Petersen et al., 2012, *Biochim Biophys Acta* 1817:814-820).

### The Rotary Mechanism in Detail

Proton flow through F₀: Protons enter from the intermembrane space via an entry half-channel in the a-subunit, neutralize a conserved carboxylate (Glu or Asp) on the c-subunit, causing rotation; protons exit to the matrix after ~360° rotation through an exit half-channel (Junge & Nelson, 2015, *Annu Rev Biochem* 84:631-657).

Torque generation: The electrostatic interaction between the protonated c-subunit carboxylate and the arginine residue in the a-subunit (aR210 in *E. coli*) generates the rotational torque (Elting et al., 2017, *PNAS* 114:15-20). Single-molecule studies measure torque of ~40-50 pN·nm (Sielaff et al., 2016, *EMBO J* 27:2684-2692).

## 2. Known Limitations

### Proton Slip

Proton slip occurs when H⁺ flow through F₀ is not coupled to rotation. This reduces the effective H⁺/ATP ratio. Slip increases at high Δp and accounts for 5-15% of total proton flux in natural systems (Petersen et al., 2012). The slip pathway likely involves transient "short circuits" through the a-c interface.

### Subunit Stoichiometry Tradeoffs

The c-ring size determines H⁺/ATP ratio, but also affects:
- **Torque per proton**: Smaller c-rings (c₈) generate higher torque per H⁺ but require higher Δp
- **Rotational smoothness**: Larger c-rings (c₁₄) give finer angular resolution but lower torque per H⁺
- **Kinetic efficiency**: The relationship between c-ring size and ATP turnover rate is complex (Feniouk et al., 2006, *FEBS Lett* 580:5047-5051)

The observed c₈ ring in mammals may represent an optimization for kinetic efficiency rather than thermodynamic efficiency.

### Transport Overhead

In natural systems, ATP/ADP exchange (ANT), Pi transport, and the Pi carrier consume part of the PMF. The cost is ~0.3-0.5 H⁺ per ATP exported (Groen et al., 1982, *Biochem J* 204:427-435). This reduces the effective ATP yield from ~3.7 to ~3.3 H⁺/ATP for the cytoplasmic ATP pool.

## 3. Unknowns

### Minimum Functional H⁺/ATP Ratio

**Status:** Open research question

What is the minimum number of protons per ATP that can drive rotary catalysis? In principle, the thermodynamic minimum is:

$$ \text{H}^+/\text{ATP} = \frac{\Delta G_{\text{ATP}}}{F\Delta p} $$

where $\Delta G_{\text{ATP}} \approx 50$ kJ/mol under physiological conditions. At $\Delta p = 220$ mV (~21.2 kJ/mol per H⁺), the minimum H⁺/ATP = 50/21.2 ≈ 2.36. However, this does not account for entropic costs or kinetic requirements.

### Zero-Slip Operation

**Status:** Open research question

Can an engineered ATP synthase achieve zero proton slip? This requires perfect gating of the proton half-channels and rigid coupling between proton transfer and rotation. Single-molecule studies suggest slip is fundamentally linked to thermal fluctuations in the a-c interface (Ernst et al., 2012, *FEBS Lett* 586:2892-2896).

### Optimal c-Ring Size

**Status:** Open research question

For the PTR-94 target of 90 ATP from 30 H⁺ (NADH) and 20 H⁺ (FADH₂), what c-ring size gives the optimal balance of efficiency, torque, and kinetics? An H⁺/ATP of 3.0 requires a c₉ ring.

## 4. Engineering Targets

| Parameter | Natural | Engineering Target |
|-----------|---------|-------------------|
| H⁺/ATP ratio | ~3.7 (functional) | 3.0 (c₉ ring) |
| Proton slip | 5-15% | <0.1% |
| Transport overhead | 0.3-0.5 H⁺/ATP | 0 (channeled) |
| Torque per H⁺ | ~12 pN·nm | ~15 pN·nm |
| ATP turnover | ~300 s⁻¹ | >500 s⁻¹ |

**Status:** Engineering hypothesis

## 5. Potential Implementation Ideas

### c₉ Ring Design

A c₉ ring gives exactly H⁺/ATP = 3. This stoichiometry has been observed in the c-ring of *Spirulina platensis* (Pogoryelov et al., 2009, *Nat Struct Mol Biol* 16:742-748) and can be considered a naturally validated design. Each c-subunit would be modified to optimize the carboxylate protonation environment:

- Mutagenesis of the essential glutamate (Glu56 in *E. coli* c-subunit) to tune pKa to ~6.5 for optimal proton capture at the entry half-channel
- Engineering of inter-subunit hydrogen bonding for increased c-ring stability
- Addition of hydrophobic patches at the membrane interface to reduce proton leak

### Stator Engineering for Zero Slip

The a-subunit half-channels must be engineered for perfect gating. The entry half-channel is open to the P-side (high H⁺); the exit half-channel is open to the N-side. The a-subunit arginine (essential for proton exclusion) must be positioned to prevent H⁺ flow past the c-subunit carboxylate without rotation.

Implement a "double gate" mechanism:
- Gate 1: approach channel to the c-ring carboxylate (open only when c-subunit is in the correct rotational state)
- Gate 2: exit channel from c-ring carboxylate (open only after 360° rotation)

Reference: Aksimentiev & Schulten, 2005, *Structure* 13:1277-1288 for computational design of ion channel gating.

### Direct ADP/ATP Channeling

Eliminate transport costs by physically coupling ATP synthase (sATP) to the adenylate translocase or by embedding the ATP binding site within a channeled system that delivers ATP directly to energy-consuming reactions. This can be achieved through:

- Fusion proteins linking ATP synthase β-subunits to ATP-utilizing enzymes
- Synthetic protein scaffolds organizing ATP production and consumption modules
- Lipid nanodomain organization to concentrate ATP synthase and ATP consumers in proximity

**Status:** Engineering hypothesis

## 6. Open Questions

1. **Thermodynamic minimum**: Can H⁺/ATP = 2.5 be achieved with a c₁₅ ring and $\Delta p = 250$ mV, or does the kinetic requirement for torque exceed this boundary?

2. **Slip-free gating**: Is zero proton slip physically achievable, or is there irreducible stochastic uncoupling due to thermal motion at the a-c interface?

3. **Tradeoff surface**: For H⁺/ATP = 3, what is the maximum achievable ATP turnover rate? Is 500 s⁻¹ compatible with zero slip?

4. **Channeling feasibility**: Can ADP/ATP channeling eliminate transport proton cost without constraining diffusion below diffusion-limited on-rates for ATP synthase?

5. **pKa tuning**: What is the optimal pKa for the c-subunit carboxylate to maximize both proton capture efficiency at the entry half-channel and release efficiency at the exit half-channel?
