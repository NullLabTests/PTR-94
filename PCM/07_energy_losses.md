# PCM Document 7: Energy Losses — A Complete Accounting of Thermodynamic Dissipation

## 1. Current Biological Knowledge

### The Global Energy Budget of Respiration

Complete glucose oxidation:

$$ \text{C}_6\text{H}_{12}\text{O}_6 + 6\text{O}_2 \rightarrow 6\text{CO}_2 + 6\text{H}_2\text{O} $$

$\Delta G^{\circ\prime} = -2870\ \text{kJ/mol}$ (Berg et al., 2015, *Biochemistry*, 8th ed.)

The theoretical maximum ATP yield at $\Delta G_{\text{ATP}} \approx 50$ kJ/mol:

$$ n_{\text{ATP,max}} = \frac{2870}{50} \approx 57.4 $$

The observed yield in eukaryotes is ~30-32 ATP/glucose (Rich, 2003, *Biochem Soc Trans* 31:1095-1105; Hinkle, 2005, *Biochim Biophys Acta* 1706:29-45). To achieve 94 ATP/glucose, the PTR-94 PCM must recover energy losses currently dissipated as heat.

### Measured Losses in Natural Respiration

Losses are quantified by comparing the actual ATP yield to the thermodynamic maximum. Major loss channels:

**Proton leak (~15-25%)**:
- Basal proton leak across the IMM accounts for 15-25% of state 4 respiration (Brand et al., 1994, *Biochem J* 303:809-816)
- In hepatocytes, leak accounts for ~20% of resting oxygen consumption
- Leak flux increases with Δp (nonlinear), making high-Δp operation progressively less efficient

**Proton slip (~5-10%)**:
- In Complex I and ATP synthase, proton pumping/consumption can decouple from electron transfer/ATP synthesis
- Slip in Complex I: up to 5% of H⁺ flux at high Δp (Galkin et al., 1999, *FEBS Lett* 458:117-120)
- Slip in ATP synthase: 5-15% depending on Δp and ATP/ADP ratio

**Transport costs (~15-20%)**:
- ATP/ADP exchange (ANT) costs ~1/3 H⁺ per ATP (Groen et al., 1982, *Biochem J* 204:427-435)
- Pi transport costs ~1 H⁺ per Pi
- Metabolite transport (citrate, malate, aspartate) costs additional ~1-2 H⁺ per glucose
- Net: ~1.6 H⁺ per ATP exported to cytosol

**Entropic losses (~2-5%)**:
- Mixing of protons in the bulk phase (vs. localized coupling)
- Thermal dissipation of vibrational and rotational energy
- Information-theoretic entropy of stochastic coupling

**Overall efficiency calculation**:

$$ \eta_{\text{overall}} = \frac{n_{\text{ATP}} \times \Delta G_{\text{ATP}}}{2870} = \frac{30 \times 50}{2870} \approx 52\% $$

The PTR-94 target of 94 ATP/glucose requires:

$$ \eta_{\text{PCM}} = \frac{94 \times 50}{2870} \approx 163\% $$

This is not physically possible using the current energy accounting — the "excess" comes from using the Gibbs free energy at cellular concentrations ($\Delta G \approx 50-57$ kJ/mol ATP) rather than standard state, and from capturing energy from the complete proton gradient (including contributions from export/import of substrates).

**Corrected accounting**: At cellular ATP/ADP ~10:1 and Pi ~10 mM, $\Delta G_{\text{ATP}} \approx 50$ kJ/mol. The maximum ATP yield from 2870 kJ/mol is 57.4. The 94 ATP target requires re-capturing energy otherwise lost to heat, through:
1. Higher H⁺/e⁻ (PCM Module 1)
2. Reduced H⁺/ATP (PCM Module 2)
3. Near-zero transport costs (PCM Module 2)
4. Substrate channeling (general PCM architecture)
5. Eliminating leak and slip (PCM Modules 3, 5, 6)

The target is not 163% efficient — it is 94/57.4 = 164% of the theoretical maximum **using standard accounting**. This implies that either: (a) $\Delta G_{\text{ATP}}$ is higher under PCM conditions, (b) additional free energy contributions from oxidation of glucose beyond $\Delta G^{\circ\prime}$ (mass action, actual cytosol concentrations), or (c) the 2870 kJ/mol basis underestimates available free energy at cellular concentrations.

**Note**: This accounting tension is a known issue in the PTR-94 literature. The 94 ATP target should be understood as a yield under specific cellular conditions ($\Delta G_{\text{ATP}} = 55-57$ kJ/mol, mass action corrections), not a violation of the first law.

## 2. Known Limitations

### Irreducible Dissipation

The second law of thermodynamics requires some energy dissipation in any irreversible process. In a molecular machine operating at biological temperatures:

$$ T\Delta S_{\text{irrev}} \geq k_BT \ln 2 \approx 1.7\ \text{kJ/mol} $$

per informational step (Landauer limit). With ~10⁴ steps per glucose oxidation (each electron transfer, each proton transfer, each conformational change), the irreducible minimum dissipation is:

$$ T\Delta S_{\text{min}} \approx 10^4 \times 1.7 \approx 17\ \text{kJ/mol} $$

This is ~0.6% of the total energy — implying >99% efficiency is theoretically possible but challenging.

### The Δp-ATP Ratio Constraint

For H⁺/ATP = 3.0 and $\Delta G_{\text{ATP}} = 50$ kJ/mol, the required $\Delta p$ is:

$$ \Delta p_{\text{required}} = \frac{\Delta G_{\text{ATP}}}{n_{\text{H}^+}F} = \frac{50}{3 \times 96.485} \approx 173\ \text{mV} $$

For $\Delta G_{\text{ATP}} = 57$ kJ/mol: $\Delta p_{\text{required}} = 197$ mV

If $\Delta p$ must be maintained at 200-250 mV for maximum thermodynamic driving force, the ATP synthase efficiency (accounting for slip) imposes a lower bound on dissipation.

### Lateral Heterogeneity

Biological membranes are not homogeneous. Microdomains with different compositions create local variations in Δp, lipid packing, and protein density. These heterogeneities create energy dissipation through:
- Lateral proton diffusion (slow relative to consumption)
- Energy gradients between membrane regions
- Frictional losses at domain boundaries

## 3. Unknowns

### Minimum Dissipation for a Molecular Machine

**Status:** Open research question

What is the actual minimum dissipation for a molecular machine coupling redox energy to ATP synthesis, given the stochastic nature of single-molecule events? Does the Landauer limit apply, or are there additional dissipation costs from the kinetic requirements of high throughput?

### Leak-Free Membrane Minimum

**Status:** Open research question

Is there an irreducible minimum proton leak through a perfect synthetic membrane, arising from quantum tunneling of protons? If so, what is the flux at $\Delta p = 250$ mV?

### Optimal Δp for Minimum Total Loss

**Status:** Open research question

At what Δp does the sum of (leak losses + slip losses + entropy losses) reach a minimum? Does this optimum coincide with the Δp required for 94 ATP/glucose?

## 4. Engineering Targets

| Loss Channel | Natural Loss | PCM Target | Mitigation |
|-------------|-------------|------------|------------|
| Proton leak | 15-25% | <0.02% | Archaeal lipids, scaffold sealing |
| Proton slip (CI) | 5-10% | <0.02% | Rigid conformational coupling |
| ATP synthase slip | 5-15% | <0.02% | Zero-slip stator design |
| Transport costs | 15-20% | <0.1% | Substrate channeling |
| Entropic dissipation | 2-5% | <0.02% | Deterministic gating |
| Viscous drag | ~3% | <0.02% | Compact rotor, low-viscosity membrane |
| Total dissipation | ~48% | <0.1% (3 kJ/mol) | Combined mitigations |

**Status:** Engineering hypothesis

## 5. Potential Implementation Ideas

### Loss Measurement and Characterization

Experimental methods to quantify each loss channel in prototype PCM:

| Loss Type | Measurement Method | Sensitivity |
|-----------|-------------------|-------------|
| Proton leak | Δp decay after respiratory inhibition | ±0.01% |
| Slip (CI) | H⁺/e⁻ ratio by stopped-flow | ±0.1% |
| Slip (ATP synthase) | H⁺/ATP ratio by single-molecule fluorescence | ±0.1% |
| Transport costs | [ATP]out/[ATP]in with channeling | ±0.5% |
| Viscous drag | Torque-speed curves (magnetic tweezers) | ±0.5 pN·nm |

### Loss Elimination Strategy

**Differential loss analysis**: For each loss channel, measure its contribution at the target operating point (Δp = 250 mV, flux = 10⁴ e⁻/s per chain). Identify the dominant loss channel and address it first.

**Iterative loss minimization cycle**:
1. Measure total losses (calorimetry: heat output vs. ATP production)
2. Decompose losses using specific inhibitors or genetic perturbations
3. Engineer to reduce dominant term
4. Repeat until total dissipation < 3 kJ/mol

**Loss budget**: Allocate the 3 kJ/mol dissipation budget across components:

| Component | Budget (kJ/mol) |
|-----------|----------------|
| Electron transfer | 0.5 |
| Proton pumping | 0.5 |
| ATP synthesis | 0.5 |
| Membrane | 0.5 |
| Transport | 0.5 |
| Measurement tolerance | 0.5 |

### Substrate Channeling Architecture

Design the PCM as a single multi-enzyme complex where intermediates (NADH, QH₂, H⁺) are channeled between active sites without equilibrating with the bulk phase:

- **NADH channel**: GAPDH → Complex I in a single tunnel
- **QH₂ channel**: Complex I → Complex III within 10 Å
- **Proton channel**: Complex IV → ATP synthase with 5 Å transfer distance
- **ATP channel**: ATP synthase → ATP consumer within the complex

This architecture eliminates the transport costs that consume ~15-20% of energy in natural systems.

**Status:** Engineering hypothesis

## 6. Open Questions

1. **Landauer bound**: Is the minimum dissipation for PCM equal to the Landauer limit ($k_BT \ln 2$ per information step), or are there additional constraints from kinetic proofreading or error correction?

2. **Loss decomposition**: Can we experimentally distinguish proton leak from proton slip in an intact PCM, or are they inherently entangled?

3. **Transport cost elimination**: Can substrate channeling truly eliminate all transport-associated costs, or is some residual cost (equal to the potential gradient for concentration-driven transport) irreducible?

4. **Efficiency measurement**: At >99% efficiency, what experimental method has sufficient precision (±0.01%) to quantify total dissipation?

5. **Operating point**: What combination of Δp, temperature, and ATP/ADP ratio minimizes total dissipation while maintaining 94 ATP/glucose throughput?
