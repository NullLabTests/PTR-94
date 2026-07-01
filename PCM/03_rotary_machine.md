# PCM Document 3: The Rotary Machine — Optimizing ATP Synthase for Maximum Efficiency

## 1. Current Biological Knowledge

### Torque Generation and Transmission

ATP synthase is a rotary molecular motor operating in reverse of the flagellar motor. Protons flowing through F₀ generate torque that is transmitted to F₁ via the central stalk (γ-subunit). The torque-generating mechanism involves:

**Electrostatic ratchet**: The conserved arginine in subunit a (aR210 in *E. coli*) creates an electrostatic barrier that only the protonated c-subunit carboxylate can cross. This generates directional bias in the Brownian motion of the c-ring (Junge, 2004, *FEBS Lett* 576:1-7; Elting et al., 2017, *PNAS* 114:15-20).

**Torque measurements by single-molecule methods:**

| Method | Torque (pN·nm) | System | Reference |
|--------|----------------|--------|-----------|
| Magnetic tweezers | ~40 | *E. coli* ATP synthase | Sielaff et al., 2016, *EMBO J* 27:2684-2692 |
| Optical tweezers | ~50 | *Bacillus* PS3 | Okuno et al., 2011, *Nature Commun* 2:326 |
| Fluorescence imaging | ~45 | Mammalian | Nakanishi et al., 2018, *J Biol Chem* 293:13400-13410 |

**Energy per 360° rotation:**

$$ E = 2\pi\tau $$

For $\tau = 45$ pN·nm = $45 \times 10^{-21}$ J, the energy per rotation = $2\pi \times 45 \times 10^{-21} = 2.83 \times 10^{-19}$ J = 170 kJ/mol. This is approximately 3× the $\Delta G_{\text{ATP}}$ ≈ 50-55 kJ/mol, consistent with 3 ATP per rotation.

### Elastic Power Transmission

The central stalk (γ-subunit) and the peripheral stalk (b₂δ in bacteria, OSCP/b/d/F₆ in mammals) act as elastic energy buffers. The 120° stepping of the γ-subunit occurs in three substeps (80°, 40°) with elastic power transmission between rotary elements (Junge et al., 2009, *Nature* 459:364-370; Sielaff et al., 2018, *Nat Commun* 9:3884).

The elastic compliance (bending flexibility) of the central stalk is essential:
- It stores torque during the 40° substeps
- It allows the F₁ head to "wait" for proton-driven rotation in F₀
- It smooths the otherwise discontinuous rotation

### Efficiency of Natural ATP Synthase

Direct calorimetric measurements indicate ATP synthase operates at 80-90% thermodynamic efficiency (Petersen et al., 2012, *Biochim Biophys Acta* 1817:814-820). The remaining 10-20% of energy is dissipated as heat through:

- **Viscous drag** in the membrane (rotating c-ring vs lipid bilayer)
- **Frictional losses** at bearing surfaces (γε interface with αβ hexamer)
- **Slippage** at the a-c interface (proton leak without rotation)
- **Conformational relaxation** in elastic elements

## 2. Known Limitations

### Viscous Drag Losses

The rotating c-ring experiences drag proportional to its radius ($r$), membrane viscosity ($\eta_{\text{mem}} \approx 10^3$ Pa·s), and angular velocity ($\omega$):

$$ \tau_{\text{drag}} \approx 4\pi\eta_{\text{mem}} r^2 d \cdot \omega $$

For $r = 5$ nm, $d = 8$ nm (membrane thickness), $\omega = 2000$ rad/s (100 μm/s tip speed):

$$ \tau_{\text{drag}} \approx 4\pi \times 10^3 \times (5 \times 10^{-9})^2 \times 8 \times 10^{-9} \times 2000 \approx 5 \times 10^{-21}\ \text{N·m} = 5\ \text{pN·nm} $$

This represents ~10% of the total torque, consistent with observed efficiency losses.

### Brownian Ratchet Infidelity

The ATP synthase operates in a Brownian regime where thermal fluctuations are comparable to the energy barriers between rotational states. Each 120° step involves surmounting a ~20-25 $k_BT$ barrier (ATP binding and product release steps). Thermal noise can cause backward steps or slip, reducing net efficiency (Okuno et al., 2011, *Nature Commun* 2:326; Noji et al., 2017, *Elife* 6:e28188).

### Bearing Friction

The α-β-ε interface acts as a molecular bearing. The ε-subunit in bacteria acts as a "brake" at low Δp (Feniouk et al., 2005, *Biochemistry* 44:9849-9862). In the mammalian system, the δ-subunit and associated proteins (IF1) regulate activity. Parasitic friction from these regulatory elements can consume 5-15 pN·nm of torque.

## 3. Unknowns

### Maximum Achievable Mechanical Efficiency

**Status:** Open research question

What is the maximum possible efficiency of a protein-based rotary motor? Structural proteins (actin/myosin) achieve ~40% efficiency in muscle. Kinesin achieves ~60%. ATP synthase at 80-90% is already exceptionally efficient. Can >99% be achieved with engineered components?

### Rigid vs. Elastic Coupling

**Status:** Open research question

For maximum efficiency, should the rotor-stator coupling be perfectly rigid (no elastic losses, but risk of mechanical jamming) or elastic (with energy storage but hysteresis losses)? The optimal design may depend on whether the system is modeled as a "power-stroke" or "Brownian ratchet" device.

### Minimum Bearing Friction

**Status:** Open research question

What is the irreducible frictional loss at the rotor-stator interface? Diamond-like carbon bearings in macroscopic machines achieve friction coefficients $\mu < 0.001$. Are protein-protein interfaces inherently limited to $\mu \approx 0.01-0.1$?

## 4. Engineering Targets

| Parameter | Natural | Engineering Target |
|-----------|---------|-------------------|
| Mechanical efficiency | 80-90% | >99.5% |
| Torque output | 45 pN·nm | 55-60 pN·nm |
| Angular step size | 120° (3 steps) | 120° |
| Viscous drag loss | ~5 pN·nm | <0.3 pN·nm |
| Backward stepping rate | ~0.1% | <0.001% |
| Rotor speed at max ATP | ~300 rpm | ~500 rpm |

**Status:** Engineering hypothesis

## 5. Potential Implementation Ideas

### Low-Friction Bearing Design

Replace the natural αβ hexamer - γ subunit interface with engineered low-friction bearings:

- **Protein surface engineering**: Mutagenesis of the γ-subunit surface to expose hydrophilic, low-adhesion residues at the αβ-γ interface. Glycosylation of the γ-subunit surface creates a hydration layer that reduces friction by 10-100× (similar to mucin-based lubrication; Lee & Spencer, 2005, *Langmuir* 21:10872-10881).

- **Encapsulated rotor design**: Embed the γ-subunit in a ring-shaped chaperonin-like cage (e.g., GroEL-derived) that acts as a roller bearing, converting sliding friction to rolling friction.

### Rigid Power Transmission

Design the stator (peripheral stalk) with increased rigidity using:

- Coiled-coil extensions of the b-subunits with increased hydrophobic packing
- Cross-linking of b-subunits with engineered disulfide bonds or isopeptide bonds (spyTag/spyCatcher; Zakeri et al., 2012, *PNAS* 109:E690-E697)
- Proline mutations to eliminate hinge regions (based on analysis of b-subunit flexibility: Del Rizzo et al., 2006, *Biophys J* 91:62-77)

### Optimized Rotor Hydrodynamics

Minimize viscous drag by reducing c-ring radius while maintaining 9 c-subunits. Natural c-rings have outer diameter ~7 nm. A compacted c₉ ring with optimized subunit interfaces could achieve ~5.5 nm diameter, reducing viscous drag by ~30% ($\propto r^2$).

Rotary speed optimization: The optimal tip speed ($v_{\text{tip}} = \omega r$) balances drag ($\propto v_{\text{tip}}$) against ATP turnover rate. For maximum efficiency, operate at tip speed that minimizes $E_{\text{drag}}/E_{\text{total}}$.

### Synthetic Rotary Coupling

For extreme efficiency, consider a completely synthetic rotary coupler:

- **Carbon nanotube rotor**: A short (5 nm) carbon nanotube functionalized at each end with ATP-binding and H⁺-binding domains. CNT bearings exhibit near-zero friction ($\mu < 0.001$) in MD simulations (Cumings & Zettl, 2000, *Science* 289:602-604).

- **DNA origami stator**: A rigid DNA barrel housing the rotor, with atomically precise bearing surfaces. DNA bearings can achieve friction coefficients of $\mu \approx 10^{-4}$ in controlled environments (Douglas et al., 2009, *Nature* 459:414-418).

**Status:** Engineering hypothesis

## 6. Open Questions

1. **Efficiency ceiling**: What is the absolute maximum efficiency of a molecular machine operating at $k_BT$ energy scales, given the requirements of Landauer's principle (erasure costs)?

2. **Power-efficiency tradeoff**: At what torque and rotation speed does the power-efficiency product peak for the engineered ATP synthase?

3. **Bearing evolution**: Can directed evolution produce ATP synthase variants with <1% frictional loss, or are the natural surfaces already near-optimal?

4. **Elasticity optimization**: Is an elastic or rigid transmission element optimal for the PCM, given the stochastic nature of both proton-driven rotation and ATP synthesis?

5. **Synthetic feasibility**: Can a hybrid protein-synthetic (CNT/DNA) rotary machine interface with biological proton gradients and ATP synthesis?
