# PCM Document 5: Membrane Constraints — Maintaining Integrity Under Extreme Proton Motive Force

## 1. Current Biological Knowledge

### Natural Membrane Composition

The inner mitochondrial membrane (IMM) is among the most protein-rich biological membranes (~70% protein by mass). Its lipid composition is distinctive:

| Lipid Class | IMM (mol%) | Function | Reference |
|-------------|-----------|----------|-----------|
| Phosphatidylcholine | ~40% | Bilayer structure | da Costa et al., 2014, *J Biol Chem* 289:30212-30220 |
| Phosphatidylethanolamine | ~35% | Curvature, fusion | Horvath & Daum, 2013, *Prog Lipid Res* 52:260-272 |
| Cardiolipin (CL) | ~15-20% | ETC stabilization | Paradies et al., 2014, *Biochim Biophys Acta* 1837:1037-1047 |
| Phosphatidylinositol | ~5% | Signaling | |

Cardiolipin is uniquely important for ETC function: CL forms non-bilayer phases that accommodate the conical shapes of integral membrane proteins. Complex I, III, IV, and ATP synthase all require CL for optimal activity (Schlame, 2008, *J Lipid Res* 49:1765-1781).

### Proton Permeability of Natural Membranes

The proton permeability coefficient ($P_{\text{H}^+}$) of lipid bilayers is remarkably low but not zero:

$$ P_{\text{H}^+} \approx 10^{-4}\ \text{cm/s} \text{ (pure lipid)} $$

$$ P_{\text{H}^+} \approx 10^{-2}\ \text{cm/s} \text{ (mitochondrial membrane)} $$

The ~100× increase in biological membranes is due to protein-mediated leak pathways (proteins acting as weak protonophores) and the presence of proton-conducting defects at protein-lipid interfaces (Benz & McLaughlin, 1983, *Biophys J* 41:381-398; Garlid et al., 1996, *J Biol Chem* 271:2797-2803).

The proton leak flux at $\Delta p = 180$ mV is approximately:

$$ J_{\text{leak}} = P_{\text{H}^+} \cdot [\text{H}^+]_{\text{p-side}} \cdot \left(1 - e^{-F\Delta\psi/RT}\right) $$

This accounts for ~10-20% of basal respiration (Brand, 2000, *Exp Gerontol* 35:157-164; Divakaruni & Brand, 2011, *Methods Mol Biol* 810:103-119).

### Membrane Capacitance and Dielectric Properties

The IMM has ~1 $\mu$F/cm² capacitance, typical of biological membranes. The energy stored in the membrane capacitor is:

$$ E_{\text{cap}} = \frac{1}{2} C V^2 \approx 0.5 \times 10^{-6} \times (0.18)^2 \approx 1.6 \times 10^{-8}\ \text{J/cm}^2 $$

This is negligible compared to continuous proton flux, but the membrane must withstand the electric field:

$$ \mathcal{E} = \frac{\Delta\psi}{d} \approx \frac{0.18\ \text{V}}{4\ \text{nm}} \approx 45\ \text{MV/m} $$

For comparison: air dielectric breakdown occurs at ~3 MV/m; typical insulators fail at 10-100 MV/m.

## 2. Known Limitations

### Proton Leak Increases with Δp

The proton conductance ($G_{\text{H}^+}$) of biological membranes is not ohmic — it increases nonlinearly with Δp (Nobes et al., 1992, *Biochem J* 282:219-224). At $\Delta p > 200$ mV, leak rates can increase by an order of magnitude. This is likely due to:

- Field-induced structural defects in the lipid bilayer
- Increased access to proteinaceous leak pathways
- Electrostriction (compression of the membrane under high field)

This intrinsic leak sets a practical maximum for Δp in natural systems.

### Dielectric Breakdown Risk

Above ~300 mV (field strength >75 MV/m), biological membranes undergo irreversible dielectric breakdown (electroporation). Transient pores form, causing massive ion leakage, loss of Δp, and cell death. The IMM operates at ~50-60% of this breakdown threshold — which already limits the maximal PMF.

### Lipid Peroxidation

High electron flux through the ETC produces ROS (see Document 6), which initiate lipid peroxidation chain reactions. Polyunsaturated fatty acids (PUFAs) — abundant in mitochondrial membranes — are especially vulnerable. Peroxidized lipids increase membrane permeability, creating a positive feedback loop (Hauet et al., 2003, *J Biol Chem* 278:4091-4100).

## 3. Unknowns

### Maximum Sustainable Δp

**Status:** Open research question

What is the maximum Δp that can be sustained indefinitely without dielectric breakdown or excessive leak? Synthetic lipid bilayers with archaeal lipids can maintain ~250 mV with low leak (<10% of basal). Can a protein-functionalized membrane achieve similar stability?

### Electrostriction Limits

**Status:** Open research question

At $\Delta\psi = 250$ mV, the electric field (~60 MV/m) compresses the membrane by an estimated 5-15% (electrostriction). Does this compression increase proton leak through the protein-lipid interface? At what field does the compression begin to distort integral membrane protein structure?

### Protein-Lipid Interface Leak

**Status:** Open research question

The major proton leak pathway may occur at the protein-lipid boundary rather than through the bulk lipid. Can engineered lipid-protein interface sealing (via lipidic "O-rings" or hydrophobic matching) reduce this to negligible levels?

## 4. Engineering Targets

| Parameter | Natural | Engineering Target |
|-----------|---------|-------------------|
| Δψ (steady state) | -180 mV | -250 mV |
| ΔpH (matrix alkaline) | 0.5 units | 1.5 units |
| Δp total | ~210 mV | ~340 mV |
| Proton leak (% of flux) | 10-20% | <0.1% |
| Membrane capacitance | 1 μF/cm² | 0.5 μF/cm² (lower) |
| Dielectric strength | ~300 kV/cm | >500 kV/cm |
| Membrane thickness | ~4 nm | 5-6 nm (tunable) |

**Status:** Engineering hypothesis

## 5. Potential Implementation Ideas

### Archaeal Lipid Analogs

Archaeal membranes use tetraether lipids (glycerol dialkyl glycerol tetraethers, GDGTs) that form a monolayer instead of a bilayer. These provide:

- **Extreme proton tightness**: $P_{\text{H}^+}$ is 10-100× lower than conventional bilayers (Komatsu & Chong, 1998, *Cell Biochem Biophys* 35:15-25)
- **High dielectric strength**: GDGT membranes withstand >400 mV without breakdown
- **Thermal stability**: Stable to >100°C

Engineering a synthetic membrane system based on C-20 to C-40 GDGT analogs with tailored headgroups for ETC protein insertion.

### Synthetic Lipid Bilayers with Reinforced Scaffolding

Embed the PCM membrane in a synthetic polymer scaffold:

- **Block copolymer membranes**: Poly(ethylene oxide)-poly(butadiene) (PEO-PBD) polymersomes exhibit 100× lower proton permeability than lipid bilayers (Meier et al., 2014, *Langmuir* 30:11645-11652)
- **Peptide-polymer hybrids**: Amphiphilic peptide-polymer conjugates form mechanically robust membranes with tunable dielectric properties
- **Nanodisc stabilization**: Use MSP1-derived nanodiscs as individual PCM reaction platforms, eliminating lateral proton leak

### Protein-Lipid Interface Sealing

Design ETC proteins with "O-ring" domains — concentric rings of hydrophobic residues at the membrane interface that create a tight seal. Rosetta-based design of interface residues to match membrane thickness and provide a ~5 Å hydrophobic mismatch seal.

### Active Leak Compensation

If zero leak is unattainable, implement a feedback system:

- Leak sensors: pH-sensitive GFP (pHluorin) or voltage-sensitive dyes in the membrane
- Active compensatory pumping: Variant of bacteriorhodopsin or other light-driven pump to compensate residual leak by ~0.1%

**Status:** Engineering hypothesis

## 6. Open Questions

1. **Breakdown threshold**: What is the absolute dielectric breakdown strength of a protein-functionalized lipid bilayer — does it differ from pure lipid bilayers?

2. **Archaeal compatibility**: Can integral membrane proteins (Complex I, ATP synthase) function in GDGT monolayer membranes, given the different hydrophobic thickness and methyl-branched chains?

3. **Scaffold toxicity**: Do polymer scaffolds compatible with low proton leak exhibit any toxic effects on embedded proteins (denaturation, altered conformational dynamics)?

4. **Lateral proton conduction**: At high Δp, does lateral proton diffusion along the membrane surface (the "proton wire" or "proton-collecting antenna" effect) limit local Δp? Can this be mitigated by spatial organization?

5. **Leak baseline**: Is there an irreducible minimum proton leak dictated by quantum proton tunneling through the membrane, or can thermodynamic considerations eliminate it entirely?
