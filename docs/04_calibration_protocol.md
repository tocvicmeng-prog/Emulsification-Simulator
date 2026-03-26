# Calibration Wet-Lab Experiment Protocol

**Version**: 1.0
**Date**: 2026-03-26
**Purpose**: Calibrate simulation constants against experimental data
**Total preparations**: ~30 | **Estimated timeline**: 2–3 weeks
**Estimated materials budget**: ~$1,500–2,500 (excluding capital equipment)

---

## Overview

This protocol calibrates 10 simulation constants across 5 independent studies:

| Study | Constants | Preparations | Priority |
|-------|-----------|-------------|----------|
| 1. Interfacial tension | K_L, Γ∞ | 6 | Highest |
| 2. Chitosan viscosity | [η]_chitosan | 5 | High |
| 3. Viscous breakage | C3 | 6 | High |
| 4. Pore structure | A, α, β (empirical) | 9 | Medium |
| 5. IPN mechanics | η_coupling | 6 | Medium |

**Execute in order** — each study's output feeds the next calibration.

---

## Study 1: Interfacial Tension (K_L, Γ∞)

### Purpose
Calibrate the Szyszkowski-Langmuir adsorption isotherm parameters for Span-80 at the paraffin oil / aqueous polymer phase interface. These control the predicted interfacial tension σ(c_Span80, T), which directly determines droplet breakage in L1.

### Materials

| Item | Specification | CAS / Cat. No. | Supplier | Qty |
|------|--------------|-----------------|----------|-----|
| Liquid paraffin (light) | Ph. Eur. grade, density 0.83–0.86 g/mL | 8012-95-1 | Sigma-Aldrich 76235 | 500 mL |
| Span-80 (sorbitan monooleate) | HLB 4.3, ≥60% oleic acid | 1338-43-8 | Sigma-Aldrich S6760 | 100 mL |
| Agarose | Standard melting, gel point 36±1.5°C | 9012-36-6 | Sigma-Aldrich A9539 | 25 g |
| Chitosan | DD ≥90%, MW 100–300 kDa | 9012-76-4 | Sigma-Aldrich 448877 | 10 g |
| Acetic acid (glacial) | ACS grade ≥99.7% | 64-19-7 | Fisher A38-212 | 100 mL |
| Milli-Q water | 18.2 MΩ·cm | — | In-house | 2 L |

### Equipment
- **Pendant-drop tensiometer** with heated cell (e.g., Krüss DSA100, Biolin Theta, or Dataphysics OCA): must support T up to 95°C
- **Heated syringe** or thermostated capillary for dispensing hot aqueous phase
- **Temperature controller**: ±0.5°C stability at 90°C
- **Analytical balance**: 0.001 g resolution

### Procedure

**Aqueous phase preparation (Day 1):**
1. Dissolve 0.6 g chitosan in 100 mL of 1% v/v acetic acid at 50°C with stirring (2 h).
2. Add 1.4 g agarose to the chitosan solution.
3. Heat to 100°C with stirring for 15 min to dissolve agarose completely.
4. Degas under vacuum for 5 min.
5. Transfer to a heated reservoir maintained at 90°C.

**Oil phase preparation:**
6. Prepare 6 oil samples (50 mL each) with Span-80 at 0, 0.25, 0.5, 1.0, 2.0, 4.0% w/v:
   - 0%: pure paraffin oil
   - 0.25%: 0.125 g Span-80 in 50 mL oil
   - 0.5%: 0.25 g
   - 1.0%: 0.50 g
   - 2.0%: 1.00 g
   - 4.0%: 2.00 g
7. Stir each at 90°C for 10 min to ensure complete dissolution.

**Measurements (Day 2):**
8. Set tensiometer cell temperature to 90.0°C. Allow 30 min equilibration.
9. Fill the cuvette with oil sample #1 (0% Span-80).
10. Form a pendant drop of the hot aqueous polymer phase (~5–10 µL) inside the oil.
11. Record drop shape at 1-second intervals for 300 seconds.
12. Extract equilibrium IFT from the final 60 seconds (plateau).
13. Repeat 3 times per oil concentration. Clean syringe between samples.
14. Repeat steps 9–13 for oil samples #2–#6 in ascending Span-80 order.

### Data Analysis

Fit the Szyszkowski-Langmuir equation to the 6 data points:

```
σ(c) = σ₀ - R·T·Γ∞·ln(1 + K_L·c_mol)
```

where c_mol = c_Span80 / 428.6 × 1000 [mol/m³].

- **σ₀**: from the 0% Span-80 measurement (typically 40–50 mN/m at 90°C)
- **Free parameters**: Γ∞ [mol/m²] and K_L [m³/mol]
- **Initial guesses**: Γ∞ = 3.5×10⁻⁶, K_L = 0.75
- **Software**: Python scipy.optimize.curve_fit or Origin nonlinear fit
- **Goodness of fit**: R² > 0.98 required

### Error Treatment
- **Temperature drift**: ±1°C causes ±0.1 mN/m IFT change. Ensure cell is equilibrated.
- **Drop detachment**: Discard runs where drop detaches before 300 s.
- **Surfactant depletion**: Use fresh oil for each measurement; Span-80 adsorption depletes bulk.
- **Polymer gelation**: If drop gels (T drops below ~40°C), data is invalid. Monitor T.
- **Acceptance**: CV < 5% across 3 replicates for each concentration.

### Expected Results
- σ₀ (no surfactant) ≈ 40–50 mN/m at 90°C
- σ at 2% Span-80 ≈ 3–8 mN/m
- Γ∞ ≈ 2–5 × 10⁻⁶ mol/m²
- K_L ≈ 0.3–2.0 m³/mol

### Safety
- Hot oil at 90°C: use heat-resistant gloves, splatter guard.
- Acetic acid: fume hood for preparation. Eye protection.

---

## Study 2: Chitosan Intrinsic Viscosity ([η]_chitosan)

### Purpose
Determine the intrinsic viscosity of the chitosan used, at the emulsification temperature (90°C), to replace the literature estimate of 800 mL/g.

### Materials
- Chitosan: same lot as Study 1
- Acetic acid solution: 1% v/v (same as Study 1)
- Milli-Q water

### Equipment
- **Rotational rheometer** with concentric cylinder (Couette) geometry (e.g., Anton Paar MCR, TA Instruments DHR, Malvern Kinexus)
- Temperature-controlled Peltier jacket (up to 95°C)
- Solvent trap to prevent evaporation

### Procedure

1. Prepare chitosan solutions at 0, 0.1, 0.25, 0.5, 1.0% w/v in 1% acetic acid (50 mL each).
2. Heat each to 90°C on a hotplate with stirring (30 min).
3. Load the solvent blank (1% acetic acid) into the rheometer at 90°C.
4. Measure steady-shear viscosity at γ̇ = 10, 50, 100, 500 s⁻¹. Record zero-shear plateau.
5. Repeat for each chitosan concentration (ascending order).
6. Between samples: rinse geometry 3× with hot 1% acetic acid.

### Data Analysis

Calculate specific viscosity: η_sp = (η_solution - η_solvent) / η_solvent

Plot η_sp/c vs. c. Fit the Huggins equation:

```
η_sp/c = [η] + k_H·[η]²·c
```

- **Intercept** = [η] (intrinsic viscosity, mL/g)
- **Slope** = k_H·[η]² (Huggins coefficient k_H typically 0.3–0.7)

### Error Treatment
- **Evaporation**: Use solvent trap. Discard if meniscus drops visibly.
- **Shear thinning**: Use the zero-shear plateau, not the shear-thinned value.
- **Temperature equilibration**: Wait 5 min after loading before measuring.
- **Acceptance**: Huggins plot R² > 0.95. If k_H < 0.2 or > 0.8, suspect aggregation.

### Expected Results
- [η] for high-MW, high-DD chitosan at 90°C: 300–1200 mL/g
- The existing estimate (800 mL/g) is within this range

### Safety
- Hot solutions at 90°C. Use insulated sample holders.

---

## Study 3: Viscous Breakage Constant (C3)

### Purpose
Calibrate the Alopaeus breakage kernel viscous resistance parameter C3 by comparing predicted vs. measured d32 under breakage-dominant conditions.

### Materials
- Paraffin oil + 3% Span-80 (from Study 1)
- Aqueous polymer phases at 3 viscosities:
  - **Low**: 2% agarose in water (µ_d ≈ 0.05–0.1 Pa·s at 90°C)
  - **Medium**: 4.2% agarose + 1.8% chitosan (standard formulation, µ_d ≈ 0.15–0.3 Pa·s)
  - **High**: 6% agarose + 2% chitosan (µ_d ≈ 0.5–1.0 Pa·s)
- Hexane (ACS grade, CAS 110-54-3) for washing

### Equipment
- **Rotor-stator homogeniser** (e.g., IKA Ultra-Turrax T25, Silverson L5M): with 25 mm rotor
- **Heated vessel**: 250 mL jacketed beaker with temperature control
- **Laser diffraction particle sizer** (e.g., Malvern Mastersizer 3000, Beckman Coulter LS13320)
- **Rheometer** (from Study 2) for confirming µ_d values
- **Optical microscope** with heated stage (backup sizing method)

### Procedure

**For each of the 3 aqueous phases × 2 RPMs (8,000 and 12,000) = 6 runs:**

1. Heat 200 mL paraffin oil + 3% Span-80 to 90°C in the jacketed beaker.
2. Prepare 2 mL of hot aqueous polymer phase (≥ 85°C).
3. Confirm µ_d of the aqueous phase at 90°C using the rheometer (γ̇ = 100 s⁻¹).
4. Set homogeniser to the target RPM.
5. Pour 2 mL aqueous phase into the oil (1 vol% holdup).
6. Homogenise for exactly 15 seconds. Start timer on contact.
7. Immediately take a 1 mL aliquot from the high-shear zone.
8. Dilute 1:100 in cold paraffin oil (20°C) to quench gelation and freeze the size distribution.
9. Measure droplet size distribution by laser diffraction within 30 minutes.
10. Record d10, d32, d50, d90, and the full distribution curve.

### Data Analysis

For each of the 6 conditions, run the PBE simulation with the same parameters (RPM, µ_d, σ from Study 1, t=15s, phi_d=0.01). Scan C3 from 0 to 0.5 in steps of 0.05. Minimise:

```
χ² = Σᵢ [(d32_sim(C3) - d32_exp)ᵢ / d32_exp,ᵢ]²
```

over all 6 conditions simultaneously.

- **Software**: scipy.optimize.minimize_scalar or grid search
- **Expected C3**: 0.1–0.4 for this viscosity range

### Error Treatment
- **Coalescence during mixing**: Keep holdup at 1% and mix time ≤15 s.
- **Gelation during sizing**: Quench immediately in cold oil. Verify by microscopy that droplets are still liquid spheres.
- **Spatial inhomogeneity**: Sample from the same location (near stator gap) every time.
- **RPM accuracy**: Verify with tachometer. ±100 RPM tolerance.
- **Acceptance**: |d32_sim - d32_exp| / d32_exp < 20% for all 6 conditions with the fitted C3.

### Expected Results
- d32 at 8,000 RPM: 5–20 µm depending on viscosity
- d32 at 12,000 RPM: 3–12 µm
- Higher viscosity → larger d32 (viscous resistance)

### Safety
- Rotating machinery: keep fingers, hair, loose clothing away from the rotor.
- Hot oil splash risk: wear face shield and lab coat.
- Hexane: flammable, use in fume hood with no ignition sources.

---

## Study 4: Pore Structure Coefficients (A, α, β)

### Purpose
Calibrate the empirical pore-size model d_pore = A · c^α · (dT/dt)^β against SEM measurements of agarose gel microstructure.

### Materials
- Agarose (same lot as Study 1)
- Milli-Q water
- Liquid nitrogen (for cryo-fracture)
- Gold/palladium sputtering target (for SEM coating)
- Carbon tape, SEM stubs

### Equipment
- **Programmable cooling block** or Peltier-controlled stage (e.g., Linkam PE120, Echotherm IC22)
- **Cryo-SEM** or conventional **SEM** with cryo-fracture capability (e.g., JEOL JSM-7800F, Hitachi SU8010)
- **Gold sputter coater** (e.g., Leica EM ACE200)
- **Image analysis software**: ImageJ/FIJI with DiameterJ plugin or BoneJ

### Procedure

**Gel preparation (3 concentrations × 3 cooling rates = 9 gels):**

| Gel # | Agarose (% w/v) | Cooling rate (°C/min) |
|-------|-----------------|----------------------|
| 1 | 2.0 | 2 |
| 2 | 2.0 | 10 |
| 3 | 2.0 | 20 |
| 4 | 4.0 | 2 |
| 5 | 4.0 | 10 |
| 6 | 4.0 | 20 |
| 7 | 6.0 | 2 |
| 8 | 6.0 | 10 |
| 9 | 6.0 | 20 |

For each gel:
1. Dissolve agarose in Milli-Q water by heating to 100°C with stirring (15 min).
2. Pour 2 mL into a cylindrical mold (ID 10 mm) pre-heated to 90°C.
3. Place on the programmable cooling block. Set the target cooling rate.
4. Cool from 90°C to 20°C at the specified rate. Record actual T(t) profile.
5. Allow to equilibrate at 20°C for 30 min.

**Cryo-SEM imaging:**
6. Fracture the gel cylinder with a razor blade in liquid nitrogen.
7. Mount on SEM stub with carbon tape, fractured face up.
8. Sputter-coat with 5 nm Au/Pd.
9. Image at 3 magnifications per sample: 10k×, 25k×, 50k×.
10. Acquire at least 5 fields of view per gel at 25k×.

**Image analysis:**
11. Threshold the SEM images (Otsu or manual).
12. Measure pore diameters using DiameterJ or manual measurement of 50+ pores per image.
13. Report: mean pore diameter, standard deviation, distribution histogram.

### Data Analysis

Fit the power law: d_pore = A · c^α · (dT/dt)^β

Take logarithm: ln(d_pore) = ln(A) + α·ln(c) + β·ln(dT/dt)

This is a linear regression in 3 unknowns: ln(A), α, β. Use the 9 data points.

- **Software**: Python numpy.linalg.lstsq or statsmodels OLS
- **Initial guesses**: A ≈ 600 nm, α ≈ -0.7, β ≈ -0.2
- **Goodness of fit**: R² > 0.90 required

### Error Treatment
- **Gel drying during SEM**: Cryo-fracture preserves hydrated structure. Conventional drying destroys pores.
- **Cooling rate deviation**: Log actual T(t) and compute realised rate from the slope in the 60–30°C range.
- **Image selection bias**: Use systematic random sampling (every 3rd field of view).
- **Pore measurement subjectivity**: Use automated analysis (DiameterJ) with manual spot-check.
- **Acceptance**: Power-law fit R² > 0.85. Residuals should not show systematic trend with c or dT/dt.

### Expected Results

| Agarose % | Pore size (2°C/min) | Pore size (10°C/min) | Pore size (20°C/min) |
|-----------|--------------------|--------------------|---------------------|
| 2% | 300–500 nm | 200–350 nm | 150–300 nm |
| 4% | 150–250 nm | 100–200 nm | 80–150 nm |
| 6% | 80–150 nm | 50–100 nm | 40–80 nm |

### Safety
- Liquid nitrogen: cryogenic burns, asphyxiation in confined spaces. Use in ventilated area.
- SEM: trained operator only. Follow facility SOP.

---

## Study 5: IPN Coupling Coefficient (η_coupling)

### Purpose
Calibrate the double-network coupling term η by comparing individual network moduli with the combined DN gel modulus.

### Materials
- Agarose, chitosan (same lots as Studies 1–2)
- Genipin (CAS 6902-77-8, Sigma-Aldrich G4796 or Challenge Bioproducts)
- DMSO (if needed to dissolve genipin, CAS 67-68-5)
- PBS buffer pH 7.4

### Equipment
- **Oscillatory rheometer** with 20 mm parallel plate geometry (e.g., Anton Paar MCR, TA DHR)
- **Water bath** at 37°C for crosslinking
- **Gel casting molds**: 20 mm diameter, 2 mm height (matching rheometer geometry)

### Procedure

**Prepare 6 gels (2 formulations × 3 gel types):**

| Gel # | Formulation | Type | Components |
|-------|-------------|------|------------|
| 1A | Standard (4.2% agar + 1.8% chit) | Agarose-only | 4.2% agarose in water |
| 1B | Standard | Chitosan+genipin | 1.8% chitosan + 2 mM genipin in 1% AcOH |
| 1C | Standard | DN gel | 4.2% agarose + 1.8% chitosan + 2 mM genipin |
| 2A | High-chit (3.5% agar + 2.5% chit) | Agarose-only | 3.5% agarose in water |
| 2B | High-chit | Chitosan+genipin | 2.5% chitosan + 3 mM genipin |
| 2C | High-chit | DN gel | 3.5% agarose + 2.5% chitosan + 3 mM genipin |

For each gel:
1. Prepare the polymer solution as in Study 1.
2. For genipin-containing gels: add genipin stock (100 mM in DMSO) at 50°C.
3. Pour into 20 mm × 2 mm cylindrical molds at 50°C.
4. Cool to 20°C at 10°C/min on the programmable block.
5. For crosslinked gels (B, C): incubate at 37°C for 24 h in sealed container.
6. For agarose-only (A): store at 4°C overnight.

**Rheological measurement:**
7. Load gel disc onto the rheometer parallel plate (20 mm).
8. Set gap to gel thickness (typically 1.8–2.0 mm).
9. Apply 0.1 N normal force to ensure contact.
10. **Strain sweep**: 0.01–10% strain at 1 Hz, 37°C. Identify the linear viscoelastic region (LVR).
11. **Frequency sweep**: 0.1–100 rad/s at strain within LVR (typically 0.5%), 37°C.
12. Record G' (storage modulus) at 1 Hz as the plateau modulus.

### Data Analysis

For each formulation (1 and 2), extract:
- G_agarose = G' of gel A
- G_chitosan = G' of gel B
- G_DN = G' of gel C

Compute η_coupling for each formulation:

```
η = (G_DN - G_agarose - G_chitosan) / √(G_agarose · G_chitosan)
```

Report the mean η across both formulations.

### Error Treatment
- **Gel slippage**: Use sandpaper-coated plates or apply a thin adhesive layer.
- **Dehydration during measurement**: Use a solvent trap or mineral oil edge seal.
- **Strain amplitude**: Verify LVR for each gel type. Over-strain destroys DN structure.
- **Crosslinking variability**: Use identical genipin stock, timing, and temperature.
- **Acceptance**: G' plateau should be frequency-independent (±10%) in the 1–10 rad/s range. η should be consistent between the two formulations (±0.1).

### Expected Results
- G' agarose-only (4.2%): 30–80 kPa
- G' chitosan+genipin (1.8%): 0.5–5 kPa
- G' DN gel: 25–85 kPa (typically slightly less than agarose-only due to swelling constraint)
- η_coupling: -0.3 to 0.0 (antagonistic for sequential IPN without sacrificial bonds)

### Safety
- Genipin: stains skin blue, irritant. Wear gloves. Non-toxic at concentrations used.
- DMSO: penetrates skin and carries dissolved substances through. Double-glove.

---

## Inputting Calibrated Values

After completing all 5 studies, update the simulation constants:

### Option A: Edit MaterialProperties directly (Python)
```python
from emulsim.datatypes import MaterialProperties
props = MaterialProperties(
    sigma=...,                 # From Study 1 (or update K_L/Gamma_inf in interfacial.py)
    # For K_L and Gamma_inf: edit src/emulsim/properties/interfacial.py lines 49-50
    breakage_C3=...,           # From Study 3
    eta_coupling=...,          # From Study 5
)
```

### Option B: Edit the source defaults (permanent)
| Constant | File | Line | Field |
|----------|------|------|-------|
| K_L | `src/emulsim/properties/interfacial.py` | 50 | `K_L = ...` |
| Γ∞ | `src/emulsim/properties/interfacial.py` | 49 | `Gamma_inf = ...` |
| [η]_chitosan | `src/emulsim/properties/viscosity.py` | 169 | `eta_intr_chit = ...` |
| C3 | `src/emulsim/datatypes.py` | 221 | `breakage_C3: float = ...` |
| A (pore) | `src/emulsim/level2_gelation/solver.py` | ~548 | `600e-9` |
| α (pore) | `src/emulsim/level2_gelation/solver.py` | ~548 | `-0.7` |
| β (pore) | `src/emulsim/level2_gelation/solver.py` | ~558 | `-0.2` |
| η_coupling | `src/emulsim/datatypes.py` | 174 | `eta_coupling: float = ...` |

### Option C: TOML configuration (planned)
A future update will allow all calibrated constants to be set via `configs/default.toml`.

---

## Quality Assurance

After calibrating all constants, run a **blind validation**:
1. Choose 1–2 formulations NOT used in any calibration study.
2. Prepare microspheres using the full workflow (emulsification → cooling → crosslinking).
3. Measure d32, pore size (SEM), and G' (rheometry).
4. Compare with simulation predictions.
5. **Acceptance**: all 3 predictions within 25% of measured values.

---

## References

- Zhao et al. (2020) Eng. Life Sci. 20:504–513. DOI: 10.1002/elsc.202000023
- Opawale & Burgess (1998) J. Colloid Interface Sci. 197:142–150
- Butler et al. (2003) J. Polym. Sci. A 41:3941–3953
- Aymard et al. (2001) Biopolymers 59:131–144
- Normand et al. (2000) Biomacromolecules 1:730–738
- Pernodet et al. (1997) Electrophoresis 18:55–58

---

> **Disclaimer**: This protocol is provided for research guidance. All experiments should be
> conducted by trained personnel following institutional safety procedures. The author is an
> AI assistant; validate all procedures through appropriate peer review before implementation.
