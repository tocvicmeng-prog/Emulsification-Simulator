# Wet-Lab Calibration Protocol for the Emulsification Simulation System

**Document:** 04 -- Calibration Protocol
**Version:** 1.0
**Date:** 2026-03-26
**Author:** Scientific Advisor (computational-experimental liaison)
**Cross-reference:** [01 -- Scientific Advisor Report](01_scientific_advisor_report.md)

---

## Table of Contents

1. [Overview and Study Design](#1-overview-and-study-design)
2. [Study 1 -- Interfacial Tension (K\_L, Gamma\_inf)](#2-study-1----interfacial-tension)
3. [Study 2 -- Chitosan Intrinsic Viscosity (eta\_intr\_chit)](#3-study-2----chitosan-intrinsic-viscosity)
4. [Study 3 -- Viscous Breakage Constant (C3)](#4-study-3----viscous-breakage-constant)
5. [Study 4 -- Pore Structure Empirical Coefficients](#5-study-4----pore-structure-empirical-coefficients)
6. [Study 5 -- IPN Coupling Coefficient (eta\_coupling)](#6-study-5----ipn-coupling-coefficient)
7. [Inputting Calibrated Values into the Simulation](#7-inputting-calibrated-values-into-the-simulation)
8. [Appendix A -- Reagent Master List](#appendix-a----reagent-master-list)
9. [Appendix B -- Equipment Summary](#appendix-b----equipment-summary)

---

## 1. Overview and Study Design

### 1.1 Objective

This protocol defines five wet-lab calibration studies that supply the material-specific constants required by the emulsification simulation system. Each study targets one or more parameters that cannot be obtained from the literature with sufficient accuracy for the specific material grades and conditions used in this project (4.2% w/v agarose + 1.8% w/v chitosan dispersed in Span-80/paraffin oil at 90 degC). The calibrated values feed directly into the `MaterialProperties` dataclass (`src/emulsim/datatypes.py`) and the TOML property database (`data/properties.toml`).

### 1.2 Summary of Calibration Targets

| Study | Parameters | Simulation Field(s) | Count of Preparations |
|-------|-----------|---------------------|-----------------------|
| 1. Interfacial Tension | K\_L, Gamma\_inf | `interfacial.K_L`, `interfacial.gamma_inf` | 6 (+ triplicates = 18 measurements) |
| 2. Chitosan Intrinsic Viscosity | [eta]\_chit, k\_H | `eta_intr_chit` in `viscosity.py` L169 | 5 solutions |
| 3. Viscous Breakage Constant | C3 | `MaterialProperties.breakage_C3` | 6 emulsification runs |
| 4. Pore Structure | A, alpha, beta in d\_pore = A * c^alpha * (dT/dt)^beta | Empirical pore model coefficients | 9 gels |
| 5. IPN Coupling | eta\_coupling | `MaterialProperties.eta_coupling` | 6 gels (3 triplets x 2 formulations) |

**Total preparations:** ~30 primary samples (plus replicates as specified per study).

### 1.3 Recommended Execution Order

The studies have partial dependencies and should be executed in the order listed:

1. **Study 1** (Interfacial Tension) -- needed before Study 3 (breakage fitting requires sigma).
2. **Study 2** (Chitosan Viscosity) -- needed before Study 3 (breakage fitting requires mu\_d).
3. **Study 3** (Breakage Constant) -- depends on Studies 1 and 2.
4. **Study 4** (Pore Structure) -- independent; can run in parallel with Studies 1-3.
5. **Study 5** (IPN Coupling) -- depends on crosslinking protocol; can run in parallel with Study 3.

### 1.4 Timeline Estimate

| Phase | Duration | Notes |
|-------|----------|-------|
| Reagent procurement | 2-3 weeks | Long-lead items: genipin, high-purity agarose |
| Study 1 | 3-4 days | Including equilibration time |
| Study 2 | 2 days | Rheometer time |
| Study 3 | 3-4 days | Emulsification + sizing |
| Study 4 | 2-3 weeks | Gel preparation (1-2 days) + cryo-SEM queue (1-2 weeks) |
| Study 5 | 5-7 days | Gel curing (24-48 h) + rheology (1 day) |
| Data analysis | 3-5 days | Fitting, cross-validation |
| **Total** | **~5-6 weeks** | With parallel execution of independent studies |

### 1.5 Budget Estimate

| Category | Estimated Cost (USD) |
|----------|---------------------|
| Reagents and consumables | $800 -- $1,200 |
| Instrument time (pendant drop, rheometer, particle sizer, cryo-SEM) | $2,000 -- $4,000 |
| Disposables (cuvettes, syringes, SEM stubs) | $200 -- $400 |
| **Total** | **$3,000 -- $5,600** |

Costs assume institutional access to major instruments. External cryo-SEM service can add $1,500-$3,000.

### 1.6 General Laboratory Practices

- All aqueous solutions should be prepared with ultrapure water (resistivity >= 18.2 MOhm cm).
- Glassware must be cleaned with chromic acid or piranha solution, rinsed exhaustively, and dried before interfacial tension work.
- Record ambient temperature, relative humidity, and barometric pressure at the start of each session.
- Maintain a laboratory notebook with cross-references to this protocol by study and step number.

---

## 2. Study 1 -- Interfacial Tension

### 2.1 Purpose

Calibrate the **Szyszkowski-Langmuir** isotherm parameters **K\_L** (Langmuir adsorption constant, m^3/mol) and **Gamma\_inf** (maximum surface excess concentration, mol/m^2) for Span-80 at the paraffin oil / aqueous polymer phase interface at 90 degC. These two parameters control the interfacial tension model in `src/emulsim/properties/interfacial.py` (function `interfacial_tension_span80`), which is the primary restoring force in the population balance equation for droplet breakup.

The model equation is:

    sigma = sigma_0(T) - R*T*Gamma_inf * ln(1 + K_L * c_mol)

where sigma\_0(T) is the bare (surfactant-free) interfacial tension, R = 8.314 J/(mol K), T is temperature in K, and c\_mol is the molar concentration of Span-80 in the oil phase.

### 2.2 Materials

| Material | Grade / Specification | CAS Number | Supplier (example) | Quantity |
|----------|----------------------|------------|---------------------|----------|
| Light liquid paraffin | Ph. Eur. / USP, kinematic viscosity 25-80 mm^2/s at 40 degC | 8012-95-1 | Sigma-Aldrich (76235) or Merck | 500 mL |
| Span-80 (sorbitan monooleate) | >= 60% oleic acid, HLB 4.3 | 1338-43-8 | Sigma-Aldrich (S6760) | 50 g |
| Agarose | Low-EEO, gel point 36 +/- 1.5 degC, standard melting | 9012-36-6 | Sigma-Aldrich (A9539) or Bio-Rad | 25 g |
| Chitosan | MW 100-300 kDa, degree of deacetylation >= 90% | 9012-76-4 | Sigma-Aldrich (448877) | 10 g |
| Glacial acetic acid | ACS reagent grade, >= 99.7% | 64-19-7 | Fisher Scientific | 100 mL |
| Ultrapure water | >= 18.2 MOhm cm | 7732-18-5 | In-house Milli-Q | 1 L |

### 2.3 Equipment

| Instrument | Recommended Model | Key Specification |
|------------|-------------------|-------------------|
| Pendant-drop tensiometer | Kruss DSA100 or Biolin Scientific Theta Flex | Temperature-controlled measurement cell to 100 degC, automated drop shape analysis, image resolution >= 1024x768, frame rate >= 25 fps |
| Temperature-controlled cell | Kruss TC40 or equivalent Peltier/circulating bath cell | Temperature stability +/- 0.5 degC at 90 degC, optical access on two sides |
| Analytical balance | Mettler Toledo XSR205 or equivalent | Readability 0.01 mg for Span-80 weighing |
| Hotplate/stirrer | IKA RCT basic or equivalent | Temperature up to 100 degC for agarose dissolution |
| Gastight syringes | Hamilton 500 uL, blunt-tip needle | For pendant-drop formation in heated cell |

### 2.4 Procedure

**A. Preparation of Oil Phase Solutions (Day 1)**

1. Weigh the following amounts of Span-80 into six clean 50 mL glass vials:
   - Vial 1: 0 mg (pure paraffin oil control)
   - Vial 2: 62.5 mg (0.25% w/v in 25 mL)
   - Vial 3: 125.0 mg (0.5% w/v in 25 mL)
   - Vial 4: 250.0 mg (1.0% w/v in 25 mL)
   - Vial 5: 500.0 mg (2.0% w/v in 25 mL)
   - Vial 6: 1000.0 mg (4.0% w/v in 25 mL)
2. Add 25.0 mL of liquid paraffin to each vial using a graduated cylinder.
3. Seal vials with PTFE-lined caps and place on an orbital shaker at 40 degC for 2 hours to ensure complete dissolution. Span-80 is oil-soluble and should dissolve readily; verify visual clarity.
4. Allow solutions to equilibrate to room temperature overnight, then store at room temperature protected from light.

**B. Preparation of Aqueous Polymer Phase (Day 1)**

5. Prepare 1% v/v acetic acid solution: add 5.0 mL glacial acetic acid to 495 mL ultrapure water.
6. Dissolve 1.8 g chitosan in 100 mL of the 1% acetic acid solution by stirring at room temperature for 4-6 hours until a clear, slightly viscous solution is obtained. If necessary, warm to 40 degC to accelerate dissolution. Filter through a 100 um mesh to remove undissolved particles.
7. In a separate vessel, add 4.2 g agarose to 100 mL ultrapure water. Heat on a hotplate with stirring to 95-100 degC until the agarose fully dissolves (solution becomes clear).
8. Immediately combine the hot agarose solution with the chitosan/acetic acid solution while the agarose is still above 85 degC. Mix vigorously for 30 seconds. The final composition is 4.2% w/v agarose + 1.8% w/v chitosan in approximately 0.5% acetic acid.
9. Keep the combined solution in a water bath at 90 degC (+/- 2 degC) until measurements begin. Do not allow to cool below 70 degC or agarose will begin to gel.

**C. Pendant-Drop Measurements (Day 2)**

10. Pre-heat the tensiometer measurement cell to 90.0 +/- 0.5 degC. Allow 30 minutes for thermal equilibration.
11. Fill the optical cuvette with the aqueous polymer phase (maintained at 90 degC). The aqueous phase is the bulk (outer) phase.
12. Load the first oil-phase solution (Vial 1, pure paraffin) into the gastight syringe.
13. Form an inverted pendant drop of oil (approximately 5-10 uL) at the tip of a J-shaped needle submerged in the aqueous phase. The drop hangs upward because rho\_oil < rho\_aq.
14. Wait for the drop shape to stabilise (typically 30-120 seconds for surfactant-free oil, up to 300-600 seconds for higher Span-80 concentrations). Monitor the real-time IFT value -- equilibrium is reached when IFT variation is < 0.1 mN/m over 60 seconds.
15. Record the equilibrium IFT value. Capture at least 10 consecutive frames at equilibrium for averaging.
16. Repeat steps 12-15 for a minimum of **3 independent drops** per oil-phase solution (fresh drop each time).
17. Flush the syringe with the next oil-phase solution (Vial 2), discard the first drop, then repeat steps 13-16.
18. Continue through all six oil-phase concentrations in order of increasing Span-80 concentration to minimise cross-contamination.
19. After completing all six concentrations, repeat the 0% and 2% w/v measurements as internal consistency checks.

### 2.5 Measurements

| Parameter | Value to Record | Replicates | Data Format |
|-----------|----------------|------------|-------------|
| Equilibrium IFT (sigma\_eq) | mN/m, 3 significant figures | 3 drops per concentration | CSV: `c_span80_wv, c_span80_mol_m3, sigma_mN_m, sigma_sd, T_K, drop_volume_uL` |
| Temperature | degC, +/- 0.5 degC | Continuous logging | Appended to CSV |
| Drop volume | uL | Each drop | Appended to CSV |
| Dynamic IFT trace | mN/m vs time | At least 1 full trace per concentration | Separate CSV per run |

**Unit conversions for fitting:**
- c\_span80 in mol/m^3 = (c\_wv in g/mL) / (M\_Span80 in g/mol) * 1e6, where M\_Span80 = 428.6 g/mol
- sigma in N/m = sigma\_mN\_m * 1e-3

### 2.6 Data Analysis

**Fitting procedure:**

1. Tabulate the six equilibrium IFT values (mean +/- SD) against Span-80 molar concentration c\_mol (mol/m^3).
2. The bare IFT sigma\_0 is directly measured from the 0% Span-80 data point. Compare to the model prediction `sigma_bare(363.15)` in `interfacial.py`, which gives approximately 0.043 N/m. If the measured value differs by > 20%, update the `sigma_bare` function parameters.
3. Fit the Szyszkowski-Langmuir equation to the remaining 5 data points by nonlinear least squares:

       sigma(c) = sigma_0 - R * T * Gamma_inf * ln(1 + K_L * c)

   where sigma\_0 is fixed at the measured value, T = 363.15 K, R = 8.314 J/(mol K).

4. **Software:** Python `scipy.optimize.curve_fit` or MATLAB `lsqcurvefit`.
5. **Initial guesses:** Gamma\_inf = 3.5e-6 mol/m^2, K\_L = 0.75 m^3/mol (current simulation defaults).
6. **Bounds:** Gamma\_inf in [1e-7, 1e-4] mol/m^2; K\_L in [0.01, 100] m^3/mol.
7. Report the best-fit values with 95% confidence intervals from the covariance matrix.
8. Compute R^2 and the root-mean-square residual. R^2 should be > 0.98 for a well-behaved Langmuir isotherm.
9. Plot: measured sigma vs c\_span80 with fitted curve and 95% prediction bands.

**Sample Python code:**

```python
import numpy as np
from scipy.optimize import curve_fit

R = 8.314
T = 363.15

def szyszkowski(c_mol, Gamma_inf, K_L):
    return sigma_0 - R * T * Gamma_inf * np.log(1.0 + K_L * c_mol)

popt, pcov = curve_fit(szyszkowski, c_data, sigma_data,
                       p0=[3.5e-6, 0.75],
                       bounds=([1e-7, 0.01], [1e-4, 100.0]))
Gamma_inf_fit, K_L_fit = popt
```

### 2.7 Error Treatment

| Source of Error | Mitigation | Acceptance Criterion |
|----------------|------------|---------------------|
| Drop detachment / instability | Use J-needle with controlled dispensing rate; discard drops that detach prematurely | Drop must be stable >= 60 s |
| Temperature gradient in cell | Verify with thermocouple probe inside cell; allow 30 min equilibration | T within +/- 0.5 degC of setpoint |
| Agarose gelation during measurement | Maintain aqueous phase at >= 85 degC; complete measurements within 2 hours of preparation | Visual clarity of solution; no turbidity |
| Span-80 concentration error | Weigh on analytical balance (0.01 mg); prepare volumetrically | Mass error < 1% of target |
| Cross-contamination between concentrations | Flush syringe with 3 volumes of next solution; work low-to-high concentration | Repeat 0% check at end should match initial within 5% |
| Adsorption kinetics (non-equilibrium) | Wait until IFT drift < 0.1 mN/m per 60 s | Dynamic trace must plateau |

**When to repeat:** If the coefficient of variation across triplicates exceeds 10% for any concentration, prepare fresh solutions and repeat that concentration. If the sigma\_0 end-check differs from the initial by > 5%, repeat the entire series.

### 2.8 Safety

- **Liquid paraffin:** Low hazard. Avoid inhalation of mist. Wear safety glasses.
- **Span-80:** Low toxicity. Skin and eye irritant. Wear nitrile gloves.
- **Glacial acetic acid:** Corrosive (GHS05). Handle in fume hood. Wear nitrile gloves, safety glasses, lab coat. Keep away from heat sources.
- **Hot solutions (90 degC):** Burn hazard. Use heat-resistant gloves when handling. Ensure tensiometer cell has splash guards.
- **General:** Standard laboratory PPE at all times. Dispose of oil/surfactant mixtures as organic waste per institutional guidelines.

### 2.9 Expected Results

- **sigma\_0 (bare, 90 degC):** 35-45 mN/m. The polymer solution may reduce this by 5-15% relative to pure water/paraffin (expect 30-43 mN/m).
- **sigma at 2% Span-80:** 3-8 mN/m. The current simulation default is 5 mN/m.
- **sigma at 4% Span-80:** 2-5 mN/m (approaching plateau above CMC).
- **Gamma\_inf:** Expected 2e-6 to 8e-6 mol/m^2 (current default 3.5e-6).
- **K\_L:** Expected 0.1 to 5.0 m^3/mol (current default 0.75).

**Red flags:**
- sigma\_0 < 20 mN/m would suggest contamination of the paraffin oil or glassware.
- sigma not decreasing monotonically with Span-80 concentration: check solution preparation.
- sigma at 4% Span-80 > 15 mN/m: Span-80 may have degraded (check expiry, storage conditions).
- Very poor fit (R^2 < 0.9): consider whether a Frumkin isotherm (with lateral interactions) is more appropriate.

---

## 3. Study 2 -- Chitosan Intrinsic Viscosity

### 3.1 Purpose

Calibrate the **intrinsic viscosity** [eta]\_chit (mL/g) and **Huggins coefficient** k\_H for chitosan in 1% acetic acid at the emulsification temperature (90 degC). These values determine the dispersed-phase viscosity through the Huggins equation in `src/emulsim/properties/viscosity.py`. The current simulation uses a hardcoded estimate of [eta]\_chit = 800 mL/g (line 169), which is appropriate for a ~300 kDa chitosan at 25 degC but has not been verified at 90 degC or for the specific grade used.

The Huggins equation is:

    eta_sp / c = [eta] + k_H * [eta]^2 * c

where eta\_sp = (eta\_solution / eta\_solvent) - 1 is the specific viscosity, c is the polymer concentration in g/mL, [eta] is the intrinsic viscosity, and k\_H is the Huggins constant.

### 3.2 Materials

| Material | Grade / Specification | CAS Number | Supplier (example) | Quantity |
|----------|----------------------|------------|---------------------|----------|
| Chitosan | MW 100-300 kDa, DD >= 90%, same lot as Study 1 | 9012-76-4 | Sigma-Aldrich (448877) | 5 g |
| Glacial acetic acid | ACS reagent, >= 99.7% | 64-19-7 | Fisher Scientific | 50 mL |
| Sodium chloride | ACS reagent, >= 99.0% | 7647-14-5 | Sigma-Aldrich (S9888) | 10 g |
| Ultrapure water | >= 18.2 MOhm cm | 7732-18-5 | In-house Milli-Q | 500 mL |

**Note:** Intrinsic viscosity of polyelectrolytes depends strongly on ionic strength. Use 0.2 M NaCl + 1% acetic acid as the solvent to screen electrostatic interactions and obtain a reproducible [eta]. This is standard practice for chitosan characterisation (Rinaudo et al. 1993).

### 3.3 Equipment

| Instrument | Recommended Model | Key Specification |
|------------|-------------------|-------------------|
| Rotational rheometer | TA Instruments DHR-3, Anton Paar MCR 302/502 | Torque resolution <= 1 nN m for low-viscosity solutions |
| Concentric cylinder geometry | DIN standard, bob diameter ~25 mm, cup diameter ~27 mm | Double-gap geometry preferred for low viscosity |
| Temperature control | Peltier-heated cup or circulating bath | Stability +/- 0.1 degC at 90 degC, range to 100 degC |
| Solvent trap / evaporation cover | Supplied with geometry | Essential at 90 degC to prevent evaporation |
| 0.45 um syringe filters | PVDF membrane | For filtering stock solutions |

### 3.4 Procedure

**A. Solvent Preparation**

1. Prepare the solvent: 1% v/v acetic acid + 0.2 M NaCl in ultrapure water.
   - Dissolve 11.69 g NaCl in 900 mL ultrapure water.
   - Add 10.0 mL glacial acetic acid.
   - Make up to 1000 mL. Mix thoroughly.

**B. Chitosan Solution Preparation (Day 1)**

2. Prepare a stock solution of 1.0% w/v chitosan:
   - Weigh 1.000 g chitosan (dry basis; correct for moisture content if known).
   - Add to 100 mL of the solvent prepared in step 1.
   - Stir at room temperature for 6-12 hours (or overnight) until fully dissolved.
   - Filter through a 0.45 um PVDF syringe filter to remove any undissolved material. Record the mass retained on the filter.
3. Prepare dilutions from the stock:
   - **Solution A (solvent blank):** pure solvent from step 1.
   - **Solution B (0.1% w/v):** dilute 10.0 mL stock to 100 mL with solvent.
   - **Solution C (0.25% w/v):** dilute 25.0 mL stock to 100 mL with solvent.
   - **Solution D (0.5% w/v):** dilute 50.0 mL stock to 100 mL with solvent.
   - **Solution E (1.0% w/v):** use stock directly.

**C. Viscometry Measurements (Day 2)**

4. Install the concentric cylinder geometry on the rheometer. Load the solvent trap.
5. Set the temperature to 90.0 degC. Allow 20 minutes for thermal equilibration.
6. Load Solution A (solvent blank) into the geometry (sufficient to cover the bob, typically 8-15 mL depending on geometry).
7. Apply an evaporation cover (solvent trap with low-volatility liquid seal or silicone oil seal).
8. Perform a steady-state flow sweep:
   - Shear rate range: 10 to 1000 s^-1 (logarithmic, 10 points per decade).
   - At each shear rate, wait for steady state (torque variation < 2% over 10 s).
   - Record the viscosity at each shear rate.
9. Identify the Newtonian plateau region. For dilute chitosan solutions at 90 degC, expect Newtonian behaviour across most of this range. Record the zero-shear viscosity eta\_solvent from the plateau.
10. Clean the geometry thoroughly with solvent, then ultrapure water, then dry.
11. Repeat steps 6-10 for Solutions B through E, in order of increasing concentration.
12. After Solution E, repeat Solution A as a control check.
13. Perform each measurement in **triplicate** (fresh loading each time).

### 3.5 Measurements

| Parameter | Value to Record | Replicates | Data Format |
|-----------|----------------|------------|-------------|
| Viscosity vs shear rate curve | Pa s at each shear rate | 3 per concentration | CSV: `c_g_mL, shear_rate_1_s, viscosity_Pa_s, T_K` |
| Zero-shear viscosity (eta\_0) | Pa s, from Newtonian plateau | 3 per concentration | Summary CSV: `c_g_mL, eta_0_Pa_s, eta_0_sd, T_K` |
| Specific viscosity (eta\_sp) | Dimensionless, calculated | -- | Derived: eta\_sp = eta\_0 / eta\_solvent - 1 |
| Reduced viscosity (eta\_sp/c) | mL/g | -- | Derived |

### 3.6 Data Analysis

**Fitting procedure:**

1. For each concentration c (in g/mL), compute the specific viscosity:

       eta_sp = eta_solution / eta_solvent - 1

2. Compute the reduced viscosity: eta\_red = eta\_sp / c (units: mL/g if c is in g/mL and we multiply by 1000 to convert).
3. Plot eta\_sp/c vs c. For the Huggins equation, this should be linear:

       eta_sp / c = [eta] + k_H * [eta]^2 * c

4. Perform a linear least-squares fit of eta\_sp/c vs c:
   - **Intercept** = [eta] (intrinsic viscosity, mL/g)
   - **Slope** = k\_H * [eta]^2
   - Therefore k\_H = slope / [eta]^2

5. Additionally, plot ln(eta\_sp/c) vs c for the Kraemer equation as a consistency check. The two extrapolations to c = 0 should converge to the same [eta].

6. **Software:** Python with `numpy.polyfit` (degree 1) or `scipy.optimize.curve_fit`.
7. **Initial guesses:** [eta] ~ 400-800 mL/g (chitosan at high temperature will be lower than 25 degC value); k\_H ~ 0.3-0.5.
8. Report [eta] with 95% confidence interval from the regression.

**Sample Python code:**

```python
import numpy as np

# c_data in g/mL, eta_sp_over_c in mL/g
coeffs = np.polyfit(c_data, eta_sp_over_c, 1)
slope, intercept = coeffs
eta_intrinsic = intercept  # mL/g
k_H = slope / eta_intrinsic**2
```

### 3.7 Error Treatment

| Source of Error | Mitigation | Acceptance Criterion |
|----------------|------------|---------------------|
| Evaporation at 90 degC | Use solvent trap; keep measurement time < 10 min per sample | Mass loss < 1% (weigh cup before/after) |
| Incomplete dissolution | Filter stock; record filtrate mass; re-dissolve if > 5% retained | < 2% mass on filter |
| Shear-thinning at high c | Use zero-shear plateau, not single-point viscosity | Plateau region spans at least 0.5 decade of shear rate |
| Thermal lag | Allow 5 min equilibration after loading before measurement | First data point discarded if > 5% above subsequent |
| Polyelectrolyte expansion | Use 0.2 M NaCl to screen charges | Reduced viscosity plot should be linear (R^2 > 0.95) |

**When to repeat:** If the Huggins plot is non-linear (curvature visible, R^2 < 0.95), the solutions may be above the dilute regime. Prepare additional dilutions at 0.05% and 0.02% w/v and repeat. If the solvent blank control at the end differs from the initial by > 3%, suspect evaporation artefacts and repeat the entire series.

### 3.8 Safety

- **Acetic acid (1% solution):** Mild irritant. Standard PPE. Low hazard at this dilution.
- **Chitosan powder:** Dust irritant. Weigh in fume hood or with dust mask.
- **Hot rheometer geometry (90 degC):** Burn hazard. Do not touch metal parts without heat-resistant gloves. Allow geometry to cool before cleaning.
- **NaCl:** Non-hazardous.

### 3.9 Expected Results

- **eta\_solvent (1% AcOH + 0.2 M NaCl at 90 degC):** approximately 0.3-0.4 mPa s (pure water at 90 degC is 0.315 mPa s; salt increases slightly).
- **[eta]\_chit at 90 degC:** Expected 300-600 mL/g. At elevated temperature, chain flexibility increases and excluded volume decreases relative to 25 degC, so [eta] should be 30-60% lower than room-temperature values (~800 mL/g).
- **k\_H:** Expected 0.3-0.7 for a good solvent / theta-solvent transition.
- **eta\_sp at 1.8% w/v (emulsification concentration):** Will be very large (>> 1), confirming that the Martin equation rather than Huggins should be used at the actual process concentration. The intrinsic viscosity from this study feeds into the Martin equation for the full simulation.

**Red flags:**
- [eta] < 100 mL/g: chitosan may be severely degraded. Check MW by GPC.
- [eta] > 1500 mL/g: insufficient charge screening. Increase NaCl concentration to 0.3 M.
- k\_H > 1.0: abnormal aggregation or association; check solution clarity.
- Non-monotonic eta\_sp/c vs c: likely measurement error at the lowest concentration (instrument near torque limit).

---

## 4. Study 3 -- Viscous Breakage Constant

### 4.1 Purpose

Calibrate the **Alopaeus viscous breakage constant C3** in the breakage rate kernel used by the Level 1 PBE solver (`src/emulsim/level1_emulsification/kernels.py`, function `breakage_rate_alopaeus`). C3 controls how strongly the dispersed-phase viscosity resists droplet breakup:

    g(d) = C1 * sqrt(epsilon/nu_c) * exp(-C2*sigma/(rho_c*epsilon^(2/3)*d^(5/3)) - C3*mu_d/sqrt(rho_c*sigma*d))

The current default is C3 = 0.0 (viscous correction disabled). For the highly viscous agarose/chitosan dispersed phase (mu\_d >> mu\_c), C3 is expected to be in the range 0.1-0.3 and has a dominant effect on the predicted droplet size.

### 4.2 Materials

| Material | Grade / Specification | CAS Number | Supplier | Quantity |
|----------|----------------------|------------|----------|----------|
| Light liquid paraffin | Same lot as Study 1 | 8012-95-1 | Sigma-Aldrich | 2 L |
| Span-80 | Same lot as Study 1 | 1338-43-8 | Sigma-Aldrich | 100 g |
| Agarose | Same lot as Study 1 | 9012-36-6 | Sigma-Aldrich | 50 g |
| Chitosan | Same lot as Study 1 | 9012-76-4 | Sigma-Aldrich | 25 g |
| Glacial acetic acid | ACS reagent | 64-19-7 | Fisher Scientific | 50 mL |
| Ultrapure water | >= 18.2 MOhm cm | 7732-18-5 | In-house | 2 L |

### 4.3 Equipment

| Instrument | Recommended Model | Key Specification |
|------------|-------------------|-------------------|
| Rotor-stator homogeniser | IKA Ultra-Turrax T25, Silverson L5M-A, or IKA Magic LAB | Variable speed 5,000-25,000 RPM; known rotor/stator geometry for epsilon estimation |
| Stator generator head | Fine-screen generator (IKA S25N-25F or equivalent) | Known gap width and slot dimensions |
| Laser diffraction particle sizer | Malvern Mastersizer 3000 or Horiba LA-960V2 | Range 0.01-3500 um; wet dispersion unit |
| Water bath | Large enough for mixing vessel | Stability +/- 1 degC at 90 degC |
| Jacketed glass vessel | 250-500 mL, with lid | For emulsification at controlled temperature |
| Digital tachometer | Handheld, contact or optical | To verify RPM independently |
| Stopwatch | -- | Timing of emulsification |

### 4.4 Procedure

**A. Experimental Design**

The goal is to isolate the viscous breakage mechanism. To do this, coalescence is suppressed by using:
- Low dispersed-phase holdup (phi\_d = 1 vol%, so coalescence is negligible)
- High Span-80 concentration (4% w/v, well above CMC)
- Short emulsification time (10-15 seconds, so steady-state is approached but coalescence has minimal effect)

Three dispersed-phase viscosities are obtained by varying the agarose/chitosan ratio while keeping total polymer at 6% w/v:
- **Low viscosity (L):** 2% agarose + 0% chitosan (no chitosan)
- **Medium viscosity (M):** 4.2% agarose + 1.8% chitosan (standard formulation)
- **High viscosity (H):** 3% agarose + 3% chitosan (high chitosan)

Each viscosity is tested at two RPMs:
- **RPM-1:** 8,000 RPM
- **RPM-2:** 15,000 RPM

This gives 3 x 2 = **6 runs**.

**B. Dispersed-Phase Preparation (Day 1)**

1. Prepare the three aqueous dispersed phases as described in Study 1, steps 5-8, but with the following polymer concentrations:
   - **Phase L:** 2.0 g agarose in 100 mL water (dissolve at 95 degC). No chitosan. Add 1 mL acetic acid after dissolution for pH consistency.
   - **Phase M:** 4.2 g agarose + 1.8 g chitosan in 100 mL (1% acetic acid). Standard preparation per Study 1.
   - **Phase H:** 3.0 g agarose + 3.0 g chitosan in 100 mL (1% acetic acid). Note: dissolve chitosan first (6-12 h), then add agarose and heat to dissolve.
2. Measure the viscosity of each phase at 90 degC using the rheometer (Study 2 setup) at a shear rate of 1000 s^-1 (representative of the homogeniser gap). Record mu\_d for each phase.
3. Maintain all three phases at 90 degC in a water bath.

**C. Oil Phase Preparation (Day 1)**

4. Prepare 2 L of continuous phase: 4.0% w/v Span-80 in liquid paraffin.
   - Weigh 80.0 g Span-80 into a 2 L bottle.
   - Add 2000 mL paraffin oil.
   - Mix at 40 degC for 2 hours.

**D. Emulsification Runs (Day 2)**

5. For each run, set up the jacketed vessel with 200 mL of oil phase at 90 degC.
6. Pre-heat the rotor-stator head by running it briefly in the hot oil (5 seconds at low speed).
7. Set the homogeniser to the target RPM. Verify with tachometer.
8. Using a pre-heated syringe (or pipette with cut tip), add 2.0 mL of the designated aqueous phase to the oil (phi\_d ~ 1 vol%). Start timing simultaneously.
9. Emulsify for exactly 15 seconds.
10. **Quench immediately:** Remove the homogeniser and plunge the vessel into an ice-water bath. The agarose will gel within 30-60 seconds as the temperature drops below ~38 degC, freezing the droplet size distribution.
11. Remove a 5 mL aliquot of the quenched emulsion and dilute into 50 mL of fresh paraffin oil + 2% Span-80 (room temperature) for particle sizing. Gentle inversion to mix; do not sonicate.
12. Repeat steps 5-11 for all 6 conditions (L-8k, L-15k, M-8k, M-15k, H-8k, H-15k).
13. Perform each run in **triplicate** (18 total runs).

**E. Particle Sizing (Day 2-3)**

14. Measure the droplet/particle size distribution of each sample using the laser diffraction particle sizer.
15. Use the wet dispersion unit with paraffin oil as the carrier liquid.
16. Set optical model parameters: refractive index of dispersed phase = 1.34 (hydrogel, close to water), absorption index = 0.001; refractive index of continuous phase = 1.48 (paraffin oil).
17. Add the diluted sample dropwise until the laser obscuration reaches 5-15%.
18. Run 5 consecutive measurements per sample with gentle stirring. Record the d32 (Sauter mean diameter), d10, d50, d90, and the full volume-weighted size distribution.
19. Compute the mean d32 and standard deviation across the 5 measurements and 3 replicates.

### 4.5 Measurements

| Parameter | Value to Record | Replicates | Data Format |
|-----------|----------------|------------|-------------|
| d32 (Sauter mean) | um | 3 runs x 5 readings = 15 per condition | CSV: `run_id, viscosity_label, RPM, mu_d_Pa_s, d32_um, d32_sd, d10_um, d50_um, d90_um, span` |
| Full PSD | Volume % vs diameter | 1 representative per condition | CSV: `diameter_um, volume_percent` |
| mu\_d at 90 degC, 1000 s^-1 | Pa s | 3 per phase | Separate CSV |
| Actual RPM | rev/min | Each run | Logged |
| sigma at 4% Span-80, 90 degC | N/m | From Study 1 | Cross-referenced |

### 4.6 Data Analysis

**Fitting procedure:**

1. For each of the 6 conditions, compute the mean experimental d32 (convert to metres).
2. For each condition, compute the simulation inputs:
   - epsilon (energy dissipation rate): from mixer geometry and RPM, using the `energy.py` functions in the codebase.
   - sigma: from Study 1 calibrated isotherm at c\_span80 = 4% w/v, T = 363.15 K.
   - rho\_c: measured or from `properties.toml` (850 kg/m^3).
   - mu\_d: measured at 90 degC for each phase.
   - rho\_aq: approximately 1020 kg/m^3.
3. Run the Level 1 PBE solver for each condition with C3 as a free parameter (C1 = 0.986, C2 = 0.0115 held fixed).
4. Define the objective function:

       J(C3) = sum_i [ (d32_sim_i(C3) - d32_exp_i)^2 / d32_exp_i^2 ]

   summed over all 6 conditions (i = 1..6).

5. Minimise J(C3) using `scipy.optimize.minimize_scalar` with bounds [0.0, 1.0].
6. **Initial guess:** C3 = 0.15.
7. Report the best-fit C3 with confidence interval estimated from the curvature of J(C3) around the minimum (approximate 1-sigma from Delta-J = 1/N).
8. Validate: run the PBE solver with the fitted C3 for all 6 conditions and plot d32\_sim vs d32\_exp. Points should fall on the 1:1 line with < 20% relative error.

**Sample Python code:**

```python
from scipy.optimize import minimize_scalar
from emulsim.level1_emulsification.solver import solve_pbe

def objective(C3):
    total_error = 0.0
    for cond in conditions:
        result = solve_pbe(rpm=cond['rpm'], sigma=cond['sigma'],
                           mu_d=cond['mu_d'], rho_c=rho_c,
                           breakage_C3=C3, ...)
        total_error += ((result.d32 - cond['d32_exp']) / cond['d32_exp'])**2
    return total_error

res = minimize_scalar(objective, bounds=(0.0, 1.0), method='bounded')
C3_fit = res.x
```

### 4.7 Error Treatment

| Source of Error | Mitigation | Acceptance Criterion |
|----------------|------------|---------------------|
| Coalescence during emulsification | Low holdup (1%), high surfactant, short time | If d32 increases with holdup (test at 0.5% vs 1%), coalescence is significant -- reduce holdup further |
| Incomplete agarose gelation during quench | Use ice bath; verify gelation by gently pressing particles (they should be solid) | Particles should not re-disperse upon re-homogenisation at room temperature |
| Particle sizing artefacts | Verify with optical microscopy on a subset; check for aggregation (bimodal peaks) | Laser diffraction and microscopy d32 should agree within 30% |
| RPM uncertainty | Verify with independent tachometer | Within +/- 5% of target |
| Temperature variation | Use jacketed vessel with circulating bath; log temperature | Within +/- 2 degC of 90 degC during the 15 s emulsification |
| epsilon estimation error | This is the largest systematic uncertainty; compare power-draw method and gap-based estimation | Use the same epsilon model as the simulation (consistency principle) |

**When to repeat:** If the fitted C3 yields residuals > 30% for any single condition, or if the overall RMSE > 20%, investigate that condition independently (repeat experiment, check mu\_d measurement, verify RPM).

### 4.8 Safety

- **Rotor-stator homogeniser:** Entanglement hazard. Never insert fingers near the running head. Ensure the generator is fully seated before operation. Wear safety glasses (splashing at high RPM).
- **Hot oil at 90 degC:** Severe burn hazard. Wear heat-resistant gloves, face shield, and lab coat. Keep water away from hot oil to prevent splattering.
- **Ice-water quench:** Use care when plunging hot vessel into ice bath; thermal shock can crack glass. Use borosilicate or stainless steel vessel.
- **Acetic acid in dispersed phase:** Low hazard at 1% but will be aerosolised during high-speed mixing. Use fume hood.
- **Paraffin oil mist:** Use local exhaust ventilation near the homogeniser.

### 4.9 Expected Results

- **d32 at 8,000 RPM, Phase M:** 5-30 um (depending on C3).
- **d32 at 15,000 RPM, Phase M:** 2-10 um.
- **d32 should decrease with increasing RPM** (approximately as RPM^(-1.2) from Kolmogorov scaling).
- **d32 should increase with increasing mu\_d** (higher chitosan content = larger drops).
- **C3 expected range:** 0.05-0.30. If C3 > 0.5, the viscous resistance is unrealistically high; check mu\_d measurements. If C3 < 0.01, the viscosity correction is negligible and the standard Coulaloglou-Tavlarides model (C3 = 0) may suffice.

**Red flags:**
- d32 independent of RPM: homogeniser may not be generating expected shear (check generator head, verify no cavitation at high RPM).
- d32 independent of mu\_d: viscosity range may be too narrow; add a very high viscosity condition (e.g., 5% chitosan).
- Bimodal PSD: could indicate two breakup mechanisms or droplet flocculation. Verify with microscopy.

---

## 5. Study 4 -- Pore Structure Empirical Coefficients

### 5.1 Purpose

Calibrate the empirical power-law relationship for mean pore diameter in agarose-only bulk gels as a function of polymer concentration and cooling rate:

    d_pore = A * c^alpha * (dT/dt)^beta

where A, alpha, and beta are empirical constants. The current simulation uses default coefficients (A = 600, alpha = -0.7, beta = -0.2) derived from the literature, but these must be validated for the specific agarose grade and cooling protocol used. The pore structure determines the chromatographic separation performance (Ogston partition coefficient K\_av) computed in Level 4 of the simulation.

### 5.2 Materials

| Material | Grade / Specification | CAS Number | Supplier | Quantity |
|----------|----------------------|------------|----------|----------|
| Agarose | Same lot as Studies 1 and 3 | 9012-36-6 | Sigma-Aldrich | 20 g |
| Ultrapure water | >= 18.2 MOhm cm | 7732-18-5 | In-house | 500 mL |
| Liquid nitrogen | Grade 4.6 | 7727-37-9 | BOC / Airgas | ~5 L (for cryo-SEM) |
| Glutaraldehyde | 25% aqueous, EM grade | 111-30-8 | Electron Microscopy Sciences (16220) | 50 mL |
| Ethanol | 200 proof, anhydrous | 64-17-5 | Sigma-Aldrich | 500 mL |
| tert-Butanol | >= 99.5%, ACS grade | 75-65-0 | Sigma-Aldrich (360538) | 200 mL |

**Note on pore visualisation approach:** Cryo-SEM is the primary method. If cryo-SEM is unavailable, an alternative is critical-point drying followed by conventional SEM, but this introduces shrinkage artefacts (typically 10-30% underestimation of pore size). The glutaraldehyde, ethanol, and tert-butanol are needed only for the alternative CPD route.

### 5.3 Equipment

| Instrument | Recommended Model | Key Specification |
|------------|-------------------|-------------------|
| Programmable cooling block / thermal cycler | TA Instruments DSC 2500 (for small samples), or custom Peltier block with PID controller | Controllable cooling rate from 0.5 to 50 degC/min; sample capacity >= 0.5 mL |
| Cryo-SEM preparation system | Quorum PP3010T or Leica EM VCT500 | Fracture stage, sputter coater, cryo-transfer shuttle |
| Field-emission SEM | JEOL JSM-7800F, FEI Quanta 250 FEG, or Zeiss Sigma 300 | Resolution < 5 nm; cryo-stage capability |
| SEM image analysis software | ImageJ/FIJI (free) with BoneJ or DiameterJ plugin | Pore segmentation and chord length measurement |
| Thermocouple data logger | Omega TC-08 or equivalent | K-type thermocouple, sampling >= 1 Hz |
| Silicone moulds | Cylindrical, 5 mm diameter x 3 mm height | For gel casting |

### 5.4 Procedure

**A. Experimental Design**

| Gel ID | Agarose Concentration (% w/v) | Cooling Rate (degC/min) |
|--------|-------------------------------|------------------------|
| G1 | 2.0 | 1 |
| G2 | 2.0 | 10 |
| G3 | 2.0 | 50 |
| G4 | 4.0 | 1 |
| G5 | 4.0 | 10 |
| G6 | 4.0 | 50 |
| G7 | 6.0 | 1 |
| G8 | 6.0 | 10 |
| G9 | 6.0 | 50 |

Each gel is prepared in triplicate = **27 gels total**.

**B. Agarose Solution Preparation (Day 1)**

1. Prepare three agarose stock solutions:
   - **Stock A:** 2.0 g agarose in 100 mL ultrapure water.
   - **Stock B:** 4.0 g agarose in 100 mL ultrapure water.
   - **Stock C:** 6.0 g agarose in 100 mL ultrapure water.
2. Dissolve each by heating to 95-100 degC with stirring until optically clear.
3. Hold at 90 degC in a water bath.

**C. Controlled Cooling and Gelation (Day 1)**

4. Transfer approximately 0.2 mL of the hot agarose solution into a silicone mould. Immediately insert a fine thermocouple (0.1 mm wire) into the centre of the gel.
5. Place the mould on the programmable cooling block.
6. Program the cooling profile:
   - **Start temperature:** 90 degC
   - **Cooling rate:** As per the experimental design table above (1, 10, or 50 degC/min)
   - **End temperature:** 4 degC
   - **Hold at 4 degC for 10 minutes** to ensure complete gelation.
7. Log the thermocouple temperature at >= 1 Hz throughout the cooling process.
8. For the 50 degC/min cooling rate: if the programmable block cannot achieve this rate, immerse the mould in a pre-cooled (4 degC) water-ethylene glycol bath and record the actual cooling rate from the thermocouple data. Use the measured average cooling rate through the gelation window (45-30 degC) for the fitting.
9. Repeat for all 9 conditions x 3 replicates = 27 gels.
10. Store gelled samples in sealed containers at 4 degC, immersed in a small volume of ultrapure water to prevent drying. Use within 48 hours.

**D. Cryo-SEM Imaging (Day 2-3)**

11. Mount a gel sample on the cryo-SEM stub using a thin layer of OCT compound or carbon adhesive.
12. Plunge-freeze the stub into liquid nitrogen slush (sub-cooled LN2 at ~-210 degC) to vitrify the water in the gel pores.
13. Transfer under vacuum to the cryo-preparation chamber.
14. Fracture the sample using the cold knife to expose a fresh cross-section.
15. Sublimate (etch) the surface at -95 degC for 3-5 minutes to reveal the pore structure by removing a thin layer of vitrified ice from the pore surfaces.
16. Sputter-coat with Pt/Pd (5-10 nm thickness) at cryo temperature.
17. Transfer to the SEM cryo-stage at -140 degC.
18. Image at multiple magnifications:
    - **Low magnification (500x-1000x):** Overview of pore network.
    - **Medium magnification (5,000x-10,000x):** For pore size measurement (aim for >= 50 pores per image).
    - **High magnification (20,000x-50,000x):** To resolve fine fibrillar structure (if relevant).
19. Capture at least **5 non-overlapping fields of view** per sample at the medium magnification.
20. Repeat for all 27 samples (or at minimum 1 replicate per condition = 9 samples, using the other replicates only if initial results have high variance).

**E. Image Analysis**

21. Open images in ImageJ/FIJI.
22. Apply a median filter (radius 2 pixels) to reduce noise.
23. Threshold to binarise the image into pore phase (dark) and fibre/wall phase (bright). Use Otsu's method for automatic thresholding. Verify visually and adjust manually if needed.
24. Measure pore sizes using one of:
    - **Chord length distribution:** Draw random lines across the binarised image and measure intercept lengths in the pore phase. Use the DiameterJ plugin.
    - **Equivalent circle diameter:** Segment individual pores using watershed; measure the area of each and convert to equivalent diameter d\_pore = sqrt(4*Area/pi).
25. For each sample, compute the **mean pore diameter** and **standard deviation** from >= 200 individual pore measurements (combined from 5 fields of view).

### 5.5 Measurements

| Parameter | Value to Record | Replicates | Data Format |
|-----------|----------------|------------|-------------|
| Mean pore diameter | nm | 3 gels per condition | CSV: `gel_id, c_agarose_pct, cooling_rate_C_min, d_pore_nm, d_pore_sd_nm, n_pores_measured` |
| Pore size distribution | nm | 1 per condition | Histogram CSV: `d_pore_nm, count` |
| Actual cooling rate | degC/min (average through 45-30 degC) | Each gel | From thermocouple log |
| Porosity | % (from image thresholding) | Each image | Derived from binarised image |

### 5.6 Data Analysis

**Fitting procedure:**

1. Tabulate the 9 condition-averaged mean pore diameters (in nm) with their concentrations (% w/v, equivalent to g/100mL) and cooling rates (degC/min).
2. Take the natural logarithm of the power-law model:

       ln(d_pore) = ln(A) + alpha * ln(c) + beta * ln(dT/dt)

3. Perform a multiple linear regression of ln(d\_pore) vs [ln(c), ln(dT/dt)] to obtain ln(A), alpha, and beta.

4. **Software:** Python `sklearn.linear_model.LinearRegression` or `numpy.linalg.lstsq`.
5. **Initial expectations:** A ~ 600 nm (at c = 1% w/v, dT/dt = 1 degC/min), alpha ~ -0.7, beta ~ -0.2.
6. Report the fitted A (in nm), alpha, and beta with standard errors from the regression.
7. Compute R^2. Expect R^2 > 0.90 for a well-controlled experiment.
8. Plot measured vs predicted d\_pore on a log-log scale with the 1:1 line.

**Sample Python code:**

```python
import numpy as np
from sklearn.linear_model import LinearRegression

# X: columns are [ln(c), ln(dT/dt)]
# y: ln(d_pore)
X = np.column_stack([np.log(c_data), np.log(cooling_rate_data)])
y = np.log(d_pore_data)

reg = LinearRegression().fit(X, y)
A = np.exp(reg.intercept_)
alpha = reg.coef_[0]
beta = reg.coef_[1]
```

### 5.7 Error Treatment

| Source of Error | Mitigation | Acceptance Criterion |
|----------------|------------|---------------------|
| Ice crystal formation (non-vitreous freezing) | Use LN2 slush (sub-cooled); keep samples thin (< 3 mm); plunge rapidly | No large (> 500 nm) ice crystal voids visible |
| Sublimation artefact (over-etching) | Optimise etch time (3-5 min at -95 degC); compare etched vs non-etched | Pore walls should be distinct but not eroded |
| Gel shrinkage during preparation | Cryo-SEM avoids this; if using CPD, apply a 1.2x correction factor | Cryo-SEM preferred |
| Non-uniform cooling | Use thin gel discs (3 mm); verify with thermocouple that centre temperature tracks the programmed profile within +/- 2 degC | dT/dt through gelation window within 20% of target |
| Subjectivity in thresholding | Use automated Otsu method; report threshold value; have a second analyst independently threshold 3 images | Inter-analyst d\_pore agreement within 15% |
| Sampling bias in pore measurement | Measure >= 200 pores per sample across >= 5 fields of view | Coefficient of variation of mean < 20% across fields |

**When to repeat:** If the power-law fit gives R^2 < 0.85, examine the residuals for systematic deviations. Common causes: (a) actual cooling rate differs from programmed rate (use actual thermocouple rate in the fit), (b) agarose gelation is too fast at high concentrations for slow cooling rates (the gel sets before the target rate is achieved). If a single condition is an outlier, repeat it. If the model systematically fails at high concentrations, consider a more complex model (e.g., including a saturation term).

### 5.8 Safety

- **Liquid nitrogen:** Cryogenic hazard. Wear cryo-gloves, face shield, and lab coat. Work in ventilated area (oxygen displacement risk in confined spaces). Do not seal LN2 in closed containers.
- **Glutaraldehyde (if using CPD route):** Toxic (GHS06), sensitiser. Handle in fume hood only. Wear double nitrile gloves. Dispose as hazardous waste.
- **tert-Butanol:** Flammable (GHS02), irritant. Keep away from heat. Fume hood.
- **Ethanol:** Flammable. Standard precautions.
- **Hot agarose solutions (95-100 degC):** Burn hazard. Use heat-resistant gloves.
- **SEM:** Electron beam. Follow institutional SEM safety training. No specific chemical hazards during imaging.

### 5.9 Expected Results

- **d\_pore at 2% agarose, 1 degC/min:** 200-500 nm (large pores, slow cooling allows extensive fibre aggregation).
- **d\_pore at 6% agarose, 50 degC/min:** 30-80 nm (small pores, dense network, rapid nucleation).
- **d\_pore at 4% agarose, 10 degC/min:** 80-200 nm (the reference condition closest to the simulation default).
- **alpha:** Expected -0.5 to -0.9 (negative; higher concentration gives smaller pores).
- **beta:** Expected -0.1 to -0.4 (negative; faster cooling gives smaller pores).
- **A:** Expected 300-1000 nm (depending on unit convention; here c in % w/v and dT/dt in degC/min).

**Red flags:**
- Pore size increasing with agarose concentration: sample preparation error or imaging artefact.
- Pore size insensitive to cooling rate: cooling rate may not have been achieved (check thermocouple data).
- Bimodal pore size distribution: could indicate phase separation during gelation (unusual for pure agarose/water but possible if contaminants are present).
- Mean pore size < 10 nm: likely a thresholding error or resolution artefact.

---

## 6. Study 5 -- IPN Coupling Coefficient

### 6.1 Purpose

Calibrate the **IPN coupling coefficient eta\_coupling** in the double-network modulus model used in `src/emulsim/level4_mechanical/solver.py` (function `double_network_modulus`):

    G_DN = G_agarose + G_chitosan + eta_coupling * sqrt(G_agarose * G_chitosan)

eta\_coupling captures the synergistic (positive) or antagonistic (negative) interaction between the two interpenetrating networks. The current default is eta\_coupling = -0.15 (mildly antagonistic, from mutual swelling constraint without sacrificial bonds).

The calibration strategy is to measure the shear modulus of three matched gels -- agarose-only, chitosan+genipin-only, and the double-network (DN) gel -- and solve for eta\_coupling algebraically.

### 6.2 Materials

| Material | Grade / Specification | CAS Number | Supplier | Quantity |
|----------|----------------------|------------|----------|----------|
| Agarose | Same lot as previous studies | 9012-36-6 | Sigma-Aldrich | 10 g |
| Chitosan | Same lot as previous studies | 9012-76-4 | Sigma-Aldrich | 10 g |
| Genipin | >= 98% HPLC, from Gardenia jasminoides | 6902-77-8 | Challenge Bioproducts (Taiwan) or Wako (078-03021) | 1 g |
| Glacial acetic acid | ACS reagent | 64-19-7 | Fisher Scientific | 20 mL |
| Sodium hydroxide | ACS reagent, pellets | 1310-73-2 | Sigma-Aldrich | 20 g |
| Phosphate-buffered saline (PBS) | pH 7.4 tablets or 10x concentrate | -- | Sigma-Aldrich (P4417) | 1 L |
| Ultrapure water | >= 18.2 MOhm cm | 7732-18-5 | In-house | 1 L |

### 6.3 Equipment

| Instrument | Recommended Model | Key Specification |
|------------|-------------------|-------------------|
| Oscillatory rheometer | TA Instruments DHR-3, Anton Paar MCR 302 | Normal force control; oscillatory mode |
| Parallel plate geometry | 20 mm diameter, sandblasted (to prevent wall slip) | Gap precision +/- 1 um |
| Peltier temperature control | Plate + hood | Range 4-100 degC, +/- 0.1 degC |
| Gel moulds | Cylindrical, 20 mm diameter x 2 mm height (matching plate geometry) | Silicone or PTFE |
| Incubator | 37 degC, humidified | For genipin crosslinking |
| pH meter | Mettler Toledo FiveEasy or equivalent | Accuracy +/- 0.01 pH units |

### 6.4 Procedure

**A. Experimental Design**

Two formulation sets are tested to provide a consistency check:

| Set | Agarose (% w/v) | Chitosan (% w/v) | Genipin (mM) | Cooling Rate (degC/min) |
|-----|-----------------|------------------|--------------|------------------------|
| **Set 1** (standard) | 4.2 | 1.8 | 2.0 | 10 |
| **Set 2** (high chitosan) | 3.0 | 3.0 | 3.0 | 10 |

For each set, three gel types are prepared:

- **Type A (agarose-only):** Agarose at the specified concentration, no chitosan, no genipin.
- **Type C (chitosan+genipin-only):** Chitosan at the specified concentration, crosslinked with genipin. No agarose. (This gel will be weak; handle carefully.)
- **Type DN (double-network):** Both agarose and chitosan, crosslinked with genipin.

Each gel type in triplicate = 2 sets x 3 types x 3 replicates = **18 gels**, but the triplet structure (A + C + DN at matched conditions) means 6 matched triplets.

**B. Gel Preparation (Day 1)**

*Type A gels (agarose-only):*

1. Dissolve agarose in ultrapure water at 95-100 degC at the target concentration.
2. Pour 0.8 mL into 20 mm silicone moulds (gives ~2 mm height).
3. Cool at 10 degC/min using the programmable cooling block (same protocol as Study 4).
4. Hold at 4 degC for 10 minutes.
5. Transfer gels (still in moulds) to PBS at room temperature. Equilibrate for 2 hours.

*Type C gels (chitosan+genipin-only):*

6. Dissolve chitosan in 1% acetic acid at the target concentration.
7. Adjust pH to 6.0 using 1 M NaOH (dropwise, with stirring and pH monitoring). This partially deprotonates chitosan and initiates physical gelation. Work quickly -- the solution will become viscous.
8. Add genipin from a 100 mM stock solution (in DMSO or warm water) to achieve the target concentration. Mix thoroughly for 30 seconds.
9. Pour 0.8 mL into 20 mm moulds.
10. Place moulds in the 37 degC incubator for 24 hours to allow genipin crosslinking. The gel will turn deep blue (genipin-chitosan reaction chromophore).
11. Transfer to PBS at room temperature. Equilibrate for 2 hours.

*Type DN gels (double-network):*

12. Dissolve chitosan in 1% acetic acid.
13. Dissolve agarose separately in water at 95-100 degC.
14. Combine the two solutions at > 85 degC (as in Study 1, step 8).
15. Add genipin while the solution is still hot (> 80 degC). Mix for 30 seconds. Note: at 80-90 degC, genipin reaction is slow, so there is a working window of ~5 minutes before significant crosslinking occurs.
16. Pour 0.8 mL into moulds.
17. Cool at 10 degC/min to set the agarose network (first network).
18. Transfer to 37 degC incubator for 24 hours for genipin crosslinking of chitosan (second network).
19. Transfer to PBS at room temperature. Equilibrate for 2 hours.

**C. Oscillatory Rheology (Day 3)**

20. Set up the rheometer with the 20 mm sandblasted parallel plates.
21. Set temperature to 25.0 degC. Allow 15 minutes for thermal equilibration.
22. Carefully demould a gel disc. Measure its diameter and height with calipers.
23. Place the gel on the lower plate. Lower the upper plate until the normal force reads 0.05-0.1 N (light contact, < 5% compression). Record the gap.
24. Apply a thin film of light mineral oil around the exposed edge to prevent drying during measurement.
25. **Amplitude sweep (strain sweep):**
    - Frequency: 1 Hz (6.28 rad/s)
    - Strain range: 0.01% to 10% (logarithmic, 10 points per decade)
    - Purpose: identify the linear viscoelastic region (LVER)
26. **Frequency sweep:**
    - Strain: within the LVER (typically 0.1-1% for these gels)
    - Frequency range: 0.1 to 100 rad/s (logarithmic, 5 points per decade)
    - Record G' (storage modulus) and G'' (loss modulus) at each frequency.
27. Report G' at 1 Hz (6.28 rad/s) as the representative shear modulus for fitting.
28. Repeat for all 18 gels.

### 6.5 Measurements

| Parameter | Value to Record | Replicates | Data Format |
|-----------|----------------|------------|-------------|
| G' at 1 Hz | Pa | 3 per gel type per set | CSV: `set_id, gel_type, replicate, G_prime_Pa, G_prime_sd, G_double_prime_Pa, strain_pct, T_K` |
| Full frequency sweep | G', G'' vs omega | 1 per gel type per set | CSV: `omega_rad_s, G_prime_Pa, G_double_prime_Pa` |
| Amplitude sweep | G', G'' vs strain | 1 per gel type per set | CSV: `strain_pct, G_prime_Pa, G_double_prime_Pa` |
| Gel dimensions | mm | Each gel | Appended to main CSV |
| Gel colour | Visual description | Each gel | Notes (blue intensity for chitosan gels indicates crosslinking extent) |

### 6.6 Data Analysis

**Fitting procedure:**

1. For each set, extract the mean G' values at 1 Hz for the three gel types:
   - G\_agar = G' of Type A gel (agarose-only)
   - G\_chit = G' of Type C gel (chitosan+genipin-only)
   - G\_DN = G' of Type DN gel (double-network)

2. Solve for eta\_coupling algebraically from the model equation:

       eta_coupling = (G_DN - G_agar - G_chit) / sqrt(G_agar * G_chit)

3. Compute eta\_coupling independently for each of the 2 formulation sets. If the two values agree within experimental uncertainty, report the weighted mean. If they differ significantly, the coupling model may need refinement (e.g., concentration-dependent coupling).

4. Propagate uncertainty using standard error propagation:

       delta_eta = eta * sqrt( (delta_G_DN/G_DN)^2 + ... )

   (full expression involves partial derivatives with respect to G\_agar, G\_chit, and G\_DN).

5. **Software:** Python or spreadsheet calculation.
6. Report eta\_coupling with 95% confidence interval.

**Sample Python code:**

```python
import numpy as np

# For each set:
G_agar = np.mean(G_agar_replicates)
G_chit = np.mean(G_chit_replicates)
G_DN = np.mean(G_DN_replicates)

eta_coupling = (G_DN - G_agar - G_chit) / np.sqrt(G_agar * G_chit)
```

### 6.7 Error Treatment

| Source of Error | Mitigation | Acceptance Criterion |
|----------------|------------|---------------------|
| Gel slip on plates | Use sandblasted plates; apply light normal force; verify no slip by checking G' independence of gap | G' should not change by > 5% when gap is varied +/- 10% |
| Incomplete genipin crosslinking | Allow full 24 hours at 37 degC; verify by deep blue colour of gel | Gel should be uniformly dark blue throughout |
| Gel drying during measurement | Apply mineral oil seal; keep measurement time < 15 min per gel | G'' / G' ratio should be < 0.15 (mostly elastic) |
| Inhomogeneous gel structure | Pour solutions quickly; ensure uniform cooling; visual inspection | Gel should be optically uniform (no visible gradients) |
| Normal force relaxation | Allow 60 s after plate contact before starting measurement | Normal force stable within +/- 10% |
| Temperature variation between gel types | Run all three types of a matched triplet on the same day | Temperature recorded within +/- 0.5 degC |

**When to repeat:** If eta\_coupling from the two sets differ by more than a factor of 2, or if the standard deviation of eta\_coupling within a set (from propagated errors across triplicates) exceeds 50% of the mean value, prepare and measure fresh gels. Check that Type C gels have achieved sufficient crosslinking (G\_chit should be > 100 Pa for 1.8% chitosan with genipin).

### 6.8 Safety

- **Genipin:** Moderate toxicity. Causes skin staining (blue). Wear double nitrile gloves. Avoid inhalation of powder (weigh in fume hood). Not classified as acutely toxic but limited toxicological data -- treat with caution.
- **Sodium hydroxide (1 M):** Corrosive (GHS05). Wear gloves and safety glasses.
- **DMSO (if used for genipin stock):** Rapidly penetrates skin and carries dissolved solutes through. Wear chemical-resistant gloves (not standard nitrile -- use butyl rubber if handling concentrated DMSO/genipin).
- **Hot agarose solutions:** Burn hazard. Standard precautions.
- **Mineral oil for rheometer edge seal:** Low hazard. Avoid ingestion.

### 6.9 Expected Results

- **G\_agar (4.2% w/v agarose):** 30-80 kPa (literature range for standard agarose at this concentration).
- **G\_chit (1.8% chitosan + 2 mM genipin, 24 h):** 0.5-5 kPa (chitosan gels are generally much weaker than agarose gels).
- **G\_DN:** If additive (eta = 0), G\_DN ~ G\_agar + G\_chit. If antagonistic (eta < 0), G\_DN < G\_agar + G\_chit. If synergistic (eta > 0), G\_DN > G\_agar + G\_chit.
- **eta\_coupling expected:** -0.3 to +0.2. The current default (-0.15) predicts mild antagonism. True DN gels with sacrificial bond mechanisms can show eta > 0 and even very large positive values, but a sequential IPN without intentional sacrificial bonds is expected to be weakly antagonistic or neutral.

**Red flags:**
- G\_DN < G\_agar (the DN gel is weaker than agarose alone): severe antagonism. Possible cause: chitosan solution at acidic pH partially degrades the agarose network, or genipin interferes with agarose helix aggregation. This would give eta\_coupling << -1.
- G\_chit ~ 0: genipin crosslinking failed. Check genipin stock freshness (should be stored at -20 degC, desiccated). Gels should be blue.
- G\_DN >> G\_agar + G\_chit (eta > 1): strong synergy. Unusual for this system but would be a scientifically interesting finding. Verify by repeating.
- G' frequency-dependent (slope > 0.1 on log-log): gel may not be fully formed. Extend curing time.

---

## 7. Inputting Calibrated Values into the Simulation

After completing all five studies and performing the data analysis, update the simulation with the calibrated values as described below.

### 7.1 TOML Property Database

Update `data/properties.toml` with the fitted values:

```toml
[interfacial]
# From Study 1
sigma_bare = { value = <measured_sigma_0>, unit = "N/m", source = "Study 1, pendant drop at 90 degC" }
gamma_inf = { value = <fitted_Gamma_inf>, unit = "mol/m2", source = "Study 1, Szyszkowski-Langmuir fit" }
K_L = { value = <fitted_K_L>, unit = "m3/mol", source = "Study 1, Szyszkowski-Langmuir fit" }
```

### 7.2 Source Code Updates

**Interfacial tension model** (`src/emulsim/properties/interfacial.py`):

Update the hardcoded values in `interfacial_tension_span80()` at lines 49-50:

```python
Gamma_inf = <fitted_value>   # mol/m2, from Study 1
K_L = <fitted_value>         # m3/mol, from Study 1
```

Also update `sigma_bare()` if the measured bare IFT differs from the current linear model by more than 10%.

**Chitosan intrinsic viscosity** (`src/emulsim/properties/viscosity.py`):

Update the hardcoded value at line 169:

```python
eta_intr_chit = <fitted_value>  # mL/g, from Study 2 (at 90 degC)
```

Also update `k_H` if it differs from the default of 0.4 (used in the `huggins_viscosity` function).

### 7.3 MaterialProperties Dataclass

Update `src/emulsim/datatypes.py`, class `MaterialProperties`:

```python
# From Study 3
breakage_C3: float = <fitted_C3>    # [-] Alopaeus viscous correction, Study 3

# From Study 5
eta_coupling: float = <fitted_eta>  # IPN coupling coefficient, Study 5
```

### 7.4 Pore Structure Model

The pore structure coefficients (A, alpha, beta from Study 4) are currently embedded in the Level 2 gelation analysis. After calibration, update the relevant pore diameter prediction function. If no explicit function exists in the current codebase, add one to `src/emulsim/level2_gelation/pore_analysis.py`:

```python
def empirical_pore_diameter(c_agarose_pct: float, cooling_rate_C_min: float,
                             A: float = <fitted_A>,
                             alpha: float = <fitted_alpha>,
                             beta: float = <fitted_beta>) -> float:
    """Empirical mean pore diameter [nm] from Study 4 calibration.

    d_pore = A * c^alpha * (dT/dt)^beta
    """
    return A * c_agarose_pct**alpha * cooling_rate_C_min**beta
```

And add the coefficients to `data/properties.toml`:

```toml
[pore_structure]
A_pore = { value = <fitted_A>, unit = "nm", source = "Study 4, power-law fit" }
alpha_pore = { value = <fitted_alpha>, unit = "-", source = "Study 4" }
beta_pore = { value = <fitted_beta>, unit = "-", source = "Study 4" }
```

### 7.5 Default Configuration

Update `configs/default.toml` if any process parameters (e.g., RPM, Span-80 concentration) are revised based on the calibration results. The formulation block:

```toml
[formulation]
c_span80 = 20.0   # [kg/m3] -- verify sigma at this concentration matches Study 1
```

### 7.6 Validation Checklist

After updating all values, run the full simulation pipeline and verify:

1. [ ] The predicted sigma at 2% Span-80, 90 degC matches Study 1 within 10%.
2. [ ] The predicted mu\_d for the standard formulation at 90 degC is consistent with the value used in Study 3.
3. [ ] The predicted d32 at the standard conditions (10,000 RPM, 2% Span-80, phi\_d = 5%) is in the expected range (2-20 um).
4. [ ] The predicted pore size at 4.2% agarose, 10 degC/min cooling is consistent with Study 4.
5. [ ] The predicted G\_DN for the standard formulation is consistent with Study 5.
6. [ ] All unit tests pass (`pytest tests/`).

### 7.7 Uncertainty Propagation

For each calibrated parameter, record the 95% confidence interval. These uncertainties should be entered into the `uncertainty` field of the `PropertyValue` entries in `properties.toml` and used in the uncertainty analysis module (`src/emulsim/uncertainty.py`) for Monte Carlo propagation through the simulation.

---

## Appendix A -- Reagent Master List

| # | Reagent | CAS | Grade | Quantity Needed | Storage |
|---|---------|-----|-------|----------------|---------|
| 1 | Light liquid paraffin | 8012-95-1 | Ph. Eur. / USP | 2.5 L | RT, sealed |
| 2 | Span-80 (sorbitan monooleate) | 1338-43-8 | >= 60% oleic acid | 150 g | RT, sealed, dark |
| 3 | Agarose (low-EEO, standard melting) | 9012-36-6 | Molecular biology grade | 100 g | RT, dry |
| 4 | Chitosan (high DD) | 9012-76-4 | MW 100-300 kDa, DD >= 90% | 50 g | RT, dry, desiccated |
| 5 | Glacial acetic acid | 64-19-7 | ACS reagent, >= 99.7% | 200 mL | RT, fume hood |
| 6 | Genipin | 6902-77-8 | >= 98% HPLC | 1 g | -20 degC, desiccated, dark |
| 7 | Sodium chloride | 7647-14-5 | ACS reagent | 20 g | RT |
| 8 | Sodium hydroxide (pellets) | 1310-73-2 | ACS reagent | 20 g | RT, sealed (hygroscopic) |
| 9 | PBS tablets | -- | Cell culture grade | 10 tablets | RT |
| 10 | Ethanol (200 proof) | 64-17-5 | Anhydrous | 500 mL | RT, flammables cabinet |
| 11 | tert-Butanol | 75-65-0 | >= 99.5% | 200 mL | RT, flammables cabinet |
| 12 | Glutaraldehyde (25%) | 111-30-8 | EM grade | 50 mL | 4 degC, sealed |
| 13 | Liquid nitrogen | 7727-37-9 | Grade 4.6 | ~5 L | Dewar |

Items 10-13 are only needed if using the CPD alternative for Study 4 or for cryo-SEM.

## Appendix B -- Equipment Summary

| Instrument | Studies Used In | Estimated Instrument Time |
|------------|----------------|--------------------------|
| Pendant-drop tensiometer (heated cell) | 1 | 8-12 hours |
| Rotational rheometer (concentric cylinder) | 2 | 4-6 hours |
| Rotational rheometer (parallel plates, oscillatory) | 5 | 6-8 hours |
| Rotor-stator homogeniser | 3 | 4-6 hours |
| Laser diffraction particle sizer | 3 | 4-6 hours |
| Programmable cooling block | 4, 5 | 8-12 hours |
| Cryo-SEM (with cryo-preparation system) | 4 | 12-20 hours |
| Incubator (37 degC) | 5 | Passive (24 h per batch) |
| Thermocouple data logger | 4 | Passive during cooling |
| Analytical balance (0.01 mg) | All | Shared |
| pH meter | 5 | 1 hour |

---

*End of Protocol*

*Cross-reference: This calibration protocol supplies the material constants described in Sections 1A.2 (interfacial tension), 1A.1 (viscous breakup), 1A.4 (PBE kernels), 1B.1 (agarose gelation/pore structure), and 1B.3 (IPN architecture) of [01 -- Scientific Advisor Report](01_scientific_advisor_report.md).*
