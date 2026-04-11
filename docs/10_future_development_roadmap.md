# EmulSim Future Development Roadmap

**Date:** 2026-04-11
**Author:** Scientific Advisor (Claude Opus 4.6)
**Status:** Strategic planning document

---

## 1. Current State Assessment

### What EmulSim Does Well

- **End-to-end L1-L4 pipeline** mirroring the real fabrication process
- **Chemistry breadth:** 8 crosslinkers across 4 mechanistic families, 6 surfactants
- **Trust infrastructure:** 3-category uncertainty propagation + mode-aware trust gates
- **Ternary CH solver** with semi-implicit Eyre splitting, no-flux BCs, gelation arrest
- **Bayesian optimization** with qLogEHVI and feasibility constraints
- **No comparable open-source competitor** covers this full pipeline

### Remaining Limitations

- L1 kernel constants not calibrated against experimental data
- L2 default is empirical correlation (mechanistic CH available but slower)
- No embedded validation datasets
- Single-droplet simulation (no batch variability)

---

## 2. Near-Term Enhancements (1-3 months)

### 2.1 Validation Dataset Architecture
Build a `data/validation/` directory with calibration datasets:
- d32 vs RPM (L1), pore vs agarose concentration (L2), G vs concentration (L4)
- p_final vs genipin concentration (L3), pore vs cooling rate (L2)
- Implement regression tests comparing model output to data

### 2.2 L1 RPM Recalibration
Refit C1/C2/C3 against published d32-vs-RPM data for agarose-in-oil emulsification.
Refs: Alopaeus et al. (2002) Chem. Eng. Sci. 57:1815; Zhou et al. (2005) Biotechnol. Bioeng. 91:609

### 2.3 Avrami alpha_final in Empirical L2
Replace hardcoded alpha_final=0.999 with Avrami kinetics (~15 LOC).

### 2.4 DDA as User-Adjustable Parameter
Expose degree of deacetylation (0.70-0.95) in UI — critically affects NH2 concentration.

---

## 3. Medium-Term Upgrades (3-6 months)

### 3.1 Accelerated Ternary Solver
- FFT-based spectral operators (O(N^2 log N) vs O(N^4) for direct solve)
- Adaptive time stepping
- Target: N=64 in <5s, enabling mechanistic L2 as practical default

### 3.2 Radial/Spherical Geometry for L2
Solve on axisymmetric domain to capture confinement-induced morphology for small beads.
Ref: Lee et al. (2014) Macromolecules 47:6667

### 3.3 Batch Variability Model
Run L2-L4 for each droplet size class from L1 DSD. Report batch-level property distributions.
Critical for chromatography applications requiring uniform beads.

### 3.4 ML Surrogate Models
Expose BoTorch GP surrogates for fast screening outside the optimization loop.

### 3.5 Thermal-Mass Transfer Coupling (L1↔L2)
Model incipient gelation during emulsification via Newton cooling + viscosity increase.
Ref: Gu et al. (2014) Chem. Eng. Sci. 117:387

---

## 4. Long-Term Research Directions (6-12 months)

### 4.1 Other Microsphere Platforms
- **Alginate:** Ca2+ diffusion gelation, closest to current framework
- **PLGA:** W/O/W double emulsion + solvent extraction (fundamentally different L2)
- **Cellulose:** NIPS instead of TIPS (ternary CH applicable)

### 4.2 Digital Twin Capability
Real-time process data assimilation via Ensemble Kalman filter.
Online Bayesian parameter estimation + predictive RPM control.

### 4.3 Inverse Design
Given target specs (d32, pore, G_DN, Kav), find optimal formulation via constrained BO.

### 4.4 Robust Optimization Under Uncertainty
Optimize E[f(x + delta)] instead of f(x) — find formulations insensitive to perturbations.

### 4.5 MD Parameter Estimation
Use MARTINI coarse-grained MD for chi_ac, kappa_CH, M_0, f_bridge.
Ref: Lopez et al. (2015) J. Chem. Theory Comput. 11:2714

---

## 5. Experimental Validation Priorities

### Minimal Validation Campaign (24 formulations, ~3 weeks)

| Study | Experiments | Measures | Purpose |
|-------|------------|----------|---------|
| 1: RPM sweep | 5 | d32, DSD | L1 calibration + monotonicity |
| 2: Agarose conc | 4 | pore, G_agar | L2 + L4 calibration |
| 3: Cooling rate | 3 | pore, porosity | L2 beta exponent |
| 4: Genipin conc | 4 | p_final, G_DN | L3 + L4 kinetics |
| 5: Time course | 4 | p_final, G_DN | L3 kinetics validation |
| 6: Surfactant conc | 4 | d32, sigma | L1 IFT model |

### Key Techniques
- Laser diffraction (d32), cryo-SEM (pore), oscillatory rheometry (G),
  ninhydrin assay (p_final), pendant drop tensiometry (sigma)

---

## 6. Publication Opportunities

| # | Title Concept | Target Journal | Effort |
|---|---------------|---------------|--------|
| 1 | EmulSim framework paper | Comput. Chem. Eng. / CES | 3-4 months |
| 2 | Ternary CH with gelation arrest | Phys. Rev. E / Soft Matter | 4-6 months |
| 3 | Multi-objective BO for hydrogel beads | AIChE J. / JCIM | 4-5 months |
| 4 | Crosslinker chemistry screening guide | J. Chromatogr. A | 5-6 months |

---

## Strategic Principle

> **Validate before extending, calibrate before predicting, be honest about what the model can and cannot do.**

EmulSim's greatest risk is not missing features — it is that existing features may be
trusted beyond their validated domain. The near-term priority should be making the
default path scientifically defensible, not adding more physics.

---

> **Disclaimer**: This roadmap is provided for research planning purposes only.
> All experimental designs should be reviewed by qualified domain experts.
