# EmulSim Modules 2 & 3: Integrated System Design

**Date:** 2026-04-11
**Synthesised from:** Scientific Advisor, Dev-Orchestrator, Architect
**Scope:** Complete architecture for expanding EmulSim from 1-module to 3-module system

---

## Executive Summary

Three specialist roles independently designed the next-phase architecture. Key consensus:

- **Module 2 heavily reuses existing L3 infrastructure** — same ODE systems, Arrhenius kinetics,
  Thiele modulus. New code is mainly ACS tracking and ligand coupling.
- **Module 3 requires fundamentally new solvers** — 1D advection-dispersion PDE (method of lines),
  adsorption isotherms (Langmuir, SMA), detection models (Beer-Lambert, ESI-MS).
- **Single orchestrator extended** with opt-in Module 2/3 (backward compatible).
- **ACSProfile is the central bookkeeping object** bridging Module 1 → Module 2.
- **FunctionalMicrosphere bridges Module 2 → Module 3**.

### Scale

| Metric | Module 2 | Module 3 | Integration | Total |
|--------|----------|----------|-------------|-------|
| New sub-modules | 8 | 16 | 4 | 28 |
| Estimated LOC | ~1,450 | ~2,500 | ~700 | ~4,650 |
| Opus-tier modules | 1 | 3 | 0 | 4 |
| New dataclasses | 8 | 10 | 0 | 18 |
| New reagent profiles | ~8 | ~12 | 0 | ~20 |
| New trust checks | 5 | 5 | 0 | 10 |

---

## Architecture Overview

```
Module 1 (existing)         Module 2 (NEW)              Module 3 (NEW)
┌──────────────────┐    ┌───────────────────┐    ┌──────────────────────┐
│ L1: Emulsification│    │ ACS Initialization │    │ Column/Bed Setup     │
│ L2: Gelation      │───>│ Modification Loop  │───>│ Transport PDE        │
│ L3: Crosslinking  │    │  (N sequential     │    │ Isotherm Models      │
│ L4: Mechanical    │    │   steps)           │    │ Detection Simulation │
└──────────────────┘    └───────────────────┘    └──────────────────────┘
     FullResult         FunctionalMicrosphere    ChromatogramResult /
                                                  CatalyticResult
```

### Data Flow at Module Boundaries

**Module 1 → Module 2:**
- d50/2 (bead radius), pore_size_mean, porosity → surface area, ACS accessibility
- p_final (crosslink conversion) → residual NH2 = NH2_0 × (1 - p_final)
- xi_final (mesh size) → diffusion constraint for modification reagents
- G_DN, E_star → mechanical baseline

**Module 2 → Module 3:**
- Ligand type + density → isotherm selection + q_max
- Updated pore size (after crosslinking/coupling) → D_eff for mass transfer
- G_DN (updated) → max operating pressure
- Enzyme loading × f_activity → V_max,obs for catalytic simulation

---

## 5 Real-World Use Cases

| # | Application | Key Input | Key Output |
|---|------------|-----------|------------|
| 1 | **IMAC** (His-tag purification) | NTA density, Ni²⁺ loading, imidazole gradient | Breakthrough curve, elution peak, DBC₁₀% |
| 2 | **Ion Exchange** (protein separation) | DEAE/SP density, NaCl gradient, protein pI values | Multi-peak chromatogram, resolution Rs |
| 3 | **Protein A** (mAb capture) | Protein A density, pH elution step | Full chromatogram, yield, HCP clearance |
| 4 | **Enzyme reactor** (lipase catalysis) | CALB loading, substrate conc, flow rate | Conversion, effectiveness factor, activity decay |
| 5 | **ACS modification workflow** | Multi-step: ECH→EDA→NHS coupling | ACS waterfall chart, conversion per step |

---

## Key Scientific Models

### Module 2: Chemical Modification

| Model | Equation | Application |
|-------|----------|-------------|
| ACS consumption | d[ACS]/dt = -k·[reagent]·[ACS] | All modification steps |
| Sequential ACS | [NH2]_remaining = [NH2]₀·(1-p₁) | After primary crosslinking |
| DS model | DS(t) = DS_max·(1-exp(-k·[reagent]·t)) | ACS conversion (OH→COOH) |
| Ligand coupling | d[L_coupled]/dt = k·[active_sites]·[L_free] | Ligand immobilization |
| Activity retention | f = f_coupling · f_orientation · f_conformational | Enzyme immobilization |

### Module 3: Performance Simulation

| Model | Equation | Application |
|-------|----------|-------------|
| Lumped Rate Model | ε∂C/∂t = -u·∂C/∂z + D_ax·∂²C/∂z² - (1-ε)·k_f·(C-C_p) | Column transport |
| Langmuir isotherm | q* = q_max·K_L·C/(1+K_L·C) | Affinity/IMAC binding |
| SMA (IEX) | q* = f(K_eq, z, σ, C_salt) | Ion exchange |
| Michaelis-Menten | v = V_max·[S]/(K_m+[S]) | Enzyme catalysis |
| Effectiveness factor | η = (1/Φ)·(1/tanh(3Φ) - 1/(3Φ)) | Diffusion-limited catalysis |
| Beer-Lambert | A = ε·c·l | UV detection |
| ESI charge states | m/z = (M + z·1.008)/z | MS detection |

---

## File Structure (New)

```
src/emulsim/
  module2_functionalization/
    __init__.py
    solver.py              # ACS init, crosslink, conversion, coupling solvers
    orchestrator.py        # Sequential step execution
  module3_performance/
    __init__.py
    solver.py              # Top-level chromatography/catalysis dispatch
    transport.py           # PDE discretization (method of lines + BDF)
    isotherms.py           # Langmuir, SMA, Freundlich, competitive
    detection.py           # UV, fluorescence, conductivity, MS simulation
  reagent_library_m2.py    # 8 activation/conversion reagents (CNBr, ECH, DVS, NHS...)
  reagent_library_m3.py    # 12 functional ligands (Protein A, DEAE, NTA, lipase...)
```

---

## Build Order (5 Phases)

```
Phase A: Data Structures + Libraries (Haiku, ~5 modules)
  datatypes extensions, reagent_library_m2, reagent_library_m3,
  protein_library, ligand_library

Phase B: Module 2 Core (Sonnet, ~5 modules)
  ACS tracker, activation solver, spacer solver,
  ligand coupling solver, modification orchestrator

Phase C: Module 3 Core (Opus-heavy, ~8 modules)
  Column model, gradient generator, isotherm models [Opus],
  transport solver [Opus], catalytic solver [Opus],
  detection models (UV/fluorescence/MS/conductivity),
  chromatogram generator, performance orchestrator

Phase D: Visualization (Sonnet, ~3 modules)
  Module 2 plots (ACS waterfall, conversion charts)
  Module 3 plots (chromatograms, breakthrough, Michaelis-Menten)

Phase E: Integration (Sonnet, ~4 modules)
  Extended pipeline orchestrator, UI tabs for M2/M3
```

---

## Trust Gate Extensions (10 new checks)

**Module 2:** ACS deficit (blocker), pore accessibility, pH compatibility,
sequential consistency, unreacted sites warning

**Module 3:** Mechanical pressure limit, pore exclusion, under-resolved column,
mass balance, Thiele modulus (catalysis)

---

## Detection Output Formats

| Detector | Model | Plot Type |
|----------|-------|-----------|
| UV 280nm | Beer-Lambert: A = ε·c·l | mAU vs time/CV with gradient overlay |
| Fluorescence | F = Φ·I₀·ε·c·l | RFU vs time |
| ESI-MS | Charge envelope at m/z = (M+z·1.008)/z | TIC + mass spectrum panel |
| Conductivity | κ = Σ(λᵢ·cᵢ) | mS/cm vs time (gradient monitoring) |
| pH | Henderson-Hasselbalch | pH vs time (elution monitoring) |

---

> Full detailed reports from each specialist role are available in the task outputs.
> This document synthesizes the consensus architecture for implementation.
