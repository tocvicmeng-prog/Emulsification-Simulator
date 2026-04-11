# Audit Report: Module 2 & Module 3 Integrated System Design

Date: 2026-04-11

Reviewed document: [`docs/11_module2_module3_system_design.md`](11_module2_module3_system_design.md)

Reviewer perspective: computational simulation science, polymer/bioconjugation chemistry, chromatographic purification science, immobilized enzyme reactor engineering

## Executive Verdict

The proposed three-module expansion is directionally correct: Module 1 generates base microspheres, Module 2 modifies or functionalizes them, and Module 3 predicts use-performance in purification or catalytic process environments.

However, the current design document is still a **high-level architecture sketch**, not yet an implementation-ready scientific plan. It correctly identifies several important building blocks, such as ACS bookkeeping, ligand coupling, column transport, adsorption isotherms, detection models, and catalytic performance. But it does not yet define enough state variables, chemistry-specific stoichiometry, surface-area models, transport limitations, validation strategy, numerical safeguards, or trust gates to support reliable simulation of real functionalized chromatographic or catalytic microspheres.

The largest weakness is that Module 2 is described as “mainly ACS tracking and ligand coupling,” but the requirement actually implies a much richer physicochemical system:

- reactive-site inventory
- accessibility of internal vs external sites
- surface-normalized ACS density
- sequential consumption, conversion, activation, deactivation, and ligand immobilization
- tradeoff between reinforcement crosslinking and later functional ligand capacity
- chemistry-specific pH, buffer, ionic strength, hydrolysis, steric, diffusion, and orientation effects

The largest weakness in Module 3 is that the proposed transport and isotherm layer is under-specified for real chromatography and catalytic beds. A simple 1D advection-dispersion model plus Langmuir/SMA labels is not sufficient unless it is expanded into a mass-transfer, pressure-drop, multi-component, detector-aware column model with rigorous boundary conditions and calibration.

## Requirement Alignment

### Correctly captured

The plan correctly captures these core requirements:

- Module 1 remains the existing microsphere-preparation simulator.
- Module 2 should simulate chemical modification, reinforcement, crosslinking, activation, and ligand coupling of base microspheres.
- Module 2 should track ACS as a central concept.
- Module 2 output should feed Module 3 as a `FunctionalMicrosphere`.
- Module 3 should simulate functional microsphere behavior in columns or packed beds.
- Module 3 should support purification ligands and catalytic ligands.
- Module 3 should generate process outputs such as breakthrough curves, chromatograms, peak shapes, conversion, and catalytic product formation.

### Partially captured

The plan partially captures the requirement for:

- ACS density, but does not sufficiently define total internal plus external surface area.
- Crosslinking density, but does not distinguish theoretical ligand capacity, actual coupled ligand density, active ligand density, and inaccessible ligand density.
- Functional-property tradeoffs, but does not define how secondary crosslinking consumes ACS and reduces later immobilization capacity.
- Real operating conditions, but does not yet include adequate pressure-drop, bed-compression, flow maldistribution, detector dispersion, or temperature/pH/salt coupling.

### Missing or under-specified

The plan does not yet adequately specify:

- ACS spatial distribution within the bead.
- ACS accessibility as a function of pore size, ligand/reagent size, charge, and swelling.
- Internal surface-area estimation from pore morphology.
- Distinction between total ACS, accessible ACS, activated ACS, consumed ACS, blocked ACS, and ligand-functional ACS.
- Chemistry-specific reaction networks for CNBr, ECH, DVS, NHS, aldehyde, epoxide, maleimide, IMAC chelation, IEX ligand coupling, Protein A coupling, and enzyme immobilization.
- Multi-step process histories and carryover of reagent residues.
- Ligand orientation and retained activity, especially for proteins and enzymes.
- Column hydrodynamics, pressure drop, bed porosity, particle porosity, resin compressibility, and packing heterogeneity.
- Multi-component competitive adsorption and realistic elution gradients.
- Detector delay volume, band broadening, response factors, noise, saturation, and matrix effects.

## Major Findings

### 1. High severity: ACS density is central, but the surface-area model is missing

**Issue**

The requirement defines ACS density as:

`number of available ACS / total internal and external microsphere surface area`

The plan states that bead radius, pore size, and porosity feed surface area and ACS accessibility, but it does not specify how surface area is computed.

**Why this matters**

ACS density is the key bridge between Module 1 and Module 2. If the surface-area denominator is wrong, all downstream outputs become wrong:

- maximum ligand loading
- crosslinking density
- functional group consumption
- adsorption capacity
- enzyme loading
- catalytic productivity

The internal surface area of a hydrogel microsphere cannot be inferred from radius alone. It depends on pore topology, polymer fiber thickness, tortuosity, swelling state, and whether pores are actually accessible to the reagent or ligand.

**Required solution**

Define an `AccessibleSurfaceModel` with at least three levels:

- `external_only`: `A_ext = 4*pi*R^2`
- `empirical_pore`: estimate internal area from porosity and pore diameter, for example `S_v ~ 4*epsilon_pore/d_pore` for simplified cylindrical pore geometry, with calibration flags
- `morphology_based`: compute accessible area from a resolved phase field or reconstructed pore network

The model must output:

- external area
- internal geometric area
- reagent-accessible internal area
- ligand-accessible internal area
- uncertainty or trust level

### 2. High severity: Module 2 cannot be treated as a light extension of existing L3

**Issue**

The plan says Module 2 heavily reuses existing L3 infrastructure and that new code is mainly ACS tracking and ligand coupling.

That is too optimistic.

Existing L3 mostly describes one class of crosslinking kinetics. Module 2 requires a broader chemical process simulator:

- activation of low-reactivity groups
- conversion of one ACS type into another ACS type
- spacer-arm insertion
- ligand immobilization
- protein/enzyme coupling with orientation/activity loss
- hydrolysis and quenching
- residual reactive group accounting
- sequential reaction history

**Why this matters**

Many real functionalization chemistries are not simple second-order irreversible reactions. Examples:

- NHS esters hydrolyze in water and compete with amine coupling.
- CNBr activation is pH-sensitive and hazardous, with activated intermediates that decay.
- Epoxide and DVS coupling require alkaline conditions and have slow/fast side reactions.
- Protein ligands are sterically hindered and orientation-dependent.
- Enzyme immobilization depends on multipoint attachment, activity retention, and deactivation.

**Required solution**

Module 2 should be designed as a sequential reaction-network simulator, not only as an ACS counter.

Recommended core object:

```text
ACSState
  total_sites[type]
  accessible_sites[type]
  activated_sites[type]
  consumed_sites[type]
  blocked_sites[type]
  ligand_sites[type]
  spatial_profile[type, r] optional
```

Each modification step should update `ACSState` through chemistry-specific stoichiometry and accessibility rules.

### 3. High severity: ACS density and crosslinking density must be separated more rigorously

**Issue**

The plan correctly notes that ACS density and crosslinking density are correlated. But it risks treating them as too directly interchangeable.

**Why this matters**

High ACS density does not automatically mean high useful ligand density. Actual ligand density depends on:

- accessible ACS fraction
- activation efficiency
- ligand size
- steric shielding
- reaction time
- ligand concentration
- diffusion into pores
- reagent hydrolysis/deactivation
- local charge and pH
- whether previous crosslinking consumed the same ACS

**Required solution**

Module 2 should report separate quantities:

- `ACS_total_density`
- `ACS_accessible_density`
- `ACS_activated_density`
- `crosslinker_consumed_density`
- `ligand_coupled_density`
- `active_ligand_density`
- `remaining_ACS_density`

For protein/enzyme ligands, `active_ligand_density` must be lower than `ligand_coupled_density` due to orientation and conformational retention.

### 4. High severity: Module 3 transport model is under-specified for chromatography

**Issue**

The plan lists a lumped-rate transport equation:

`epsilon*dC/dt = -u*dC/dz + D_ax*d2C/dz2 - (1-epsilon)*k_f*(C-C_p)`

This is a good start, but insufficient as a column model.

**Why this matters**

Real resin performance depends on coupled effects:

- axial dispersion
- film mass transfer
- intraparticle diffusion
- adsorption/desorption kinetics
- bed porosity
- particle porosity
- interstitial velocity
- pressure drop
- compressibility
- non-specific binding
- gradient delay and mixing
- detector broadening

A single concentration field and one film term will not capture breakthrough or peak shape reliably across affinity, IEX, IMAC, and catalytic modes.

**Required solution**

Use a tiered chromatography solver:

- `equilibrium_dispersive_model`: fast, assumes instantaneous binding
- `lumped_rate_model`: includes film/LDF particle mass transfer
- `general_rate_model`: includes radial intraparticle diffusion for high-fidelity cases

At minimum, the implementation should define:

- mobile-phase concentration `C_i(z,t)`
- particle-phase pore concentration `C_p,i(z,t)`
- bound concentration `q_i(z,t)`
- axial dispersion `D_ax`
- film mass transfer `k_f`
- intraparticle effective diffusivity `D_eff`
- binding kinetics or equilibrium isotherm
- Danckwerts inlet/outlet boundary conditions

### 5. High severity: Column pressure and resin mechanics must be first-class Module 3 outputs

**Issue**

The plan mentions backpressure and mechanical pressure limit trust gates, but does not specify a pressure-drop model.

**Why this matters**

For chromatographic resins, mechanical strength is not only a material property. It determines:

- maximum flow rate
- bed compression
- permeability loss
- channeling risk
- usable dynamic binding capacity
- column lifetime

**Required solution**

Module 3 should include:

- Ergun or Kozeny-Carman pressure drop model
- permeability as a function of bead diameter, bed porosity, and compression
- resin compressibility linked to Module 1/2 mechanical outputs
- maximum allowable pressure and flow-rate envelope

The model should compute:

- `deltaP(flow_rate)`
- `max_safe_flow_rate`
- `bed_compression_fraction`
- `mechanical_failure_risk`

### 6. High severity: Adsorption model selection is too simple

**Issue**

The plan lists Langmuir and SMA, but does not define when they apply or how multi-component competition is handled.

**Why this matters**

Different ligand classes require different adsorption physics:

- IMAC: metal-chelate coordination, imidazole competition, pH dependence, metal loading/leaching
- IEX: salt-dependent electrostatic binding, pH/pI/protein charge, SMA or steric mass-action models
- Protein A: pH-dependent affinity, mass-transfer limitations, possible ligand leaching
- Hydrophobic interaction: ammonium sulfate dependence, protein-specific hydrophobicity
- Protein/enzyme immobilization: activity and steric accessibility, not only loading

**Required solution**

Implement an `IsothermModel` interface with chemistry-specific parameters and validity ranges:

- `LangmuirIsotherm`
- `CompetitiveLangmuirIsotherm`
- `SMAIsotherm`
- `pHDependentAffinityIsotherm`
- `IMACCompetitionIsotherm`
- `NonSpecificBindingTerm`

All isotherm outputs should include a validity warning if salt, pH, ligand, or protein properties are outside calibrated ranges.

### 7. High severity: Catalytic-bed model needs reaction-transport coupling, not just Michaelis-Menten

**Issue**

The plan lists Michaelis-Menten and an effectiveness factor, but this is not sufficient for packed-bed catalytic simulation.

**Why this matters**

Immobilized enzyme beds couple:

- convective residence time
- intraparticle diffusion
- substrate depletion
- product formation
- inhibition
- temperature dependence
- enzyme deactivation
- local pH
- cofactor concentration if relevant

The effectiveness factor equation listed in the plan is also not robustly stated. The standard first-order spherical catalyst effectiveness factor is commonly written in forms equivalent to:

`eta = 3/phi * (1/tanh(phi) - 1/phi)`

For Michaelis-Menten kinetics there is generally no simple universal closed-form effectiveness factor without approximation.

**Required solution**

The catalytic module should support:

- plug-flow packed-bed reactor model
- substrate and product PDEs
- immobilized enzyme loading from Module 2
- `Vmax_obs = active_enzyme_loading * kcat`
- Arrhenius temperature correction
- deactivation model
- optional inhibition terms
- numerical intraparticle diffusion or validated effectiveness approximation

### 8. Medium severity: Detection models are oversimplified

**Issue**

The plan lists Beer-Lambert, fluorescence, ESI charge states, conductivity, and pH.

These are useful, but real detector outputs need instrument response modeling.

**Why this matters**

Column effluent peaks are shaped by both column physics and extra-column effects:

- mixer volume
- tubing volume
- detector flow-cell volume
- detector response time
- baseline drift
- noise
- saturation
- matrix effects
- ion suppression in MS

The proposed ESI-MS model is especially simplified. A charge-state equation gives possible `m/z` values, but not ionization efficiency, suppression, isotope envelope, adducts, fragmentation, or quantitative response.

**Required solution**

Separate detector simulation into:

- ideal concentration-to-signal conversion
- instrument response convolution
- noise/baseline model
- detector-specific limitations

MS should initially be framed as qualitative or semi-quantitative unless calibrated response factors are supplied.

### 9. Medium severity: Multi-component feed chemistry is not yet represented

**Issue**

The requirement explicitly names non-target solutes such as buffers, ammonium sulfate, CHAPS, imidazole, surfactants, and related components. The plan mentions gradients and some isotherms, but does not define a general multi-component fluid composition object.

**Why this matters**

These components are not passive labels. They affect:

- protein charge
- ionic strength
- hydrophobic interaction
- IMAC competition
- enzyme activity
- viscosity
- density
- detector baseline
- adsorption selectivity

**Required solution**

Add a `MobilePhaseComposition` model:

```text
MobilePhaseComposition
  pH
  conductivity
  buffer_species
  salt_species
  additives
  target_solutes
  impurity_solutes
  viscosity
  density
```

This object should drive isotherm parameters, detector baselines, and catalytic activity corrections.

### 10. Medium severity: The file structure is too compressed for maintainability

**Issue**

The plan proposes only a few files for very broad functionality:

- `solver.py`
- `transport.py`
- `isotherms.py`
- `detection.py`

**Why this matters**

Module 3 alone spans chromatography, catalysis, detection, hydrodynamics, pressure, packed-bed mechanics, gradients, and process outputs. If implemented in a few large files, it will become difficult to test and validate.

**Required solution**

Use more explicit subpackages:

```text
module2_functionalization/
  acs.py
  surface_area.py
  reactions.py
  modification_steps.py
  ligand_coupling.py
  protein_coupling.py
  orchestrator.py

module3_performance/
  column_geometry.py
  hydrodynamics.py
  pressure_drop.py
  transport/
    equilibrium_dispersive.py
    lumped_rate.py
    general_rate.py
  isotherms/
    langmuir.py
    sma.py
    imac.py
    protein_a.py
  catalysis/
    packed_bed.py
    kinetics.py
    deactivation.py
  detection/
    uv.py
    fluorescence.py
    conductivity.py
    ms.py
  orchestrator.py
```

## Module 2 Detailed Review

### ACS initialization

The plan should initialize ACS from:

- polymer composition from Module 1 formulation
- chitosan DDA for amines
- cellulose/agarose hydroxyl availability
- carboxymethyl substitution degree if CMC is present
- residual ACS after Module 1 primary crosslinking
- accessible fraction from pore/mesh model

The current statement `residual NH2 = NH2_0 * (1 - p_final)` is only valid if `p_final` is explicitly the fraction of amines consumed. It should be guarded by stoichiometry metadata from the primary crosslinking chemistry.

### ACS accessibility

Accessibility should depend on:

- reagent hydrodynamic radius
- ligand hydrodynamic radius
- pore-size distribution
- mesh size
- charge exclusion
- tortuosity
- swelling state

Recommended relation:

- small molecules can access most pore surfaces
- proteins access only pores larger than a threshold relative to hydrodynamic radius
- immobilized proteins can block pores and reduce subsequent accessibility

### Sequential modification

The plan correctly includes a modification loop, but it should explicitly support:

- step order
- reagent excess
- buffer/pH/temperature per step
- conversion efficiency
- side reactions
- hydrolysis/deactivation
- wash/removal efficiency
- quenching
- residual activated groups

### Ligand and enzyme coupling

Protein and enzyme ligands need special handling:

- orientation factor
- conformational retention
- multipoint attachment
- steric crowding
- activity loss
- leaching/stability
- pore exclusion

The proposed `f = f_coupling * f_orientation * f_conformational` is useful but should be expanded to include calibration ranges and ligand-specific defaults.

## Module 3 Detailed Review

### Chromatography outputs

The target outputs should include:

- breakthrough curve
- chromatogram
- elution peak area
- peak height
- retention volume
- peak width
- asymmetry/tailing factor
- resolution between components
- dynamic binding capacity at 5%, 10%, and 50% breakthrough
- yield
- purity
- productivity
- buffer consumption
- pressure drop

### Column and bed geometry

Required inputs:

- column inner diameter
- bed height
- bed volume
- particle size distribution
- bed porosity
- particle porosity
- extra-column volume
- tubing/detector volume
- operating temperature
- flow rate or linear velocity
- pressure limit

### Mass balance

Every Module 3 solver must pass mass-balance tests:

- total injected mass
- effluent mass
- bound mass
- reacted mass for catalytic mode
- remaining mass in column at final time

Mass balance should be a trust gate and a test criterion.

### Numerical methods

Recommended baseline:

- conservative finite-volume discretization in `z`
- implicit or semi-implicit method-of-lines integration
- positivity preservation for concentrations and bound states
- Danckwerts boundary conditions
- adaptive timestepping
- sparse Jacobian support for stiff binding/catalysis

The plan should explicitly define numerical tolerances and acceptance tests.

## Revised Data Model Recommendation

### Module 2 core dataclasses

Recommended minimum dataclasses:

```text
ACSProfile
  site_type
  total_sites_per_particle
  total_density_per_area
  accessible_density_per_area
  activated_density_per_area
  consumed_density_per_area
  remaining_density_per_area
  accessibility_model
  uncertainty

ModificationStep
  reagent
  target_site_type
  product_site_type
  stoichiometry
  pH
  temperature
  time
  reagent_concentration
  rate_model
  side_reactions

FunctionalMicrosphere
  base_result
  acs_profiles
  ligand_profiles
  updated_mechanics
  pore_accessibility
  process_history
```

### Module 3 core dataclasses

Recommended minimum dataclasses:

```text
ColumnGeometry
  diameter
  bed_height
  bed_volume
  particle_size_distribution
  bed_porosity
  extra_column_volume

MobilePhaseComposition
  pH
  ionic_strength
  salt_concentration
  buffer_species
  additives
  solutes

ChromatographyMethod
  equilibration
  load
  wash
  elution_gradient
  regeneration
  flow_rate
  temperature

PerformanceResult
  chromatograms
  breakthrough_curves
  pressure_trace
  mass_balance
  yield
  purity
  productivity
  trust_assessment
```

## Revised Build Order

The proposed build order starts with data structures and libraries, which is correct. But the implementation should be stricter.

### Phase 0: Resolve Module 1 dependencies

Before Module 2/3 outputs are trusted, Module 1 must provide stable and clearly labeled outputs:

- bead radius or bead size distribution
- pore size distribution
- porosity
- mesh size
- residual ACS after primary crosslinking
- mechanical strength
- trust level

If Module 1 result is empirical or low-trust, Module 2 and 3 must inherit that uncertainty.

### Phase A: Data model and ACS surface-area core

Implement:

- `ACSProfile`
- `AccessibleSurfaceModel`
- `FunctionalMicrosphere`
- surface-area tests
- ACS conservation tests

### Phase B: Minimal Module 2 chemistry

Start with only two validated workflows:

- amine coupling or genipin-like secondary crosslinking on chitosan
- epoxy/DVS-like hydroxyl activation on agarose/cellulose-like polysaccharide

Do not implement broad reagent libraries before the core accounting is validated.

### Phase C: Minimal Module 3 chromatography

Start with a single-component Langmuir breakthrough model:

- one solute
- one ligand
- no gradient
- finite-volume 1D transport
- mass balance
- pressure drop
- UV detection

Only after that should the system add IMAC, IEX, Protein A, multi-component mixtures, and gradients.

### Phase D: Catalytic packed bed

Implement catalysis after chromatography transport is validated.

Start with:

- one substrate
- one product
- immobilized enzyme loading from Module 2
- plug-flow with axial dispersion
- Michaelis-Menten kinetics
- optional first-order deactivation

### Phase E: Multi-component purification and realistic detection

Add:

- competitive adsorption
- salt/pH gradients
- detector convolution
- MS semi-quantitative output
- resolution/yield/purity metrics

## Trust Gate Recommendations

### Module 2 trust gates

Required checks:

- ACS conservation across sequential steps
- ACS accessibility not exceeding total ACS
- activated sites not exceeding accessible sites
- ligand-coupled sites not exceeding activated sites
- pH compatibility for each reagent
- pore accessibility for each reagent and ligand
- protein ligand size exclusion
- excessive pore blockage warning
- secondary crosslinking reduces later ligand capacity
- residual reactive groups after quenching

### Module 3 trust gates

Required checks:

- column pressure below mechanical limit
- solute hydrodynamic radius compatible with pore access
- axial grid resolution adequate for peak width
- mass balance closure within tolerance
- non-negative concentrations and loadings
- bed Reynolds number within model range
- Peclet number within dispersion correlation range
- binding capacity not exceeding ligand density
- detector saturation warning
- unsupported isotherm/chemistry warning

## Validation Plan

### Module 2 validation

Minimum wet-lab validation targets:

- total amine density by TNBS/ninhydrin assay
- hydroxyl activation or epoxide density by titration or probe coupling
- ligand loading by elemental analysis, UV assay, BCA/Bradford for proteins, or ICP for IMAC metals
- residual ACS after crosslinking
- mechanical strength before and after secondary crosslinking
- pore-size change after modification

### Module 3 validation

Minimum process validation targets:

- tracer residence-time distribution
- pressure drop vs flow rate
- breakthrough curve for a single known protein
- DBC10% vs residence time
- salt-gradient elution peak for IEX
- imidazole elution for IMAC
- pH elution for Protein A-like affinity
- catalytic conversion vs flow rate for enzyme beds

## Implementation Risk Assessment

| Risk | Severity | Reason |
|---|---:|---|
| ACS surface area poorly defined | High | Breaks all density and capacity calculations |
| Chemistry library too broad too early | High | Produces unsupported simulations with false confidence |
| Module 3 transport under-specified | High | Peak shapes and breakthrough become unreliable |
| Pressure-drop/mechanics omitted | High | Resin may be predicted usable at impossible flow rates |
| Detector models oversold | Medium | Simulated spectra/chromatograms may look more quantitative than they are |
| Validation delayed | High | Model complexity grows without falsifiability |
| Existing Module 1 uncertainty ignored | High | Downstream results inherit upstream errors silently |

## Recommended Reframing of the Plan

The design should be reframed from:

> “Build Modules 2 and 3 with ACS tracking, ligand coupling, column PDEs, isotherms, and detection.”

to:

> “Build a validation-gated, tiered process simulator in which Module 2 converts Module 1 microsphere structure into chemically accessible functional-site inventories, and Module 3 converts those inventories into pressure-constrained chromatographic or catalytic performance predictions with explicit mass balance, transport limits, and detector response assumptions.”

This framing is more scientifically accurate and better aligned with the original requirement.

## Final Verdict

The current Module 2/3 plan is a strong conceptual starting point, but it is not yet sufficiently rigorous for implementation as a reliable scientific simulator.

The most important corrections are:

- define ACS density through a real accessible surface-area model
- make Module 2 a sequential chemistry and accessibility simulator, not just a reuse of L3
- implement Module 3 as a tiered chromatography/catalysis transport framework with pressure drop and mass balance
- constrain chemistry libraries to validated workflows before broad expansion
- add validation and trust gates from the beginning

If implemented with these corrections, the three-module EmulSim system could become a credible platform for designing functional polysaccharide microsphere resins and catalytic media. Without these corrections, it risks becoming a visually impressive but scientifically overclaimed simulator.
