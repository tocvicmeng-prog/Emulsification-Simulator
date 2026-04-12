# EmulSim Module 2 — Detailed Integration Plan

## New Ligand, Protein, and Linker Arm Candidates

**Version:** 1.0
**Date:** 2026-04-12
**Prepared by:** Scientific Advisor + Architect + Dev-Orchestrator

---

## Table of Contents

1. [Goal Statement](#1-goal-statement)
2. [Realization Method](#2-realization-method)
3. [Specific Process — Ligand Coupling Candidates](#3-specific-process--ligand-coupling-candidates)
4. [Specific Process — Protein Coupling Candidates](#4-specific-process--protein-coupling-candidates)
5. [Specific Process — Linker Arm Candidates](#5-specific-process--linker-arm-candidates)
6. [Specific Process — FunctionalMediaContract Extension](#6-specific-process--functionalmediacontract-extension)
7. [Specific Process — UI Integration](#7-specific-process--ui-integration)
8. [Implementation Roadmap](#8-implementation-roadmap)
9. [Validation and Acceptance](#9-validation-and-acceptance)

---

## 1. Goal Statement

### 1.1 Overall Goal

Expand EmulSim Module 2's functionalization capability from a **6-ligand + 2-protein demonstration library** to a **production-grade library of 22 coupling profiles + 3 spacer profiles** covering the five major chromatography modes used in industrial biopharmaceutical downstream processing:

| Chromatography Mode | Current M2 Coverage | Target Coverage |
|---|---|---|
| Ion Exchange (IEX) | 2 of 4 (DEAE, SP) | **4 of 4** (+ Q, CM) |
| Hydrophobic Interaction (HIC) | 1 of 2 (Phenyl) | **2 of 2** (+ Butyl) |
| Immobilized Metal Affinity (IMAC) | 1 of 2 (IDA) | **2 of 2** (+ NTA) |
| Affinity (Protein) | 2 of 4 (Protein A, G) | **4 of 4** (+ A/G, Streptavidin) |
| Affinity (Small molecule) | 0 of 2 | **2 of 2** (+ Glutathione, Heparin) |

Additionally, introduce **spacer arm support** so the simulation can model the effect of linker arms on protein ligand activity retention — a critical physical parameter that determines the binding capacity of affinity chromatography media.

### 1.2 Why This Matters

EmulSim models the complete fabrication-to-performance pipeline for hydrogel microspheres. Module 2 (functionalization) is the bridge between microsphere fabrication (M1) and chromatographic performance simulation (M3). Without adequate ligand/protein coverage, M2 can only simulate a narrow subset of real-world chromatography media, limiting the tool's practical utility for:

- **Process development scientists** designing chromatography steps for mAb purification
- **Resin manufacturers** screening bead formulations for new media products
- **Academic researchers** modelling the relationship between microsphere properties and column performance

### 1.3 Specific Goals

| # | Goal | Measurable Outcome |
|---|---|---|
| G1 | Complete IEX ligand coverage | User can select all 4 standard IEX types (DEAE, Q, SP, CM) in M2 |
| G2 | Add NTA for His-tag IMAC | User can simulate Ni-NTA resin preparation |
| G3 | Add mild HIC option | Butyl available as lower-hydrophobicity alternative to Phenyl |
| G4 | Add GST-tag and heparin affinity | User can simulate glutathione and heparin affinity media |
| G5 | Add broad-spectrum antibody affinity | Protein A/G fusion and streptavidin coupling available |
| G6 | Model spacer arm effects | Protein coupling activity retention reflects spacer arm improvement |
| G7 | Map new ligands to M3 capacity | FunctionalMediaContract generates q_max estimates for all new types |

---

## 2. Realization Method

### 2.1 Architectural Approach: Minimal Extension, No Restructuring

The integration follows the principle of **minimal architectural change**. The ACS terminal-state model (v2), the 5 step types, and the ODE solver infrastructure remain unchanged. New candidates are added purely as:

1. **New `ReagentProfile` entries** in the existing `REAGENT_PROFILES` dictionary
2. **New fields** on the existing `ReagentProfile` dataclass (4 fields)
3. **Extended mapping logic** in `FunctionalMediaContract` builder
4. **Extended UI dropdown options** in the Streamlit app

**Scientific reason:** The coupling chemistry for all new candidates is identical to existing candidates — nucleophilic ring-opening of epoxide by amine (or thiol). The ODE template (`_competitive_hydrolysis_rhs`) already models this correctly. What changes between candidates is the kinetic parameters (k_forward, E_a, hydrolysis_rate), the physical parameters (MW, hydrodynamic radius), and the downstream functional mode (IEX, HIC, IMAC, affinity). These are all ReagentProfile attributes, not architectural components.

### 2.2 Spacer Arm Approach: Activity Multiplier (Phase 1)

**Method:** Model spacer arms as a **scalar multiplier** on `activity_retention`, stored as a ReagentProfile attribute.

```
effective_activity = activity_retention * spacer_activity_multiplier
ligand_functional_sites = ligand_coupled_sites * effective_activity
```

**Scientific reason:** The primary experimentally measurable effect of a spacer arm is an increase in the fraction of immobilized protein that retains biological activity. Published Protein A studies show:

- Direct coupling (no spacer): 50-70% activity retention
- DADPA spacer (13 A): 65-80% activity retention
- PEG spacer (35 A): 75-85% activity retention

The ratio (with-spacer / without-spacer) = 1.15-1.30x is the spacer_activity_multiplier. This is a first-order empirical correction that:

1. Captures the dominant physical effect (reduced steric hindrance → more accessible binding domains)
2. Is within the uncertainty band of the activity_retention parameter itself (+/- 0.15)
3. Does not require new ACS states, new step types, or new ODE solvers
4. Can be upgraded to a full SPACER_ARM step type in Phase 2 without breaking existing workflows

**What this approach does NOT model (acknowledged limitations):**

- Spacer-dependent coupling kinetics (e.g., AHA's -COOH terminus requires EDC/NHS activation, a different chemistry from DADPA's -NH2 terminus coupling directly to epoxide)
- Spacer hydrophobicity effects on non-specific binding
- Multi-step spacer-then-couple workflow sequencing

These are deferred to Phase 2 and documented in the confidence tier as "semi_quantitative".

### 2.3 Charge Type Approach: Explicit Field

**Method:** Replace the current string-matching heuristic for IEX anion/cation detection with an explicit `charge_type` field on ReagentProfile.

**Scientific reason:** The current code uses `"anion" in installed_ligand.lower() or "deae" in installed_ligand.lower()` to determine charge polarity. This is fragile — "Q" (quaternary ammonium) is an anion exchanger but the string "Q" doesn't match "anion" or "deae". Adding an explicit `charge_type: str = ""` field with values `"anion"`, `"cation"`, or `""` eliminates the heuristic and makes charge assignment unambiguous for all current and future IEX ligands.

---

## 3. Specific Process — Ligand Coupling Candidates

### 3.1 Q (Quaternary Ammonium) — Strong Anion Exchanger

**Goal:** Enable simulation of Q-type strong anion exchange media preparation.

**Scientific reason for inclusion:** Quaternary ammonium groups carry a permanent positive charge (non-titratable) across the entire working pH range (2-12). This makes Q media suitable for binding negatively charged proteins at high pH where weak exchangers like DEAE lose their charge (pKa ~11.5). Q Sepharose is one of the most widely used industrial IEX media for mAb polishing and viral clearance.

**Coupling chemistry:** The quaternary ammonium group is introduced by reacting an epoxide-activated agarose surface with a trimethylaminoethyl reagent. The epoxide ring-opening by the tertiary amine produces a stable secondary amine linkage with the quaternary ammonium at the distal end:

```
Agarose-O-CH2-CH(OH)-CH2-Epoxide + (CH3)3N-CH2-CH2-NH2
→ Agarose-O-CH2-CH(OH)-CH2-NH-CH2-CH2-N+(CH3)3
```

**Key parameters:**

| Parameter | Value | Basis |
|---|---|---|
| k_forward | 6e-5 m3/(mol*s) | Estimated from general epoxide-amine kinetics; quaternary amine reagents are moderately nucleophilic |
| E_a | 50,000 J/mol | Typical for epoxide ring-opening by amines at pH 10-11 |
| k_hydrol | 1e-5 /s | Epoxide hydrolysis rate at alkaline pH (competing reaction) |
| pH_opt | 10.5 | Optimal nucleophilicity of the amine without excessive hydrolysis |
| charge_type | "anion" | Permanent positive charge binds anions |

### 3.2 CM (Carboxymethyl) — Weak Cation Exchanger

**Goal:** Enable simulation of CM-type weak cation exchange media preparation.

**Scientific reason for inclusion:** Carboxymethyl groups (-CH2-COO-) carry a negative charge above pH ~4 (pKa of carboxyl). CM media are used for binding positively charged proteins under mild conditions. The titratable charge allows pH-based elution strategies. CM Sepharose is the standard weak cation exchanger in downstream processing.

**Coupling chemistry:** Carboxymethylation is achieved by reacting an amine-functional surface (from epoxide ring-opening by a short spacer amine) with chloroacetic acid, or by direct coupling of iminodiacetate-type reagents. In the M2 model, the simplified route treats the carboxymethyl reagent as an amine nucleophile attacking an epoxide:

```
Agarose-Epoxide + H2N-CH2-COOH → Agarose-CH2-CH(OH)-CH2-NH-CH2-COOH
```

**Key parameters:**

| Parameter | Value | Basis |
|---|---|---|
| k_forward | 3e-5 m3/(mol*s) | Glycine/aminoacetic acid is a moderate nucleophile |
| E_a | 45,000 J/mol | Slightly lower than Q due to smaller, more accessible reagent |
| charge_type | "cation" | Negative charge above pH 4 binds cations |

### 3.3 NTA (Nitrilotriacetic Acid) — IMAC Chelator

**Goal:** Enable simulation of Ni-NTA resin preparation for His-tag protein purification.

**Scientific reason for inclusion:** NTA is a tetradentate chelator (coordinates metal via 4 valencies) compared to IDA's tridentate coordination (3 valencies). This leaves 2 free coordination sites on NTA-Ni for His-tag binding vs 3 for IDA-Ni. The result: NTA has **higher selectivity** (less non-specific binding) but slightly **lower capacity** than IDA. Ni-NTA (Qiagen) is the industry standard for His-tag protein purification.

**Coupling chemistry:** NTA couples to epoxide-activated agarose via the central amine nitrogen:

```
Agarose-Epoxide + N(CH2-COOH)3 → Agarose-CH2-CH(OH)-CH2-N(CH2-COOH)3
```

After coupling, the NTA chelator is charged with Ni2+ (or Co2+) by incubation with metal salt solution. The metal charging step is not modeled in M2 Phase 1 — it is treated as implicit (the NTA profile represents the charged form).

**Key parameters:**

| Parameter | Value | Basis |
|---|---|---|
| k_forward | 2e-5 m3/(mol*s) | NTA is bulkier than IDA; slower diffusion to epoxide sites |
| MW | 191 Da | Larger than IDA (133 Da) — slightly reduced pore accessibility |
| functional_mode | "imac_chelator" | Same as IDA |

### 3.4 Butyl (Butylamine) — Mild HIC Ligand

**Goal:** Provide a lower-hydrophobicity HIC option for mild elution conditions.

**Scientific reason for inclusion:** In HIC chromatography, the ligand hydrophobicity determines the strength of protein-surface interaction. Phenyl groups (already in M2) are medium-high hydrophobicity and can cause irreversible binding of hydrophobic proteins. Butyl groups are less hydrophobic, enabling milder binding and elution conditions. The hydrophobicity ranking is: Butyl-S < Octyl < Butyl < Phenyl. Together, Phenyl + Butyl cover the practical HIC selectivity range.

**Coupling chemistry:** Identical to Phenyl — epoxide ring-opening by a primary alkyl amine:

```
Agarose-Epoxide + H2N-(CH2)3-CH3 → Agarose-CH2-CH(OH)-CH2-NH-(CH2)3-CH3
```

The C4 alkyl chain provides the hydrophobic interaction surface.

### 3.5 Glutathione — GST-Tag Affinity Ligand

**Goal:** Enable simulation of glutathione affinity media for GST-tagged protein purification.

**Scientific reason for inclusion:** Glutathione (GSH, a tripeptide: gamma-Glu-Cys-Gly, MW 307 Da) binds specifically to Glutathione S-Transferase (GST) — one of the most widely used fusion tags in recombinant protein production. GST-tagged proteins bind to immobilized glutathione with Kd ~0.1-1 uM and are eluted with free glutathione (10-20 mM). The GST-glutathione system is the third most popular tag-based purification after His-tag and Streptavidin-biotin.

**Coupling chemistry:** Glutathione couples to epoxide-activated agarose via its free amine (on the gamma-glutamate) or thiol (on the cysteine):

```
Agarose-Epoxide + GSH (via -NH2) → Agarose-CH2-CH(OH)-CH2-NH-GSH
```

At pH 8-9, the amine route dominates. At pH 7-8, the thiol route can be used for oriented coupling (thioether bond).

**Key parameters:**

| Parameter | Value | Basis |
|---|---|---|
| MW | 307 Da | Small molecule — full pore access |
| r_h | 0.5 nm | Tripeptide |
| functional_mode | "gst_affinity" | NEW mode for FunctionalMediaContract |
| spacer_key | "aha_spacer" | AHA spacer recommended for improved GST access |
| spacer_activity_multiplier | 1.15 | Moderate improvement from 10 A spacing |

### 3.6 Heparin — Dual-Mode Affinity + IEX Ligand

**Goal:** Enable simulation of heparin affinity media for growth factor and coagulation factor purification.

**Scientific reason for inclusion:** Heparin is a naturally occurring sulfated glycosaminoglycan (GAG) polysaccharide (avg MW 12-15 kDa) that binds a wide range of biologically important proteins through a combination of:

1. **Electrostatic interactions** — heparin is strongly negatively charged (sulfate + carboxylate groups) and acts as a cation exchanger
2. **Specific affinity** — many proteins have heparin-binding domains (HBDs) with specific spatial arrangements of basic residues

This dual-mode character makes heparin one of the most versatile affinity ligands in downstream processing. Target proteins include: FGF, VEGF, antithrombin III, lipases, DNA-binding proteins, and growth factors.

**Coupling chemistry:** Heparin is immobilized on epoxide-activated agarose via its free amino groups (from unsubstituted glucosamine residues):

```
Agarose-Epoxide + Heparin-NH2 → Agarose-CH2-CH(OH)-CH2-NH-Heparin
```

Due to heparin's large size (12-15 kDa), a DADPA spacer arm is recommended to reduce steric hindrance between the agarose surface and the target protein.

**Key parameters:**

| Parameter | Value | Basis |
|---|---|---|
| MW | ~14,000 Da | Polysaccharide — intermediate between small molecule and protein |
| r_h | 3.0 nm | Estimated from hydrodynamic data for heparin in solution |
| k_forward | 1e-5 m3/(mol*s) | Slower due to larger size and diffusion limitation |
| functional_mode | "heparin_affinity" | NEW mode |
| spacer_key | "dadpa_spacer" | DADPA recommended for large ligand |

---

## 4. Specific Process — Protein Coupling Candidates

### 4.1 Protein A/G Fusion — Broadest IgG Affinity

**Goal:** Enable simulation of Protein A/G fusion affinity media for maximum antibody species and subclass coverage.

**Scientific reason for inclusion:** Recombinant Protein A/G is an engineered fusion protein containing the IgG-binding domains from both *Staphylococcus aureus* Protein A (5 IgG-Fc binding domains) and *Streptococcus* Protein G (2-3 IgG-Fc binding domains). This fusion combines the strengths of both:

| Property | Protein A | Protein G | **Protein A/G** |
|---|---|---|---|
| Human IgG1 | Strong | Strong | **Strong** |
| Human IgG2 | Moderate | Strong | **Strong** |
| Human IgG3 | Weak | Strong | **Strong** |
| Human IgG4 | Strong | Strong | **Strong** |
| Mouse IgG1 | Weak | Strong | **Strong** |
| Mouse IgG2a | Strong | Strong | **Strong** |
| Rat IgG | Weak | Moderate | **Moderate** |
| Goat/Sheep IgG | Weak | Strong | **Strong** |

**Coupling chemistry:** Identical to Protein A — epoxide ring-opening by surface lysine amine groups. The larger MW (51 kDa vs 42 kDa for Protein A) means slightly lower molar coupling density but broader functionality.

**Key parameters:**

| Parameter | Value | Scientific Basis |
|---|---|---|
| MW | 51,000 Da | Fusion of A (42 kDa) + G domains (9 kDa insert) |
| r_h | 2.8 nm | Larger than A (2.5 nm) due to additional domains |
| activity_retention | 0.55 | Lower than individual A or G due to larger size and more random orientation |
| spacer_key | "dadpa_spacer" | Recommended for 51 kDa protein |
| spacer_activity_multiplier | 1.22 | 13 A DADPA spacer reduces steric hindrance for Fc domain access |
| max_surface_density | 2e-8 mol/m2 | Similar steric jamming to Protein A |
| is_macromolecule | True | Uses ligand_accessible_area (pore exclusion) |
| confidence_tier | "ranking_only" | Activity retention varies widely with coupling conditions |

**Spacer arm scientific reason:** At 51 kDa, the A/G fusion protein occupies a larger surface footprint than Protein A alone. The DADPA spacer (13 A) provides sufficient separation from the agarose surface for the Fc-binding domains to extend into solution and engage IgG molecules (~150 kDa, r_h ~5 nm). Without a spacer, multipoint attachment to the surface restricts domain mobility, reducing the fraction of protein molecules with accessible binding sites from ~55% to ~40%.

### 4.2 Streptavidin — Biotin-Tag Affinity

**Goal:** Enable simulation of streptavidin affinity media for biotin-tagged molecule capture.

**Scientific reason for inclusion:** The streptavidin-biotin interaction is the strongest non-covalent interaction known in biology (Kd ~10^-15 M, approaching covalent bond strength). Streptavidin is a homotetrameric protein (4 x 13 kDa = 53 kDa) from *Streptomyces avidinii*, with 4 biotin-binding sites per tetramer. This enables:

1. **Biotinylated protein capture** — any protein tagged with biotin can be purified on streptavidin resin
2. **Biotinylated antibody immobilization** — oriented antibody display on streptavidin surfaces
3. **Biotinylated nucleic acid capture** — DNA/RNA purification

The extreme binding affinity means elution requires harsh conditions (8M guanidine HCl, pH 1.5, or 70C heating). For applications requiring gentle elution, engineered streptavidin mutants (mSA, Strep-Tactin) with reduced affinity (Kd ~10^-7 M) are used.

**Coupling chemistry:** Streptavidin couples to epoxide-activated agarose via surface lysine residues. The tetrameric structure means all 4 biotin-binding pockets face different directions, providing inherent multi-directional accessibility regardless of coupling orientation.

**Key parameters:**

| Parameter | Value | Scientific Basis |
|---|---|---|
| MW | 53,000 Da | Homotetramer (4 x 13 kDa) |
| r_h | 2.8 nm | Compact globular tetramer |
| activity_retention | 0.70 | Higher than Protein A because all 4 binding sites are surface-exposed regardless of orientation |
| spacer_activity_multiplier | 1.22 | DADPA spacer improves biotin accessibility for biotinylated macromolecules |
| max_surface_density | 3e-8 mol/m2 | Slightly higher than Protein A (more compact tetramer) |
| functional_mode | "biotin_affinity" | NEW mode |

---

## 5. Specific Process — Linker Arm Candidates

### 5.1 DADPA (Diaminodipropylamine) — Industry Standard Protein Spacer

**Goal:** Provide a 13 A amine-terminated spacer for protein ligand immobilization.

**Scientific reason for inclusion:** DADPA is the spacer arm used in **EAH-Sepharose** (GE Healthcare/Cytiva), one of the most commercially important activated chromatography supports. The 9-atom chain (H2N-(CH2)3-NH-(CH2)3-NH2) provides:

1. **13 angstrom extension** — sufficient to move protein binding domains beyond the surface steric exclusion zone
2. **Internal secondary amine** — adds chain flexibility at the midpoint, allowing the distal end to sweep a larger solid angle than a rigid alkyl chain
3. **Terminal primary amine** — provides a nucleophilic group for subsequent coupling to the protein via crosslinkers (glutaraldehyde, EDC/NHS) or directly to a second epoxide layer

**How it works in the M2 Phase 1 model:**

DADPA is represented as a spacer profile (`dadpa_spacer`) with `spacer_activity_multiplier = 1.22`. When a protein coupling profile references `spacer_key = "dadpa_spacer"`, the coupling solver computes:

```python
effective_activity = activity_retention * 1.22
# Example: Protein A with DADPA
# effective_activity = 0.60 * 1.22 = 0.732
# Without DADPA:
# effective_activity = 0.60
```

This 22% improvement is the midpoint of the literature range (15-30%) for medium-length amine spacers on protein A resin.

### 5.2 AHA (6-Aminohexanoic Acid) — Industry Standard Acid-Terminated Spacer

**Goal:** Provide a 10 A carboxyl-terminated spacer for ligands requiring EDC/NHS coupling chemistry.

**Scientific reason for inclusion:** AHA is the spacer arm used in **NHS-activated Sepharose HP** (Cytiva), the gold standard for amine-reactive protein immobilization. The 7-atom chain (H2N-(CH2)5-COOH) provides:

1. **10 angstrom extension** — moderate spacing for mid-sized ligands (glutathione, heparin)
2. **Terminal carboxyl group** — enables activation with EDC/NHS to create an amine-reactive NHS ester, which then couples to protein amines under mild conditions (pH 7-8, room temperature)
3. **Methylene chain flexibility** — 5 CH2 groups provide moderate conformational freedom

**Key difference from DADPA:** AHA provides a **-COOH distal group** (acidic), while DADPA provides **-NH2** (basic). The choice depends on the coupling chemistry needed for the functional ligand. In Phase 1, this distinction is captured only as a different multiplier value (AHA = 1.15 vs DADPA = 1.22), reflecting AHA's slightly shorter length and the additional activation step required.

### 5.3 DAH (1,6-Diaminohexane) — Simple Amine Spacer

**Goal:** Provide a 9 A amine-terminated spacer as a simpler alternative to DADPA.

**Scientific reason for inclusion:** DAH is the spacer arm used in **AH-Sepharose 4B** (Cytiva). The 6-atom chain (H2N-(CH2)6-NH2) is:

1. **Simpler chemistry** — no internal amine, just a straight C6 diamine
2. **9 angstrom extension** — slightly shorter than DADPA (13 A)
3. **More hydrophobic** — the uninterrupted methylene chain is more hydrophobic than DADPA's amine-interrupted chain, which can increase non-specific binding in some applications

**When to choose DAH vs DADPA:**

| Property | DADPA (13 A) | DAH (9 A) |
|---|---|---|
| Length | Longer, better for large proteins | Shorter, adequate for medium proteins |
| Flexibility | Higher (internal amine hinge) | Lower (rigid alkyl chain) |
| Hydrophobicity | Lower (amine reduces hydrophobicity) | Higher (pure methylene chain) |
| Non-specific binding | Lower | Slightly higher |
| Activity multiplier | 1.22 | 1.08 |
| Best for | Protein A, A/G, streptavidin | NTA, IDA, small proteins |

---

## 6. Specific Process — FunctionalMediaContract Extension

### 6.1 Goal

Map functional ligand density from M2 to estimated chromatographic binding capacity (q_max) in M3 for all new ligand types.

### 6.2 Scientific Basis for q_max Mapping

The binding capacity of a chromatography column depends on:

```
q_max [mol/m3 bed] = functional_ligand_density [mol/m2] 
                     * specific_surface_area [m2/m3 bed]
                     * binding_stoichiometry [-]
```

Where:
- `specific_surface_area = 6*(1-eps_bed) / d_particle` (spherical particles)
- `binding_stoichiometry` = number of target molecules bound per functional ligand

### 6.3 Mapping Table

| Ligand Type | Binding Stoichiometry | Scientific Reason |
|---|---|---|
| IEX (DEAE, Q, SP, CM) | Charge-dependent; ~1.0 per accessible site at low ionic strength | One protein molecule binds to multiple charge patches; effective stoichiometry depends on protein surface charge |
| HIC (Phenyl, Butyl) | Not directly mapped to q_max | HIC capacity depends on salt concentration; ligand density alone is insufficient. Mapped as `"not_mapped"` |
| IMAC (IDA, NTA) | ~1.0 per metal-chelator site (for His-tagged protein) | Each NTA-Ni coordinates one His-tag (6xHis) |
| Affinity (Protein A/G) | 2.0 per Protein A molecule | Protein A has 5 Fc-binding domains but steric constraints limit effective binding to ~2 IgG per Protein A |
| Affinity (Streptavidin) | 4.0 per streptavidin tetramer | 4 biotin-binding sites per tetramer, each can capture one biotinylated molecule |
| Affinity (Glutathione) | 1.0 per glutathione molecule | One GST domain binds one glutathione |
| Affinity (Heparin) | ~1.0 (approximate, variable) | Heparin-binding stoichiometry varies widely by protein; 1.0 is a first approximation |

---

## 7. Specific Process — UI Integration

### 7.1 Goal

Enable users to select all new ligands, proteins, and spacers from the M2 tab dropdowns without code changes.

### 7.2 Dropdown Extensions

**Ligand Coupling dropdown** expands from 4 to 10 options:

```
Current:  DEAE, IDA, Phenyl, Sulfopropyl
Add:      Q, CM, NTA, Butyl, Glutathione, Heparin
```

**Protein Coupling dropdown** expands from 2 to 4 options:

```
Current:  Protein A, Protein G
Add:      Protein A/G Fusion, Streptavidin
```

### 7.3 Spacer Selectbox (New UI Element)

A new optional selectbox appears when Ligand Coupling or Protein Coupling is selected:

```
Spacer Arm: [None (Direct Coupling)] [DADPA (13 A)] [AHA (10 A)] [DAH (9 A)]
```

**Scientific reason for UI placement:** The spacer choice affects the physical distance between the agarose surface and the functional ligand. For small-molecule ligands (IEX, HIC, IMAC chelators), direct coupling is typically sufficient because these molecules have full pore accessibility. For protein ligands, the spacer improves activity retention by 10-30%. Making the spacer optional (default: "None") preserves backward compatibility while giving advanced users control.

**Implementation:** When a spacer is selected, its `spacer_activity_multiplier` overrides the coupling profile's default. This is applied in the coupling solver as:

```python
effective_activity = min(activity_retention * spacer_activity_multiplier, 1.0)
```

The `min(..., 1.0)` cap prevents the physically impossible case of >100% activity retention.

---

## 8. Implementation Roadmap

### 8.1 Phase Overview

| Phase | Work Nodes | Description | Model Tier |
|---|---|---|---|
| A | WN-1 | Add 4 new ReagentProfile fields (spacer_key, spacer_length_angstrom, spacer_activity_multiplier, charge_type) | Opus |
| B | WN-2 + WN-3 (parallel) | Add 8 coupling profiles + 3 spacer profiles | Sonnet |
| C | WN-4 + WN-5 + WN-6 (parallel) | Apply spacer multiplier in solvers + extend FMC + extend UI | Sonnet |
| D | WN-7 + WN-8 + WN-9 (parallel) | Tests for profiles, solvers, and FMC | Sonnet |
| E | WN-10 | Docstring updates | Haiku |

### 8.2 File Changes

| File | Changes | LOC |
|---|---|---|
| `reagent_profiles.py` | 4 new fields + 8 coupling + 3 spacer profiles + backfill charge_type on existing | +280 |
| `modification_steps.py` | Apply spacer_activity_multiplier in 2 solvers | +12 |
| `orchestrator.py` | Extend _mode_map, add q_max branches, fix IEX detection | +45 |
| `app.py` | Extend dropdowns, add spacer selectbox | +40 |
| `tests/` | Profile validation, spacer multiplier, FMC mapping | +190 |
| **Total** | | **~567** |

---

## 9. Validation and Acceptance

### 9.1 Scientific Validation

| Check | Method |
|---|---|
| All new profiles have physically reasonable kinetic parameters | Compare k_forward with existing profiles of same chemistry class |
| Spacer multipliers within literature range | 1.0 <= multiplier <= 1.40; protein profiles use 1.08-1.30 |
| Charge_type correctly assigned | DEAE/Q = "anion"; SP/CM = "cation"; others = "" |
| FMC q_max estimates in reasonable range | Compare with published resin specifications (e.g., Protein A: 30-50 mg IgG/mL) |
| activity_retention * multiplier never exceeds 1.0 | Assertion in test |

### 9.2 Computational Validation

| Check | Method |
|---|---|
| All 25 profiles instantiate without error | Unit test |
| All profiles have required metadata (confidence_tier, calibration_source, hazard_class) | Parametrized test |
| Spacer multiplier applied correctly in ligand + protein coupling | End-to-end coupling test with and without spacer |
| FMC maps all new functional_mode values | Test each new mode produces non-zero q_max |
| Existing 14 profiles unchanged | Regression test comparing field values |
| Full M1 → M2 → M3 pipeline with new profiles | Integration test |
| 142 existing tests still pass | Full test suite run |

### 9.3 UI Validation

| Check | Method |
|---|---|
| All 10 ligand options appear in dropdown | Visual check |
| All 4 protein options appear in dropdown | Visual check |
| Spacer selectbox appears for coupling steps | Visual check |
| Spacer profiles do NOT appear in step-type dropdowns | Filter test |
| Switching reagent within same chemistry resets parameters | Widget key test |

---

> **Disclaimer**: This integration plan is provided for informational, research, and development purposes only. All kinetic parameters are order-of-magnitude estimates from general reaction chemistry principles. Spacer arm improvement factors are empirical estimates requiring experimental calibration for specific materials. The author is an AI assistant and this plan should be treated as a structured starting point for further investigation and implementation.
