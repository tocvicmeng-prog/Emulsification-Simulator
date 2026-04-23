# Appendix J — Functionalisation Wet-Lab Protocols

*First Edition supplement to the Polysaccharide-Based Microsphere Emulsification
Simulator user manual. Written for first-time users: downstream-processing
technicians or junior researchers who have never done affinity chromatography.
Every protocol is self-contained. You do not need to read the rest of the manual
to use the page in front of you.*

**How to use this appendix**

1. Decide what you want to attach: small-molecule ligand (§J.2), protein (§J.3),
   a metal ion for IMAC (§J.5). That is the "coupling" step.
2. Pick the activation chemistry that your coupling step needs (§J.1).
3. Decide whether you need a spacer arm (§J.4) before the ligand.
4. If the ligand is a protein, pretreat it (§J.6).
5. Execute the chosen activation, spacer, and coupling protocols. Wash (§J.7)
   between every step. Quench unreacted sites at the end (§J.8).
6. Store the finished gel in the recommended buffer (§J.7, end).

**Conventions used throughout**

- **"CV"** = column volume. One CV of a 1 mL gel bed = 1 mL of buffer.
- **"wet gel"** = drained polysaccharide matrix (agarose, sepharose, cellulose).
  When a protocol says "10 mL wet gel", it means gel drained to a damp cake
  with no free supernatant. Weigh or measure by volume in a graduated cylinder.
- **"Evidence tier"** mirrors the simulator's `ModelEvidenceTier` classification:
    - **VALIDATED** — protocol is reproducible in multiple peer-reviewed
      sources with expected yields quantified.
    - **SEMI_QUANTITATIVE** — protocol works reliably; numerical outcome
      (coupling density, activity) depends on matrix lot, ligand, and scale.
    - **QUALITATIVE_TREND** — directional only; yields vary by 10× or more.
- **"EmulSim M2 key"** tells you which dropdown option in the Functionalisation
  tab invokes the matching simulation model. Some entries say "not yet wired";
  that means the backend supports the chemistry, but the UI dropdown for it
  has not been added yet. Use the protocol anyway; record your results; they
  will calibrate the model in a future release.
- **Temperature** is in °C. **"RT"** = room temperature ≈ 20-25 °C unless
  the protocol specifies otherwise.
- **pH** is measured with a calibrated pH meter, not indicator strips, unless
  indicated.
- **Water** means ultrapure ≥18 MΩ·cm (e.g., MilliQ), not tap or DI water,
  for every protocol.

**Universal safety reminder** (read once, remember every time)

- Nitrile gloves, safety glasses, and a buttoned lab coat are the minimum PPE
  for every protocol in this appendix. Long trousers, closed-toe shoes.
- Work in a certified fume hood for any step that evolves fumes, handles
  volatile organics, or uses reagents marked "fume hood required" below.
- Every reagent box below has GHS pictograms and H-codes. If you cannot
  identify a pictogram, read the SDS for that reagent **before** you handle it.
- Waste: never mix halogenated, non-halogenated, aqueous, and solid waste.
  Each protocol names its waste stream. If in doubt, ask your lab's safety
  officer before disposing.
- If you spill, splash, or inhale anything: stop, rinse skin/eyes for 15
  minutes at the nearest safety shower or eyewash, notify a supervisor,
  consult the SDS, do not finish the protocol without clearance.

**How this appendix relates to the M1 Suggestion Derivation pages (v9.2.0)**

After every M1 run the simulator produces optimization suggestions
with **[📊 derivation]** hyperlinks (see §9 of the main manual). Those
derivation pages answer *why* the system recommended a parameter
change and *by how much*. This appendix answers *how to execute the
change at the bench*. Two concrete cross-references:

- If the M1 derivation page for **"adjust cooling rate"** reports a
  rate the passive bath cannot deliver, the STMP two-phase protocol
  (§J.1.7) is the canonical wet-lab strategy for decoupling reagent
  loading from reaction activation — its "cold load / hot activate"
  temperature profile sidesteps the passive-cooling ceiling.
- If the derivation page for **"increase crosslinker"** reports a
  required conversion p close to the stoichiometric ceiling, the
  hydroxyl-activation protocols (§J.1.1 – §J.1.6) are alternatives
  that target the agarose -OH network instead of the chitosan -NH₂
  network, giving you a second crosslinkable surface to draw from.

Each protocol card lists an **EmulSim M2 key** (e.g. `stmp`, `dvs`,
`ech`) that lines up with the simulator's Module 2 backend and —
when applicable — with the key shown on the derivation page URL.

---

## J.1 Hydroxyl Activation

Polysaccharide matrices (agarose, sepharose, cellulose, dextran) present -OH
groups on every sugar ring. These -OH groups are chemically inert under mild
conditions. Activation converts -OH into a leaving group or a reactive
electrophile so that downstream coupling (amine, thiol, hydroxyl, carboxyl)
can proceed at room temperature and near-neutral pH.

Pick **one** activator. More aggressive activators = higher ligand-coupling
density but also more matrix damage and more residual reactivity after
coupling (so a bigger quench step at the end).

---

### J.1.1  Hydroxyl Activation — CNBr (cyanogen bromide)

**Purpose:** Activate -OH on agarose / sepharose so that primary amines can
couple directly. The classical gold-standard protocol since 1967.

**Evidence tier:** VALIDATED

**EmulSim M2 key:** `cnbr_activation` (not yet wired in v8.3.5 UI; backend
ranking-only tier).

**Based on:** Axén, Porath & Ernback (*Nature* 1967, 214:1302); Kohn & Wilchek
(*Appl. Biochem. Biotechnol.* 1984, 9:285); Hermanson *Bioconjugate Techniques*
3rd ed. Ch. 15.3.1.

#### Safety — READ BEFORE HANDLING

CNBr is **highly toxic**. Most labs buy sepharose pre-activated with CNBr
("CNBr-activated Sepharose 4B", Cytiva) instead of activating in-house.
Follow the Cytiva protocol if you use pre-activated gel.

If you must activate in-house:
- GHS pictograms: skull-and-crossbones, corrosive, environmental hazard.
- H-codes: H300 (fatal if swallowed), H311 (toxic in contact with skin),
  H330 (fatal if inhaled), H314 (causes severe skin burns and eye damage),
  H410 (very toxic to aquatic life).
- PPE: double nitrile gloves, face shield over safety glasses, chemical-
  resistant apron, lab coat. Fume hood mandatory, never handle on the bench.
- CNBr releases HCN on hydrolysis. Keep a cyanide test kit and a supply of
  sodium thiosulfate + sodium hydroxide on hand for accidental decomposition.
- Waste: collect all CNBr-contaminated aqueous waste in a labelled carboy
  containing 1 M NaOH. Let it neutralise for ≥ 24 h before disposal via your
  institution's hazardous-chemical-waste stream. Solid waste (gloves, towels)
  goes in a sealed container labelled "CNBr contaminated".

**Strong recommendation:** use pre-activated sepharose. The in-house
activation is documented here only for lab completeness.

#### What you need (pre-activated-sepharose route)

- 10 mL CNBr-activated Sepharose 4B, supplied as freeze-dried powder
- 200 mL cold 1 mM HCl (ice bath)
- 50 mL coupling buffer: 0.1 M NaHCO₃, 0.5 M NaCl, pH 8.3
- Sintered glass funnel (porosity 2, 40-100 μm), vacuum flask, vacuum source
- End-over-end rotator
- pH meter

#### Procedure (pre-activated-sepharose route)

1. Put on the PPE listed above. Confirm the fume hood sash is at the
   recommended working height.
2. Weigh 3 g of freeze-dried CNBr-activated Sepharose 4B powder into a 50 mL
   conical tube. This swells to ≈10 mL wet gel.
3. Add 30 mL cold 1 mM HCl. Cap the tube. Invert gently 5 times to disperse.
4. Pour the suspension into the sintered glass funnel. Apply gentle vacuum to
   drain. Wash with 6 × 30 mL 1 mM HCl (total 200 mL). Keep on ice. Do **not**
   let the gel dry out — stop vacuum when the surface just loses its shine.
5. Immediately transfer the damp gel to a fresh 50 mL tube and add the ligand
   solution in coupling buffer (§J.2). Proceed to coupling within 15 min. If
   you cannot couple immediately, transfer to 1 mM HCl and store at 4 °C for
   no more than 2 h — activity decays rapidly at neutral pH.

#### Quality control

- Activation density: 12-16 μmol reactive sites per mL wet gel (specification
  from Cytiva; each lot has a Certificate of Analysis).
- Visual: gel is off-white and translucent; any brown or yellow colour
  indicates hydrolysis. Reject.

#### Troubleshooting

- **Low coupling yield (< 50 %):** gel was over-washed or sat too long in
  coupling buffer before ligand addition. Move faster.
- **Ligand activity lost after coupling:** CNBr creates an isourea bond that
  is slightly positively charged at neutral pH. For charge-sensitive ligands,
  use CDI or BDDE (§J.1.5, §J.1.3) instead.
- **Matrix turns yellow/brown:** hydrolysis. Start over with fresh gel.

#### References

Axén R, Porath J, Ernback S (1967) *Nature* 214:1302. Kohn J, Wilchek M (1984)
*Appl. Biochem. Biotechnol.* 9:285. Hermanson GT (2013) *Bioconjugate
Techniques* 3rd ed., Academic Press, Ch. 15.3.1.

---

### J.1.2  Hydroxyl Activation — Epichlorohydrin

**Purpose:** Activate -OH groups to epoxides. The epoxide can couple amines,
thiols, and hydroxyls directly. Also introduces a small 3-atom spacer.

**Evidence tier:** VALIDATED

**EmulSim M2 key:** `epichlorohydrin` (backend: `qualitative_trend`; CMR
reagent, institution approval usually required).

**Based on:** Matsumoto et al. (*J. Chromatogr.* 1980, 188:457); Sundberg &
Porath (*J. Chromatogr.* 1974, 90:87); Hermanson Ch. 15.3.4.

#### Safety — READ BEFORE HANDLING

- GHS pictograms: skull-and-crossbones, corrosive, health hazard (CMR).
- H-codes: H301 (toxic if swallowed), H311 (toxic in contact with skin),
  H331 (toxic if inhaled), H314 (severe burns), H317 (skin sensitisation),
  H350 (may cause cancer), H341 (suspected mutagen), H361 (suspected
  reproductive toxin).
- PPE: double nitrile gloves (change every 30 min; epichlorohydrin permeates
  nitrile within ~1 h), face shield, chemical-resistant apron, lab coat.
  Fume hood mandatory.
- Epichlorohydrin is a confirmed animal carcinogen and probable human
  carcinogen (IARC Group 2A). Some institutions require a written risk
  assessment and supervisor sign-off before use.
- Waste: halogenated organic aqueous waste. Collect in a dedicated labelled
  carboy. Never pour down the drain.

#### What you need (per 10 mL wet gel)

| Reagent | Amount | CAS | Storage |
|---|---|---|---|
| Epichlorohydrin | 2 mL | 106-89-8 | RT, tightly capped, fume hood cabinet |
| 2 M NaOH | 10 mL | 1310-73-2 | RT |
| Sodium borohydride (NaBH₄) | 20 mg | 16940-66-2 | Desiccator, RT |
| Water (for washes) | ≥ 200 mL | — | RT |

Equipment: 50 mL screw-cap polypropylene tube; end-over-end rotator; water
bath at 40 °C; sintered glass funnel; pH meter.

#### Procedure

1. PPE on. Open the fume hood sash. Verify airflow (test with a tissue at the
   sash face — it should draw inward).
2. Wash 10 mL wet gel with 5 × 30 mL water on a sintered glass funnel. Drain
   to a damp cake.
3. Transfer drained gel to the 50 mL screw-cap tube. Add 10 mL 2 M NaOH and
   20 mg NaBH₄ (the NaBH₄ suppresses Maillard-like side reactions).
4. In the fume hood, add 2 mL epichlorohydrin. Cap tightly. Mix by inversion
   (10 inversions).
5. Rotate end-over-end at 40 °C, 2 h.
6. Transfer to the sintered glass funnel. Wash with 10 × 30 mL water until
   the filtrate pH is ≤ 8 and the smell of epichlorohydrin is gone (confirm
   with a clean nose-wave from above the funnel — do NOT lean in; fume hood
   is still on).
7. Drain to a damp cake. The gel is now epoxide-activated and ready for
   coupling (§J.2.2).

#### Quality control

- Epoxide density: 15-25 μmol epoxide per mL wet gel. Titrate by reacting a
  100 μL aliquot with 2 mL 1.3 M Na₂S₂O₃ pH 7.0 for 30 min at RT, then
  back-titrating the liberated NaOH with 0.1 M HCl to pH 7.0 (Sundberg &
  Porath 1974).

#### Troubleshooting

- **Gel becomes yellow / brown:** over-reaction or NaBH₄ missing. Restart.
- **Low epoxide density:** epichlorohydrin bottle open too long (hydrolysed).
  Use a fresh, unopened bottle. Epichlorohydrin hydrolyses on air exposure.
- **Coupling fails after activation:** activated gel was stored. Epoxide
  hydrolyses in water over 24 h. Couple the same day.

#### References

Matsumoto I, Mizuno Y, Seno N (1980) *J. Chromatogr.* 188:457. Sundberg L,
Porath J (1974) *J. Chromatogr.* 90:87. Hermanson (2013) Ch. 15.3.4.

---

### J.1.3  Hydroxyl Activation — BDDE (1,4-butanediol diglycidyl ether, "bisoxirane")

**Purpose:** Activate -OH with a bifunctional epoxide. Adds a 12-atom
flexible spacer in the same step. The free end couples to amine / thiol /
hydroxyl ligands. Milder than epichlorohydrin, no halogen.

**Evidence tier:** SEMI_QUANTITATIVE

**EmulSim M2 key:** `bdde` (not yet wired).

**Based on:** Sundberg & Porath (*J. Chromatogr.* 1974, 90:87); Hermanson Ch. 15.3.5.

#### Safety

- GHS pictograms: exclamation mark, health hazard.
- H-codes: H315 (skin irritant), H319 (eye irritant), H317 (skin sensitiser).
- Not a carcinogen in the epichlorohydrin sense, but a known sensitiser —
  repeated exposure causes contact dermatitis.
- PPE: nitrile gloves, safety glasses, lab coat. Fume hood recommended.
- Waste: non-halogenated organic-aqueous stream.

#### What you need (per 10 mL wet gel)

| Reagent | Amount | CAS | Storage |
|---|---|---|---|
| BDDE | 1 mL | 2425-79-8 | 4 °C, tightly capped |
| 0.6 M NaOH | 10 mL | 1310-73-2 | RT |
| NaBH₄ | 20 mg | 16940-66-2 | Desiccator |

#### Procedure

1. PPE on. Fume hood recommended but not mandatory.
2. Drain 10 mL wet gel on the sintered funnel, wash with 5 × 30 mL water.
3. Transfer to a 50 mL tube. Add 10 mL 0.6 M NaOH and 20 mg NaBH₄.
4. Add 1 mL BDDE. Cap, invert 10 times.
5. Rotate at 25 °C (RT), 8 h (or overnight).
6. Wash with 10 × 30 mL water, then 3 × 30 mL acetone (to remove unreacted
   BDDE), then 3 × 30 mL water.
7. Drain to damp cake. Ready for coupling (§J.2.2).

#### Quality control

- Epoxide density: 10-20 μmol/mL wet gel (lower than epichlorohydrin but
  steric access is better thanks to the spacer).

#### Troubleshooting

- **Low density:** NaOH too dilute or too cold. Verify NaOH concentration
  with a standard HCl titration.
- **Two-sided crosslinking (gel shrinks):** too much BDDE relative to matrix.
  Halve the BDDE next time.

#### References

Sundberg & Porath (1974) *J. Chromatogr.* 90:87. Hermanson (2013) Ch. 15.3.5.

---

### J.1.4  Hydroxyl Activation — DVS (divinyl sulfone)

**Purpose:** Activate -OH with a vinyl sulfone. Reacts via Michael addition
with amines (fast at pH 9-11), thiols (fast at pH 7-9), and hydroxyls
(slow, pH 10-12). One-step coupling from vinyl sulfone.

**Evidence tier:** SEMI_QUANTITATIVE

**EmulSim M2 key:** `dvs` (backend: `ranking_only` in v8.3.5).

**Based on:** Porath et al. (*J. Chromatogr.* 1975, 103:49); Mateo et al.
(*Enzyme Microb. Technol.* 2007, 39:274); Hermanson Ch. 15.3.6.

#### Safety

- GHS pictograms: skull-and-crossbones, corrosive, health hazard.
- H-codes: H300 (fatal if swallowed), H330 (fatal if inhaled), H314 (severe
  burns), H317 (strong skin sensitiser — a small repeated exposure can
  produce a persistent allergy).
- PPE: double nitrile gloves (change hourly), face shield, chemical apron,
  lab coat. Fume hood mandatory.
- DVS is volatile. Keep tightly capped. Do not weigh on an open balance
  outside the hood.
- Waste: non-halogenated organic-aqueous. Before disposal, quench by adding
  excess ethanolamine (1 M) and letting stand ≥ 4 h.

#### What you need (per 10 mL wet gel)

| Reagent | Amount | CAS | Storage |
|---|---|---|---|
| DVS | 1 mL | 77-77-0 | 4 °C, ampoule |
| 0.5 M Na₂CO₃ pH 11 | 10 mL | — | RT (freshly prepared) |

#### Procedure

1. PPE on. Fume hood sash at working height.
2. Drain 10 mL wet gel, wash with 5 × 30 mL water, then 2 × 30 mL 0.5 M Na₂CO₃.
3. Transfer to a 50 mL tube. Add 10 mL 0.5 M Na₂CO₃ pH 11.
4. In the hood, add 1 mL DVS dropwise via 1 mL serological pipette. Cap.
5. Rotate end-over-end at RT, 70 min. Exact time matters — over-activation
   produces crosslinks that reduce coupling density.
6. Wash rapidly with 10 × 30 mL water on the sintered funnel. Drain damp.
7. Gel is now vinyl-sulfone activated. Proceed to coupling (§J.2.3) within
   30 minutes — activated sites hydrolyse over hours at neutral pH.

#### Quality control

- Vinyl-sulfone density: 20-40 μmol/mL wet gel (titrate by reaction with
  excess cysteine, then quantify free thiol with Ellman's reagent; Hermanson
  Ch. 15.3.6 details).

#### Troubleshooting

- **Gel crosslinks / shrinks during activation:** DVS concentration too high
  or reaction time too long. Reduce DVS to 0.5 mL and time to 60 min.
- **Coupling fails:** activated gel was stored > 30 min before use. Plan
  ligand to be ready before starting activation.

#### References

Porath J, Låås T, Janson J-C (1975) *J. Chromatogr.* 103:49. Mateo C et al.
(2007) *Enzyme Microb. Technol.* 39:274. Hermanson (2013) Ch. 15.3.6.

---

### J.1.5  Hydroxyl Activation — CDI (1,1'-carbonyldiimidazole)

**Purpose:** Activate -OH to an imidazole carbamate. Couples primary amines
at pH 9-10 forming a stable N-alkyl carbamate. Milder than CNBr, no charged
isourea. Good for charge-sensitive ligands.

**Evidence tier:** VALIDATED

**EmulSim M2 key:** `cdi` (not yet wired; protein-friendly successor to CNBr).

**Based on:** Hearn (*Meth. Enzymol.* 1987, 135:102); Bethell et al.
(*J. Biol. Chem.* 1979, 254:2572); Hermanson Ch. 15.3.2.

#### Safety

- GHS pictograms: exclamation mark.
- H-codes: H315 (skin irritant), H319 (eye irritant), H335 (respiratory
  irritant). No CMR classification.
- PPE: nitrile gloves, safety glasses, lab coat. Fume hood recommended for
  powder weighing.
- CDI is moisture-sensitive. Store tightly capped in a desiccator. Weigh
  quickly.
- Waste: non-halogenated organic-aqueous.

#### What you need (per 10 mL wet gel)

| Reagent | Amount | CAS | Storage |
|---|---|---|---|
| CDI | 200 mg | 530-62-1 | –20 °C, desiccated |
| Anhydrous acetone | 30 mL | 67-64-1 | RT, over molecular sieves |
| Anhydrous dioxane | 30 mL | 123-91-1 | RT; peroxide-tested within 12 months |

#### Procedure

1. PPE on.
2. Drain 10 mL wet gel on the sintered funnel. Wash with 5 × 30 mL water,
   then dehydrate with an acetone gradient: 2 × 30 mL 30 % acetone, 2 × 30 mL
   60 % acetone, 3 × 30 mL 100 % acetone, 2 × 30 mL anhydrous acetone, then
   2 × 30 mL anhydrous dioxane. **This step is essential: CDI hydrolyses
   instantly in water.** The gel must be in a water-free solvent.
4. Transfer drained dioxane-exchanged gel to a 50 mL tube. Add 10 mL
   anhydrous dioxane.
5. Add 200 mg CDI (weigh in the fume hood, transfer quickly).
6. Rotate at RT, 30-60 min. CO₂ evolves visibly (small bubbles).
7. Wash with 5 × 30 mL anhydrous dioxane, then reverse the gradient back to
   water: 2 × 30 mL 100 % acetone, 2 × 30 mL 60 % acetone, 2 × 30 mL 30 %
   acetone, 2 × 30 mL water.
8. The gel is now imidazole-carbamate-activated. Proceed to coupling (§J.2).

#### Quality control

- Activation density: 20-80 μmol/mL wet gel (titrate by release of imidazole
  at 280 nm after reaction with excess glycine).

#### Troubleshooting

- **Almost no activation:** gel was not fully dehydrated. Check the gradient
  wash step; the last acetone wash should pour off clear and the gel should
  be visibly less swollen.
- **Gel shrinks catastrophically:** rapid dehydration. Slow the acetone
  gradient next time (smaller steps).

#### References

Hearn MTW (1987) *Meth. Enzymol.* 135:102. Bethell GS, Ayers JS, Hancock WS,
Hearn MTW (1979) *J. Biol. Chem.* 254:2572. Hermanson (2013) Ch. 15.3.2.

---

### J.1.6  Hydroxyl Activation — Tresyl chloride (2,2,2-trifluoroethanesulfonyl chloride)

**Purpose:** Activate -OH to a tresyl ester. Couples primary amines under
very mild conditions (pH 7-9, RT). The linkage is a secondary amine with
no charge. Ideal for labile / charge-sensitive proteins.

**Evidence tier:** SEMI_QUANTITATIVE

**EmulSim M2 key:** `tresyl` (not yet wired).

**Based on:** Nilsson & Mosbach (*Meth. Enzymol.* 1984, 104:56); Hermanson
Ch. 15.3.3.

#### Safety

- GHS pictograms: corrosive, exclamation mark.
- H-codes: H314 (severe burns), H335 (respiratory irritant), H330 (fatal
  if inhaled — mainly from the HCl it releases on hydrolysis).
- PPE: double nitrile gloves, face shield, chemical apron, lab coat. Fume
  hood mandatory.
- Store under dry argon or nitrogen. Reacts violently with water.
- Waste: halogenated organic-aqueous. Neutralise aqueous rinses with solid
  sodium bicarbonate before pooling.

#### What you need (per 10 mL wet gel)

| Reagent | Amount | CAS | Storage |
|---|---|---|---|
| Tresyl chloride | 100 μL | 1648-99-3 | 4 °C, under dry inert gas |
| Anhydrous acetone | 30 mL | 67-64-1 | RT, over sieves |
| Pyridine | 100 μL | 110-86-1 | RT (peroxide-tested) |

#### Procedure

1. PPE on. Fume hood mandatory.
2. Dehydrate 10 mL wet gel through an acetone gradient (as in §J.1.5 steps
   2-3) ending in anhydrous acetone.
3. Transfer dehydrated gel to a 50 mL tube with 10 mL anhydrous acetone.
4. Add 100 μL pyridine (acid scavenger), then 100 μL tresyl chloride dropwise
   in the hood with gentle swirling.
5. Rotate at 0-4 °C (on ice), 10 min.
6. Wash in the funnel: 5 × 30 mL ice-cold acetone (not water — hydrolysis),
   then 5 × 30 mL ice-cold 1 mM HCl, then 3 × 30 mL cold coupling buffer.
7. Proceed to coupling (§J.2) immediately. Tresyl-activated matrix has a
   half-life of ≈ 1 h at 4 °C, pH 7.

#### Quality control

- Activation density: 10-30 μmol/mL wet gel (approximated from subsequent
  amine coupling yield; direct titration is difficult).

#### Troubleshooting

- **Very low coupling yield:** tresyl chloride bottle was opened before. It
  absorbs moisture rapidly. Use a fresh ampoule.
- **Gel turns yellow:** pyridine is discoloured (old). Use fresh pyridine.

#### References

Nilsson K, Mosbach K (1984) *Meth. Enzymol.* 104:56. Hermanson (2013) Ch. 15.3.3.

---

### J.1.7  Hydroxyl Activation — STMP (Sodium Trimetaphosphate, triggerable dual-network)

**Purpose:** Crosslink **both** agarose hydroxyls (dominant, phosphate diester)
and chitosan primary amines (secondary, phosphoramide) in the *same* bead, in a
thermally- and pH-triggerable two-phase reaction. Phase A loads STMP into the
pre-gelled bead at 4 °C / pH 7 (no reaction). Phase B raises to 60 °C / pH 11
and the ring-opened trimetaphosphate crosslinks the polymer networks. The
two-phase design gives a uniform radial crosslink profile (Thiele modulus ~0.35
for a 250 µm bead radius) — the "dip in acid TPP" ionic-gelation approach gives
a skin-core structure instead.

**Do not confuse STMP with STPP.** Sodium **Tri**metaphosphate (STMP,
Na₃P₃O₉, CAS **7785-84-4**) is the *cyclic* trimer used here — covalent,
alkaline pH, triggerable. Sodium **Tripolyphosphate** (STPP, Na₅P₃O₁₀, CAS
7758-29-4) is the *linear* ionic crosslinker in EmulSim's `tpp` entry —
different chemistry, different pH window, reversible. Check the CAS on every
reagent bottle before you weigh anything.

**Evidence tier:** SEMI_QUANTITATIVE for **both** the agarose-OH phosphate
diester path (well-established in starch/cellulose phosphorylation
literature) and, as of v9.2.2, the chitosan-NH₂ phosphoramide side-reaction
(rate constants calibrated to Seal 1996 + Salata 2015, with the parallel
NH₂ ODE track wired into the L3 solver). The simulator now reports the
two contributions to the chitosan-network modulus separately in
`CrosslinkingResult.G_chit_diester` and `G_chit_phosphoramide`.

**EmulSim keys:** `stmp` (L3 primary crosslinker); `stmp_secondary` (M2
secondary crosslinker). Both map to `mechanism="hydroxyl"` and route through
the existing hydroxyl-covalent second-order solver.

**Based on:** Lim & Seib (*Cereal Chem.* 1993, 70:137); Kasemsuwan & Jane
(*Cereal Chem.* 1996, 73:702); Lack et al. (*Carbohydr. Res.* 2004, 339:2391);
Seal (*Biomaterials* 1996, 17:1869); Salata et al. (*Int. J. Biol. Macromol.*
2015, 81:1009); SA-EMULSIM-XL-002 Rev 0.1.

#### Safety

- GHS pictograms: exclamation mark (mild irritant). Food-grade additive E452.
- H-codes: H315 (skin irritation), H319 (eye irritation). No acute-toxicity
  or carcinogenicity codes.
- PPE: single nitrile gloves, safety glasses, lab coat. Bench work is fine;
  fume hood not required for STMP itself, but Phase B uses NaOH/Na₂CO₃
  (corrosive) — wear a face shield when pipetting the pH 11 buffer hot.
- Waste: aqueous-neutral after quenching with HCl. Dispose as non-halogenated
  aqueous waste. The quenched solution contains orthophosphate — if your
  institution caps phosphate discharge, collect for chemical waste instead.

#### What you need (per 10 mL wet gel)

| Reagent | Amount | CAS | Storage |
|---|---|---|---|
| Sodium trimetaphosphate (Na₃P₃O₉) | 200 mg | 7785-84-4 | RT, dry |
| 0.1 M HEPES pH 7.0 | 20 mL | 7365-45-9 | 4 °C (freshly prepared) |
| 0.5 M Na₂CO₃ + 0.5 M NaOH pH 11.0 | 15 mL | 497-19-8 / 1310-73-2 | RT (freshly prepared, ≤ 4 h) |
| 0.1 M HCl (ice-cold) | 20 mL | 7647-01-0 | 4 °C |
| 10 mM EDTA pH 8 (wash) | 15 mL | 60-00-4 | RT |

#### Procedure

1. PPE on. Bench work.
2. Drain 10 mL wet gel on a sintered funnel, wash with 3 × 20 mL water.
3. **Phase A — cold loading.** Transfer gel into a 50 mL tube. Add 10 mL
   cold 0.1 M HEPES pH 7.0 and dissolve 200 mg STMP in the tube
   (final 2 % w/v). Rotate end-over-end at 4 °C, 30 min. Expected state:
   STMP is uniformly distributed throughout the bead volume; negligible
   reaction has occurred (pH 7, 4 °C).
4. Drain quickly on the sintered funnel. Do **not** wash at this stage —
   washing removes the uniformly-loaded STMP.
5. **Phase B — hot alkaline activation.** Pre-warm 15 mL 0.5 M Na₂CO₃ +
   0.5 M NaOH pH 11.0 buffer to 60 °C in a water bath. Transfer gel into
   the warm buffer. Rotate end-over-end in a 60 °C incubator, 2 h.
   Exact time matters — over-activation produces a brittle gel.
6. **Phase C — quench + wash.** Transfer gel into 20 mL ice-cold 0.1 M HCl
   on the sintered funnel (drops local pH to ~4, terminates phosphoramide
   formation). Wait 2 min. Then wash: 10 × 20 mL water (drain between each);
   3 × 20 mL 10 mM EDTA pH 8 (strips phosphate-chelated metals); 5 × 20 mL
   water. Drain damp.
7. Gel is now phosphate-diester crosslinked. Store in 20 % ethanol / 0.1 M
   NaCl pH 7 at 4 °C. Use within 3 months.

#### Quality control

- Bulk phosphorus (ICP-OES after acid digestion of a dried 10 mg aliquot):
  1.5–3.5 mmol P / g dry matrix (≈ 5–10 % crosslink conversion).
- FTIR (ATR, dried powder): P=O stretch at 1230 cm⁻¹ and P-O-C ester band
  at 990 cm⁻¹ both present.
- Equilibrium swelling (water, 24 h): ratio drops 20–40 % vs. uncrosslinked
  control from the same lot.
- Storage modulus (DMA, 1 Hz, swollen bead): 2–5× increase vs. control.

#### Troubleshooting

- **Skin-core structure (brittle outer shell, soft centre):** bead radius
  too large for the Phase B time. STMP is homogeneous for `d50/2 < 500 µm`
  (Thiele modulus < 1). If your L1 output shows d50/2 > 500 µm, either
  reduce bead size (raise rpm or add surfactant in L1) or shorten Phase B
  to 60 min and compensate with a second STMP cycle.
- **Gel melts / liquefies during Phase B:** temperature exceeded 80 °C
  or pH exceeded 12. Agarose double-helix is hydrolysed above these
  thresholds. Re-run with thermometer-verified 60 °C and a freshly
  pH-calibrated buffer. Do not exceed 70 °C without a pilot-scale
  validation of your specific matrix.
- **No measurable crosslinking (FTIR P-bands absent):** Phase A was
  too long (≥ 2 h) and the STMP reacted at pH 7 before transfer; or
  Phase B buffer drifted below pH 10. Freshly prepare both buffers on
  the day of use and keep Phase A at or below 30 min.
- **Downstream IMAC column loses metal:** unwashed phosphate residues
  chelate the loading ion. Extend the EDTA wash to 5 × 20 mL before
  metal charging.
- **The M1 UI shows `STMP homogeneity window exceeded` warning after
  the run:** d50/2 > 500 µm; see the skin-core entry above. The
  simulation result is still valid; the warning flags that the
  homogeneity assumption behind the kinetic fit may not hold.

#### References

Lim S, Seib PA (1993) *Cereal Chem.* 70:137. Kasemsuwan T, Jane J (1996)
*Cereal Chem.* 73:702. Lack S, Dulong V, Picton L, Le Cerf D, Condamine E
(2004) *Carbohydr. Res.* 339:2391. Seal BL (1996) *Biomaterials* 17:1869.
Salata GC, Kim JH, Chen C-H, McClements DJ (2015) *Int. J. Biol. Macromol.*
81:1009. Van Wazer JR (1958) *Phosphorus and its Compounds*, vol. I.

---

## J.2 Ligand Coupling

Attach a small-molecule ligand (dye, inhibitor, cofactor, peptide) to the
activated matrix. Select the protocol that matches the activation chemistry
you used in §J.1. The five common strategies covered here, plus click
chemistry for advanced users.

---

### J.2.1  Ligand Coupling — amine on CNBr- or CDI-activated matrix

**Purpose:** Couple a primary-amine-bearing ligand (e.g., peptide, amino
acid, aminated dye) to a CNBr-activated (§J.1.1) or CDI-activated (§J.1.5)
matrix. Forms an isourea (CNBr) or N-alkyl carbamate (CDI) linkage.

**Evidence tier:** VALIDATED

**EmulSim M2 key:** `amine_coupling`.

**Based on:** Cuatrecasas & Anfinsen (*Annu. Rev. Biochem.* 1971, 40:259);
Hermanson Ch. 16.2.

#### Safety

- GHS pictograms: none intrinsic beyond reagent-specific (the activated gel
  carries the CNBr / CDI inherited hazard until quenched).
- PPE: nitrile gloves, safety glasses, lab coat. Fume hood if the activation
  was CNBr-based.
- Waste: aqueous with trace activator — into the activator's waste stream.

#### What you need (per 10 mL activated wet gel)

| Reagent | Amount | Notes |
|---|---|---|
| Ligand (amine-bearing) | 2-20 μmol per mL wet gel target | In coupling buffer |
| Coupling buffer: 0.1 M NaHCO₃, 0.5 M NaCl, pH 8.3 | 20 mL | Freshly prepared |

#### Procedure

1. Dissolve the ligand in coupling buffer at 2-10 mg/mL. Adjust pH to 8.3
   with 1 M NaOH or 1 M HCl. Filter 0.22 μm if the ligand is clean enough.
2. Combine 10 mL activated wet gel (drained damp) with 20 mL ligand solution
   in a 50 mL tube.
3. Rotate end-over-end at 4 °C overnight (16 h) OR at RT for 2 h.
4. Drain on the sintered funnel, save the filtrate for coupling-yield
   analysis (measure ligand A₂₈₀ or your ligand's signal).
5. Wash with 3 × 30 mL coupling buffer.
6. Quench unreacted activated sites (§J.8.1). **This is critical.**

#### Quality control

- Coupling yield: compare ligand concentration in pre- and post-coupling
  solutions. 60-95 % incorporation is typical.
- Ligand density: (ligand added – ligand in wash) / gel volume.

#### Troubleshooting

- **< 30 % coupling:** ligand concentration too low, or pH drifted down.
  Verify pH at the start of coupling and at 1 h in.
- **Ligand activity lost:** check for non-specific adsorption by running a
  control coupling to an unreacted (non-activated) matrix; the ligand should
  wash off.

#### References

Cuatrecasas P, Anfinsen CB (1971) *Annu. Rev. Biochem.* 40:259. Hermanson (2013) Ch. 16.2.

---

### J.2.2  Ligand Coupling — amine / thiol / hydroxyl on epoxide matrix

**Purpose:** Couple ligands to an epoxide-activated matrix (from
epichlorohydrin §J.1.2 or BDDE §J.1.3). Epoxides react with amines (pH
9-11), thiols (pH 8-9), or hydroxyls (pH 11-12, slow).

**Evidence tier:** VALIDATED (amine, thiol); SEMI_QUANTITATIVE (hydroxyl)

**EmulSim M2 key:** `epoxide_coupling`.

**Based on:** Mateo et al. (*Enzyme Microb. Technol.* 2007, 39:274);
Hermanson Ch. 16.4.

#### Safety

- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: aqueous organic. Residual epoxide gets quenched in §J.8.2.

#### What you need (per 10 mL epoxide-activated wet gel)

| Reagent | Amount | Notes |
|---|---|---|
| Ligand | 5-50 μmol target | For amines, in 20 mL 0.5 M Na₂CO₃ pH 10. For thiols, in 20 mL 0.1 M Na phosphate pH 8.5 with 1 mM EDTA. |

#### Procedure (amine ligand)

1. Dissolve ligand in 0.5 M Na₂CO₃ pH 10 at 2-10 mg/mL.
2. Combine 10 mL activated wet gel with 20 mL ligand solution.
3. Rotate at 25 °C, 16-24 h (amine+epoxide is slow). For heat-sensitive
   ligands, extend to 48 h at RT.
4. Wash with 3 × 30 mL coupling buffer, 3 × 30 mL water.
5. Quench (§J.8.2).

#### Procedure (thiol ligand)

1. Dissolve thiolated ligand in 0.1 M Na phosphate pH 8.5, 1 mM EDTA.
   EDTA prevents oxidative disulfide formation.
2. Combine gel + ligand. Rotate at RT, 2-4 h (thiol+epoxide is fast).
3. Wash, quench (§J.8.2).

#### Quality control

- Amine coupling density: 5-25 μmol/mL wet gel typical.
- Thiol coupling density: 10-40 μmol/mL wet gel typical.

#### Troubleshooting

- **Very low amine coupling:** pH drifted. Add 1 M Na₂CO₃ to restore pH 10.
- **Very low thiol coupling:** thiols oxidised (disulfide). Include TCEP or
  freshly reduce the ligand with DTT then desalt before adding.

#### References

Mateo et al. (2007) *Enzyme Microb. Technol.* 39:274. Hermanson (2013) Ch. 16.4.

---

### J.2.3  Ligand Coupling — thiol on DVS / vinyl sulfone matrix

**Purpose:** Couple thiol-bearing ligands (cysteine-containing peptides,
thiolated probes) to a DVS-activated matrix (§J.1.4). Michael addition, fast
and clean.

**Evidence tier:** VALIDATED

**EmulSim M2 key:** `vinyl_sulfone_coupling`.

**Based on:** Morpurgo et al. (*Bioconjug. Chem.* 1996, 7:363); Hermanson Ch.
16.5.

#### Safety

- Active gel still carries unreacted vinyl sulfone — a sensitiser. Nitrile
  gloves mandatory.
- Waste: quench post-reaction (§J.8.3) before pooling.

#### What you need (per 10 mL DVS-activated wet gel)

| Reagent | Amount | Notes |
|---|---|---|
| Thiol ligand | 5-30 μmol | In 20 mL 0.1 M Na phosphate pH 8.0, 1 mM EDTA |
| TCEP | 0.5 mM | To keep ligand reduced |

#### Procedure

1. Dissolve thiol ligand. Add TCEP to 0.5 mM. Verify no precipitate.
2. Combine gel + ligand. Rotate at RT, 1-3 h. Coupling is often > 90 %
   complete in 30 min.
3. Wash 3 × 30 mL phosphate buffer.
4. Quench (§J.8.3).

#### Quality control

- Thiol coupling density: 15-40 μmol/mL wet gel.
- Free thiol on the product (should be zero): add Ellman's reagent (DTNB) to
  a 10 μL slurry in 1 mL of 0.1 M Na phosphate pH 8; measure A₄₁₂. A₄₁₂
  > 0.05 indicates unreacted thiol on surface (oxidise with iodoacetamide).

#### Troubleshooting

- **Low coupling:** DVS-activated gel aged before coupling — vinyl sulfone
  hydrolysed. Must couple within 30 min of activation.

#### References

Morpurgo M, Veronese FM, Kachensky D, Harris JM (1996) *Bioconjug. Chem.*
7:363. Hermanson (2013) Ch. 16.5.

---

### J.2.4  Ligand Coupling — hydrazone on aldehyde matrix

**Purpose:** Couple hydrazide-functional ligands to an aldehyde-bearing
matrix (from periodate oxidation of agarose-bound diols — see §J.6.3 for
matrix oxidation). The hydrazone bond is pH-labile and reversible unless
reduced with NaCNBH₃.

**Evidence tier:** SEMI_QUANTITATIVE

**EmulSim M2 key:** `hydrazone_coupling`.

**Based on:** O'Shannessy et al. (*J. Immunol. Methods* 1984, 75:11);
Hermanson Ch. 19.2.

#### Safety

- Sodium cyanoborohydride (NaCNBH₃) is highly toxic (H300+H311+H331), and
  releases HCN on contact with acid. Fume hood mandatory.
- PPE: double nitrile gloves, face shield, chemical apron.
- Waste: alkaline aqueous; add 1 M NaOH + sodium hypochlorite before
  disposal to destroy cyanide residue.

#### What you need (per 10 mL aldehyde-activated wet gel)

| Reagent | Amount | CAS |
|---|---|---|
| Hydrazide ligand | 5-30 μmol | — |
| Coupling buffer: 0.1 M Na phosphate pH 6.0 | 20 mL | — |
| NaCNBH₃ | 50 mg | 25895-60-7 |

#### Procedure

1. Dissolve ligand in coupling buffer pH 6.0 (hydrazone forms best at pH
   5-6).
2. Combine 10 mL aldehyde-activated wet gel with 20 mL ligand solution.
3. In the fume hood, add 50 mg NaCNBH₃. Cap.
4. Rotate at RT, 4-24 h.
5. Wash with 3 × 30 mL PBS.
6. Quench unreacted aldehydes with 100 mM ethanolamine or 100 mM glycine
   pH 7.4 + 50 mg NaCNBH₃, 2 h RT (§J.8.5).

#### Quality control

- Coupling yield typically 50-90 %.
- Acid stability: treat 100 μL slurry with 0.1 M HCl 10 min; the hydrazone
  should survive (< 5 % leach) because NaCNBH₃ has reduced it to a stable
  hydrazine.

#### Troubleshooting

- **Ligand leaches over time:** reduction step was skipped or NaCNBH₃
  exhausted (aged reagent). Repeat the reduction with fresh NaCNBH₃.

#### References

O'Shannessy DJ, Dobersen MJ, Quarles RH (1984) *J. Immunol. Methods* 75:11.
Hermanson (2013) Ch. 19.2.

---

### J.2.5  Ligand Coupling — click chemistry (CuAAC or SPAAC)

**Purpose:** Attach an alkyne- or azide-functional ligand to an
azide- or alkyne-functional matrix via 1,3-dipolar cycloaddition. Two
variants: CuAAC (copper-catalysed, fast, needs Cu removal) and SPAAC
(strain-promoted, Cu-free, slower but biocompatible).

**Evidence tier:** VALIDATED (CuAAC on polymeric supports);
SEMI_QUANTITATIVE (SPAAC on gel matrices).

**EmulSim M2 key:** `click_cuaac`, `click_spaac` (not yet wired).

**Based on:** Kolb, Finn, Sharpless (*Angew. Chem.* 2001, 40:2004); Meldal &
Tornøe (*Chem. Rev.* 2008, 108:2952); Hermanson Ch. 17.

#### Safety (CuAAC)

- Copper sulfate: H302 (harmful if swallowed), H315, H319, H410.
- Sodium ascorbate: low hazard.
- THPTA / BTTAA (Cu ligand): exclamation mark.
- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: Cu-containing aqueous must go to the heavy-metal waste stream.

#### Safety (SPAAC)

- DBCO reagents: H315, H319, H317. Some are expensive but lab-safe at
  μmolar scale.
- PPE: standard lab PPE.

#### What you need (per 10 mL azide-activated wet gel, CuAAC)

| Reagent | Amount | CAS |
|---|---|---|
| Alkyne ligand | 5-30 μmol | — |
| CuSO₄·5H₂O | 5 mM final | 7758-99-8 |
| Sodium ascorbate | 25 mM final | 134-03-2 |
| THPTA (copper stabiliser) | 25 mM final | 760952-88-3 |
| Coupling buffer: 0.1 M Tris pH 7.4 | 20 mL | — |

#### Procedure (CuAAC)

1. Mix pre-formed Cu(II)-THPTA complex: in a small tube add CuSO₄ and THPTA
   to final concentrations 5 mM and 25 mM in water. Mix 1 min (blue).
2. In the coupling tube, add 10 mL wet gel, 20 mL Tris buffer, the alkyne
   ligand, and the Cu-THPTA premix.
3. Add sodium ascorbate to 25 mM final (this reduces Cu(II) → Cu(I) in situ).
4. Cap, displace headspace with nitrogen or argon (click is O₂-sensitive
   because ascorbate reoxidises).
5. Rotate at RT, 2-4 h, in the dark.
6. Wash with 3 × 30 mL Tris-EDTA (10 mM EDTA) to scavenge Cu, then 5 × 30 mL
   water, then 3 × 30 mL storage buffer.

#### Procedure (SPAAC, Cu-free)

1. Dissolve DBCO-ligand in PBS pH 7.4 at 1-5 mM. DBCO is often DMSO-soluble
   only — keep DMSO ≤ 10 % v/v in the reaction.
2. Combine with 10 mL azide-activated gel + PBS to 30 mL total.
3. Rotate at RT, 4-24 h.
4. Wash with 3 × 30 mL PBS, 3 × 30 mL PBS + 0.1 % Tween-20 (removes
   non-specifically bound DBCO excess), 3 × 30 mL PBS.

#### Quality control

- Fluorescence (if fluorescent ligand) on the gel; compare to free-ligand
  spectrum.
- Coupling yield from solution depletion: 50-95 % (CuAAC), 30-80 % (SPAAC).

#### Troubleshooting

- **CuAAC dead at 30 min:** oxygen ingress. Re-seal, re-add 10 mM ascorbate.
- **SPAAC slow:** DBCO aged / hydrolysed. Use within 1 week of reconstitution.

#### References

Kolb HC, Finn MG, Sharpless KB (2001) *Angew. Chem. Int. Ed.* 40:2004. Meldal M,
Tornøe CW (2008) *Chem. Rev.* 108:2952. Hermanson (2013) Ch. 17.

---

## J.3 Protein Coupling

Immobilise an antibody, enzyme, lectin, or other protein onto an activated
matrix. Protein coupling differs from ligand coupling in three ways:

1. Proteins are large (50-200 kDa). Steric access to internal lysines or
   cysteines is limited; surface residues dominate.
2. Protein activity is sensitive to pH, temperature, and side reactions.
   Prefer mild chemistries.
3. Orientation matters. Random coupling scrambles active-site access; use
   oriented strategies when possible.

Pretreat the protein as in §J.6 before coupling.

---

### J.3.1  Protein Coupling — NHS-ester (EDC/NHS activated matrix)

**Purpose:** Couple amine-bearing proteins (essentially every protein — via
surface lysines and N-terminus) to an NHS-ester-activated matrix. Most
common, most flexible, non-oriented.

**Evidence tier:** VALIDATED

**EmulSim M2 key:** `nhs_protein`.

**Based on:** Staros (*Biochemistry* 1982, 21:3950); Hermanson Ch. 16.3.

#### Safety

- EDC: H302, H315, H319.
- NHS (N-hydroxysuccinimide): H315.
- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: aqueous organic.

#### What you need (per 10 mL wet gel)

| Reagent | Amount | CAS |
|---|---|---|
| Carboxylated matrix (or NHS-activated matrix from vendor) | 10 mL | — |
| EDC·HCl | 50 mg (for in-situ activation) | 25952-53-8 |
| Sulfo-NHS | 15 mg | 106627-54-7 |
| Activation buffer: 0.1 M MES pH 5.0 | 20 mL | — |
| Coupling buffer: 50 mM sodium phosphate pH 7.4, 150 mM NaCl | 20 mL | — |
| Protein | 1-5 mg / mL gel | Pretreated per §J.6 |

#### Procedure

1. (If using pre-activated NHS-agarose from vendor: skip to step 5.)
2. In-situ activation: combine 10 mL carboxylated wet gel + 20 mL MES pH 5
   + 50 mg EDC + 15 mg sulfo-NHS. Rotate RT, 15 min.
3. Rapidly wash with 3 × 30 mL ice-cold MES pH 5 on the funnel. Drain damp.
4. Immediately proceed (NHS-ester half-life ≈ 30 min at pH 7, 10 min at pH 8).
5. Mix drained gel with protein in coupling buffer pH 7.4 (20 mL total).
6. Rotate at 4 °C, 2-4 h (OR RT 1 h for activity-robust proteins).
7. Wash with 3 × 30 mL coupling buffer, 3 × 30 mL PBS + 0.5 M NaCl
   (to remove non-covalently bound protein), 3 × 30 mL PBS.
8. Quench (§J.8.1).

#### Quality control

- Bradford / A₂₈₀ on pre- and post-coupling solutions. 60-90 % of added
  protein typically couples.
- Activity assay on the immobilised protein (e.g., enzyme turnover in a
  packed micro-column). Retained activity is typically 30-80 % of soluble
  protein.

#### Troubleshooting

- **Low coupling:** NHS-ester expired. Re-activate or use fresh matrix.
- **Low activity:** active site lysines got coupled. Try oriented coupling
  (§J.3.4 or §J.3.5) instead.
- **Non-specific sticking:** ionic interactions with matrix. Run a
  higher-salt wash (1 M NaCl).

#### References

Staros JV (1982) *Biochemistry* 21:3950. Hermanson (2013) Ch. 16.3.

---

### J.3.2  Protein Coupling — glutaraldehyde (GA) two-step method

**Purpose:** Immobilise enzymes (esp. proteases, lipases) by first
activating an amine-matrix with glutaraldehyde, then exposing to protein.
Creates an imine / reduced-amine network; robust but non-oriented.

**Evidence tier:** VALIDATED (especially for enzyme immobilisation)

**EmulSim M2 key:** `glutaraldehyde`.

**Based on:** Monsan (*J. Mol. Catal.* 1978, 3:371); Migneault et al.
(*BioTechniques* 2004, 37:790); Hermanson Ch. 16.1.

#### Safety

- Glutaraldehyde: H301, H330, H314, H317, H334 (respiratory sensitiser).
  Fume hood mandatory.
- PPE: double nitrile gloves, face shield, chemical apron, lab coat.
- Waste: neutralise with glycine or ethanolamine before disposal into aqueous
  organic waste.

#### What you need (per 10 mL amine-matrix wet gel)

| Reagent | Amount | CAS |
|---|---|---|
| Glutaraldehyde 25 % aq. | 2 mL | 111-30-8 |
| 0.1 M Na phosphate pH 7.0 | 20 mL | — |
| Sodium cyanoborohydride (NaCNBH₃) | 50 mg | 25895-60-7 |
| Protein | 1-5 mg / mL gel | — |

#### Procedure

1. PPE on. Fume hood mandatory for glutaraldehyde.
2. Wash 10 mL amine-matrix with 3 × 30 mL 0.1 M Na phosphate pH 7.
3. Activate: combine gel + 20 mL buffer + 2 mL 25 % glutaraldehyde (final
   ≈ 2.3 % GA). Rotate RT, 1-2 h. Matrix turns yellow.
4. Wash rapidly with 5 × 30 mL phosphate buffer until filtrate is clear.
5. Immediately add protein in 20 mL phosphate buffer. Rotate 4 °C, 16 h.
6. Add 50 mg NaCNBH₃ (in the fume hood) to reduce the Schiff base to a
   stable secondary amine. Rotate 4 °C, 2 h.
7. Wash with 3 × 30 mL PBS + 0.5 M NaCl, 3 × 30 mL PBS.
8. Quench (§J.8.5).

#### Quality control

- Activity: for enzymes, measure before and after coupling. Glutaraldehyde
  often gives high activity (40-80 %) for proteases and lipases.

#### Troubleshooting

- **Gel crosslinks excessively:** GA concentration too high. Reduce to 1 %.
- **Low activity:** reduce coupling time to 4 h or reduce GA concentration.

#### References

Monsan P (1978) *J. Mol. Catal.* 3:371. Migneault I et al. (2004)
*BioTechniques* 37:790. Hermanson (2013) Ch. 16.1.

---

### J.3.3  Protein Coupling — Protein A + DMP crosslinking (oriented IgG)

**Purpose:** Immobilise an antibody in defined orientation by first binding
it via Fc region to immobilised Protein A / G, then crosslinking with
dimethyl pimelimidate (DMP). Fab domains point outward, antigen-binding
retained.

**Evidence tier:** VALIDATED

**EmulSim M2 key:** `proteinA_crosslink`.

**Based on:** Schneider et al. (*J. Biol. Chem.* 1982, 257:10766); Hermanson Ch. 20.2.

#### Safety

- DMP: H315, H319. Low hazard.
- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: aqueous.

#### What you need (per 1 mL Protein-A-agarose)

| Reagent | Amount | CAS |
|---|---|---|
| Protein A-agarose | 1 mL | (commercial) |
| IgG | 1-10 mg | — |
| PBS pH 7.4 | 10 mL | — |
| Triethanolamine buffer: 0.2 M, pH 8.2 | 10 mL | 102-71-6 |
| DMP·2HCl (dimethyl pimelimidate dihydrochloride) | 20 mg | 58537-94-3 |
| Quench buffer: 50 mM ethanolamine pH 8.2 | 10 mL | 141-43-5 |
| Storage buffer: PBS + 0.05 % NaN₃ | 20 mL | — |

#### Procedure

1. Equilibrate 1 mL Protein A-agarose with 5 × 2 mL PBS.
2. Incubate gel with 1-10 mg IgG in 5 mL PBS, rotate RT 30-60 min. Affinity
   capture.
3. Wash 3 × 10 mL PBS to remove unbound IgG.
4. Exchange to 10 × 2 mL 0.2 M triethanolamine pH 8.2. Drain.
5. Dissolve 20 mg DMP in 2 mL triethanolamine buffer (fresh — DMP hydrolyses
   in minutes). Quickly add to gel.
6. Rotate RT, 30 min.
7. Quench: drain, add 10 mL 50 mM ethanolamine pH 8.2, rotate 5 min.
8. Wash with 3 × 10 mL PBS, 3 × 10 mL 100 mM glycine pH 2.5 (strips
   uncrosslinked IgG), 3 × 10 mL PBS.
9. Store in PBS + 0.05 % sodium azide at 4 °C.

#### Quality control

- Binding capacity: test with known antigen; should be 60-80 % of
  equivalent solution-phase IgG binding.
- Leaching: incubate gel in 100 mM glycine pH 2.5 for 10 min. IgG leach
  in the supernatant should be < 10 % of immobilised.

#### Troubleshooting

- **High IgG leach on low-pH wash:** crosslinking incomplete. Repeat steps
  5-7 with fresh DMP.
- **Low antigen binding:** DMP over-crosslinked the Fab region. Reduce DMP
  to 10 mg or reduce time to 15 min.

#### References

Schneider C et al. (1982) *J. Biol. Chem.* 257:10766. Hermanson (2013) Ch. 20.2.

---

### J.3.4  Protein Coupling — sortase A (site-specific, LPXTG tag)

**Purpose:** Oriented, site-specific coupling using Sortase A enzyme. The
protein of interest must carry a C-terminal LPXTG tag (genetically engineered).
The matrix must present N-terminal polyglycine (Gly₃ or Gly₅). Sortase
transpeptidates: LPXTG cleaved, LPXT-Gly_n product covalently attached.

**Evidence tier:** VALIDATED

**EmulSim M2 key:** `sortase` (not yet wired).

**Based on:** Popp & Ploegh (*Angew. Chem. Int. Ed.* 2011, 50:5024);
Mao et al. (*J. Am. Chem. Soc.* 2004, 126:2670).

#### Safety

- Low-hazard reagents throughout. Standard lab PPE.
- Sortase A itself has no known hazards.
- Waste: aqueous.

#### What you need (per 1 mL Gly₃-agarose)

| Reagent | Amount | Notes |
|---|---|---|
| Gly₃-agarose matrix | 1 mL | See §J.2 for coupling a Gly-Gly-Gly peptide to an activated matrix |
| LPXTG-tagged target protein | 50 μmol | Pretreated per §J.6 |
| Sortase A enzyme | 10 μM | Recombinant; Addgene or commercial |
| Sortase buffer: 50 mM Tris-HCl pH 7.5, 150 mM NaCl, 10 mM CaCl₂ | 5 mL | Ca²⁺ is essential |

#### Procedure

1. Equilibrate 1 mL Gly₃-agarose with 3 × 5 mL sortase buffer.
2. Combine gel + 5 mL sortase buffer + 50 μmol LPXTG-protein + 10 μM sortase A.
3. Rotate RT, 2-6 h (monitor reaction by SDS-PAGE; the LPXTG-protein band
   should shift slightly and the soluble supernatant should deplete).
4. Wash with 3 × 10 mL sortase buffer.
5. Wash with 3 × 10 mL sortase buffer WITHOUT CaCl₂ + 5 mM EGTA to
   inactivate residual sortase.
6. Wash 3 × 10 mL storage buffer (PBS + 0.05 % NaN₃). Store 4 °C.

#### Quality control

- Incorporation: 50-90 % of added LPXTG-protein couples.
- Activity: often > 80 % retained because the C-terminus is non-functional
  for most target proteins.

#### Troubleshooting

- **No coupling:** target protein LPXTG tag not accessible (buried).
  Verify by running sortase reaction in solution with a Gly-peptide — should
  produce the cleavage product.
- **Slow coupling:** low sortase activity. Titrate sortase from 1 to 50 μM.

#### References

Popp MW, Ploegh HL (2011) *Angew. Chem. Int. Ed.* 50:5024. Mao H et al.
(2004) *J. Am. Chem. Soc.* 126:2670.

---

### J.3.5  Protein Coupling — oriented IgG via Fc-glycan oxidation + hydrazide

**Purpose:** Couple an IgG in defined orientation by oxidising the carbohydrate
diols in the Fc region (away from the antigen-binding site) to aldehydes,
then coupling to a hydrazide-agarose. Fully oriented, preserves antigen
binding.

**Evidence tier:** VALIDATED

**EmulSim M2 key:** `fc_hydrazone`.

**Based on:** O'Shannessy et al. (*J. Immunol. Methods* 1984, 75:11);
Hermanson Ch. 19.2.

#### Safety

- NaIO₄: H272 (oxidiser), H315, H319. Keep dry.
- NaCNBH₃: H300+H311+H331. Fume hood, double gloves.
- PPE: double nitrile, face shield, chemical apron.
- Waste: alkaline + bleach treatment for NaCNBH₃-containing streams (see §J.2.4).

#### What you need (per 5 mg IgG)

| Reagent | Amount | CAS |
|---|---|---|
| IgG (pretreated per §J.6.4) | 5 mg, 1 mg/mL in PBS | — |
| NaIO₄ | 10 mM final | 7790-28-5 |
| Coupling buffer: 0.1 M Na phosphate pH 6.0 | 20 mL | — |
| Hydrazide-agarose (e.g., via hydrazide §J.2.4) | 1 mL | — |
| NaCNBH₃ | 50 mg | 25895-60-7 |

#### Procedure

1. Dissolve NaIO₄ in ice-cold PBS to 100 mM. Add to 5 mg IgG (5 mL at 1 mg/mL)
   to a final 10 mM NaIO₄. **On ice, in the dark, 30 min.** Over-oxidation
   destroys activity.
2. Quench excess NaIO₄ with ethylene glycol (100 mM final, 15 min).
3. Buffer-exchange to pH 6.0 phosphate buffer via PD-10 column (§J.6.4).
4. Add oxidised IgG to 1 mL hydrazide-agarose + 50 mg NaCNBH₃ in the fume hood.
5. Rotate RT, 16-24 h.
6. Wash 3 × 10 mL PBS + 0.5 M NaCl, 3 × 10 mL PBS.
7. Quench residual aldehydes with 100 mM glycine pH 7.4 + 20 mg NaCNBH₃, 2 h RT.
8. Final wash, store.

#### Quality control

- Coupling yield: 50-80 % of oxidised IgG incorporates.
- Antigen binding: typically 60-90 % retained relative to non-oriented
  NHS-coupled IgG.

#### Troubleshooting

- **Low antigen binding:** NaIO₄ concentration or time too high. Reduce to
  5 mM / 15 min.

#### References

O'Shannessy DJ et al. (1984) *J. Immunol. Methods* 75:11. Hermanson (2013) Ch. 19.2.

---

## J.4 Spacer Arm

Spacer arms are short bifunctional molecules inserted between the matrix
and the ligand. They reduce steric hindrance so that larger targets (e.g.,
antibodies binding a small immobilised hapten) can approach the ligand.

Pick a spacer based on (a) length, (b) hydrophilicity (longer alkyl chains
increase non-specific binding), and (c) the functional groups you need at
each end. Spacers are typically coupled to an activated matrix first (§J.1
→ one end of spacer), then the free end of the spacer is activated
(usually by EDC/NHS for a carboxyl end) and reacted with the ligand (§J.2/J.3).

---

### J.4.1  Spacer Arm — Ethylenediamine (EDA, 3 Å diamine)

**Purpose:** Shortest diamine spacer. Both ends primary amines. Use when a
small ligand needs only a slight offset from the matrix surface.

**Evidence tier:** VALIDATED

**Based on:** Cuatrecasas (*J. Biol. Chem.* 1970, 245:3059); Hermanson Ch. 5.2.

#### Safety

- EDA: H302, H312, H314, H317, H334. Strong base, strong sensitiser.
- PPE: nitrile gloves, safety glasses, face shield, lab coat.
- Fume hood for weighing.
- Waste: neutralise before disposal (aqueous organic).

#### What you need (per 10 mL NHS-activated or CNBr-activated wet gel)

| Reagent | Amount | CAS |
|---|---|---|
| Ethylenediamine | 1 mL (neat) or 10 mL of 1 M in coupling buffer | 107-15-3 |
| Coupling buffer: 0.1 M NaHCO₃, 0.5 M NaCl, pH 8.3 | 20 mL | — |

#### Procedure

1. Prepare 1 M EDA in coupling buffer; adjust pH to 8.3 with HCl. **EDA is
   a strong base — add to buffer, not vice versa.**
2. Couple to activated matrix per §J.2.1. Use 20-fold excess EDA over
   activated sites so that mostly one amine of EDA reacts and the other
   stays free.
3. Wash, quench (§J.8.1 for NHS). The gel now presents terminal primary
   amines as the spacer free end.

#### QC

- Ninhydrin test (10 μL gel, add 2 mL 0.1 % ninhydrin in ethanol, heat 90 °C
  5 min): strong purple = abundant primary amines.

#### References

Cuatrecasas P (1970) *J. Biol. Chem.* 245:3059. Hermanson (2013) Ch. 5.2.

---

### J.4.2  Spacer Arm — 1,6-Diaminohexane (11 Å diamine)

**Purpose:** Longer diamine, same chemistry as EDA. Gives 11 Å of reach;
useful for ligands binding into a pocket (e.g., immobilised biotin reaching
into streptavidin's binding site).

**Evidence tier:** VALIDATED

**Based on:** Hermanson Ch. 5.2.

#### Safety

- 1,6-Diaminohexane: H302, H312, H314, H317. Similar hazard profile to EDA
  but less volatile.
- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: aqueous organic, neutralise first.

#### What you need (per 10 mL NHS-activated wet gel)

| Reagent | Amount | CAS |
|---|---|---|
| 1,6-Diaminohexane | 1 g (dissolve in 20 mL coupling buffer) | 124-09-4 |
| Coupling buffer: 0.1 M NaHCO₃, 0.5 M NaCl, pH 8.3 | 20 mL | — |

#### Procedure

Same as §J.4.1 but substitute 1,6-diaminohexane. 1 g in 20 mL = 0.43 M.

#### QC / Refs

Ninhydrin as above. Hermanson (2013) Ch. 5.2.

---

### J.4.3  Spacer Arm — 6-Aminohexanoic acid (9 Å amine-carboxyl)

**Purpose:** Heterobifunctional 9 Å spacer: amine one end, carboxyl the
other. Use when the matrix is amine-activated and the ligand to be coupled
is a carboxyl, or vice versa, via EDC/NHS chemistry.

**Evidence tier:** VALIDATED

**Based on:** Hermanson Ch. 5.3.

#### Safety

- 6-Aminohexanoic acid: low hazard (H315, H319).
- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: aqueous, routine.

#### What you need (per 10 mL NHS-activated amine-matrix)

| Reagent | Amount | CAS |
|---|---|---|
| 6-Aminohexanoic acid | 500 mg | 60-32-2 |
| Coupling buffer: 0.1 M NaHCO₃, 0.5 M NaCl, pH 8.3 | 20 mL | — |

#### Procedure

1. Dissolve 500 mg 6-aminohexanoic acid in 20 mL coupling buffer. pH 8.3.
2. Couple to NHS-activated matrix per §J.2.1. The amine of the spacer couples;
   the carboxyl is left free for later ligand coupling.
3. Wash and store.

#### Carboxyl activation of the spacer free end (when coupling a protein)

- Per §J.3.1 — in-situ EDC/NHS activation, then add protein.

#### Refs

Hermanson (2013) Ch. 5.3.

---

### J.4.4  Spacer Arm — DADPA (3,3'-diamino-N-methyldipropylamine, 11 Å triamine)

**Purpose:** 11 Å spacer with three nitrogen functional handles. Useful for
multi-point attachment or for a heterotrifunctional linker.

**Evidence tier:** SEMI_QUANTITATIVE

**Based on:** Hermanson Ch. 5.2.

#### Safety

- DADPA: H302, H314, H317. Strong base.
- PPE: nitrile gloves, face shield, lab coat.
- Waste: aqueous organic, neutralise.

#### What you need

| Reagent | Amount | CAS |
|---|---|---|
| DADPA | 500 mg | 105-83-9 |
| Coupling buffer pH 8.3 | 20 mL | — |

#### Procedure

Same as §J.4.1. 500 mg in 20 mL = 0.15 M; 20-fold over activated sites.

#### Refs

Hermanson (2013) Ch. 5.2.

---

### J.4.5  Spacer Arm — PEG bis-amine (~30-350 Å, tunable)

**Purpose:** Polyethylene glycol diamine of defined MW. Highly hydrophilic
(prevents non-specific binding). Lengths from 500 Da (~3 nm) to 5000 Da
(~30 nm).

**Evidence tier:** VALIDATED

**Based on:** Harris JM (ed.) *Poly(ethylene glycol) Chemistry: Biotechnical
and Biomedical Applications* (1992); Hermanson Ch. 18.

#### Safety

- PEG bis-amine (3400 Da): very low hazard. H303 (may be harmful if swallowed).
- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: aqueous, routine.

#### What you need (per 10 mL NHS-activated wet gel)

| Reagent | Amount | CAS |
|---|---|---|
| PEG-bisamine 3.4 kDa | 500 mg | 24991-53-5 (for 3.4 kDa) |
| Coupling buffer pH 8.3 | 20 mL | — |

#### Procedure

1. Dissolve 500 mg PEG-bisamine in 20 mL coupling buffer. Expect a slightly
   viscous solution.
2. Couple to NHS-activated matrix per §J.2.1. Overnight at 4 °C to drive
   the reaction to completion on both ends of any PEG chain that orients
   away from the matrix.
3. Wash with 5 × 30 mL coupling buffer + 0.5 M NaCl (removes non-covalently
   bound PEG), 3 × 30 mL PBS.
4. Quench.

#### QC

- Ninhydrin on end-points: should be strongly positive (many amine termini).
- MW-discrimination: 3.4 kDa PEG gives ~30 nm of reach; confirm by titrating
  with a fluorescent amine probe and measuring accessible density.

#### Refs

Harris JM (ed.) (1992). Hermanson (2013) Ch. 18.

---

### J.4.6  Spacer Arm — Jeffamine (polyether diamine, 600-2000 Da)

**Purpose:** Jeffamine is a family of polyether diamines (commercial,
Huntsman). Lengths 600-2000 Da, more hydrophobic than PEG-bisamine but
cheaper. Use for enzyme immobilisation where modest hydrophobicity is OK.

**Evidence tier:** SEMI_QUANTITATIVE

**Based on:** Product literature (Huntsman); Hermanson Ch. 18.

#### Safety

- Jeffamine: H302, H315, H318.
- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: aqueous organic.

#### Procedure

Same pattern as §J.4.5 with Jeffamine (D-2000 or D-400 depending on desired
length).

#### Refs

Hermanson (2013) Ch. 18.

---

## J.5 Metal Charging

IMAC (Immobilised Metal Affinity Chromatography). The matrix has a chelator
group (NTA, IDA, or TED). Charging adds a specific metal ion onto the
chelator. The loaded metal then binds a target with the complementary
metal-coordinating group — most commonly a His-tag on a recombinant protein.

Load the matrix the day of purification. Metals leach during elution; if
the matrix has been used, strip and reload (see final section).

---

### J.5.1  Metal Charging — Ni²⁺ on NTA

**Purpose:** The workhorse for His-tag purification. Ni²⁺-NTA has a Kd of
~10⁻⁶ M for a His₆-tag.

**Evidence tier:** VALIDATED

**EmulSim M2 key:** `ni_nta`.

**Based on:** Hochuli et al. (*Biotechnology* 1988, 6:1321); Porath et al.
(*Nature* 1975, 258:598).

#### Safety

- NiSO₄·6H₂O: H302, H315, H317, H334, H341, H350, H360, H372, H410. CMR
  reagent.
- PPE: double nitrile gloves, safety glasses, lab coat.
- Waste: Ni-containing aqueous to heavy-metal waste stream. Solid Ni salts:
  do NOT place in general chemical waste.

#### What you need (per 10 mL NTA-agarose bed)

| Reagent | Amount | CAS |
|---|---|---|
| NiSO₄·6H₂O 100 mM | 20 mL | 10101-97-0 |
| Water | ≥ 200 mL | — |
| Storage buffer: 20 % ethanol | 20 mL | 64-17-5 |

#### Procedure

1. Wash NTA-agarose with 5 CV water to remove storage buffer.
2. Load 2 CV of 100 mM NiSO₄ at ≤ 1 CV/min (gravity is fine). Matrix turns
   light blue-green.
3. Wash 10 CV water to remove unbound Ni²⁺.
4. Wash 5 CV binding buffer (usually 50 mM Na phosphate pH 7.5, 300 mM NaCl,
   10 mM imidazole).
5. Matrix is ready for His-tag loading. Do not store charged; re-load the
   day of use.

#### QC

- Colour: light blue-green. If pale, recharge.
- Capacity: 30-50 mg His-tagged protein per mL matrix (specification;
  measure by loading a known protein in excess).

#### Troubleshooting

- **Weak protein binding:** low Ni loading or leached. Strip and recharge.
- **Non-specific protein binding:** histidine-rich non-target proteins
  co-purify. Increase imidazole in binding/wash buffers to 20-40 mM.

#### References

Hochuli E, Döbeli H, Schacher A (1988) *Biotechnology* 6:1321. Porath J et
al. (1975) *Nature* 258:598.

---

### J.5.2  Metal Charging — Co²⁺ on NTA (or CMA)

**Purpose:** Higher-selectivity alternative to Ni²⁺. Binds His-tag less
tightly (Kd ~10⁻⁵ M) but co-purifies fewer endogenous His-rich contaminants.
Use when purity > yield.

**Evidence tier:** VALIDATED

**EmulSim M2 key:** `co_nta`.

#### Safety

- CoCl₂·6H₂O: H302, H317, H334, H341, H350, H360, H410. CMR reagent.
- PPE: double nitrile, safety glasses, lab coat.
- Waste: Co-containing aqueous to heavy-metal waste.

#### Procedure

As §J.5.1, substitute 100 mM CoCl₂·6H₂O. Matrix turns pink-purple.

#### Refs

Hochuli (1988); Porath (1975).

---

### J.5.3  Metal Charging — Cu²⁺ on IDA

**Purpose:** Broad specificity — binds proteins with surface histidines,
cysteines, or acidic patches. Use for screening "generic" metal-binding
proteins.

**Evidence tier:** VALIDATED

**EmulSim M2 key:** `cu_ida`.

#### Safety

- CuSO₄·5H₂O: H302, H315, H319, H410.
- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: Cu-containing aqueous to heavy-metal waste.

#### Procedure

As §J.5.1, substitute 100 mM CuSO₄·5H₂O on an IDA-agarose (not NTA; Cu²⁺
works poorly on NTA). Matrix turns deep blue.

#### QC

- Strong blue colour. Capacity 20-30 mg/mL matrix for general
  metal-binding proteins.

#### Refs

Porath (1975). Hermanson Ch. 10.

---

### J.5.4  Metal Charging — Zn²⁺ on IDA

**Purpose:** Niche use: immobilised lectins, carbonic anhydrases, some
phosphopeptide enrichment (limited performance vs Fe³⁺).

**Evidence tier:** QUALITATIVE_TREND

**EmulSim M2 key:** `zn_ida` (not yet wired).

#### Safety

- ZnCl₂: H302, H314, H410.
- PPE: nitrile gloves, face shield, lab coat.
- Waste: Zn-containing aqueous to heavy-metal waste.

#### Procedure

As §J.5.1 with 100 mM ZnCl₂ in water (no organic buffer — ZnCl₂ hydrolyses).
Adjust water pH to 4.5 with HCl before charging.

---

### J.5.5  Metal Charging — Fe³⁺ on IDA (phosphopeptide enrichment)

**Purpose:** Phosphopeptide / phosphoprotein enrichment. Fe³⁺ coordinates to
the phosphate oxygens. Use after a tryptic digest for phosphoproteomics.

**Evidence tier:** VALIDATED (phosphopeptide enrichment)

**EmulSim M2 key:** `fe_ida`.

**Based on:** Andersson & Porath (*Anal. Biochem.* 1986, 154:250).

#### Safety

- FeCl₃·6H₂O: H290, H302, H315, H318.
- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: Fe-containing aqueous; most institutions accept in general aqueous
  waste at trace levels. Check local policy.

#### What you need (per 1 mL IDA-agarose)

| Reagent | Amount | CAS |
|---|---|---|
| FeCl₃·6H₂O | 100 mM in 10 mM HCl, 2 mL | 10025-77-1 |
| 0.1 M acetic acid pH 3.0 | 20 mL | 64-19-7 |

#### Procedure

1. Wash 1 mL IDA-agarose with 5 mL water, then 5 mL 0.1 M acetic acid pH 3.0.
2. Load 2 mL 100 mM FeCl₃ in 10 mM HCl. Matrix turns orange-brown.
3. Wash 10 mL 0.1 M acetic acid pH 3.0 to remove unbound Fe.
4. Equilibrate in binding buffer (typically 0.1 % TFA in 30 % acetonitrile
   for phosphopeptide enrichment).
5. Immediately load sample.

#### Refs

Andersson L, Porath J (1986) *Anal. Biochem.* 154:250.

---

### J.5.bonus Stripping a used IMAC column

When a charged IMAC column has been used (and leached metal), strip
completely before recharging:

1. Wash 5 CV water.
2. Strip with 5 CV 50 mM EDTA, pH 8.0 (chelates and removes all metal).
3. Wash 10 CV water.
4. Recharge per §J.5.1-5.5.

Never re-load metal onto a column that still has residual metal — you'll
co-precipitate oxides.

---

## J.6 Protein Pretreatment

Before coupling a protein to a matrix, pretreat to (a) make the coupling
chemistry effective, (b) preserve activity, and (c) remove contaminants
that compete with the coupling reaction.

---

### J.6.1  Protein Pretreatment — Disulfide reduction with DTT

**Purpose:** Reduce protein disulfide bonds to free thiols. Required for
thiol-maleimide coupling or DVS-thiol coupling (§J.2.3). Not recommended
for proteins whose tertiary structure depends on disulfides (most
antibodies, most secreted proteins).

**Evidence tier:** VALIDATED

#### Safety

- DTT (dithiothreitol): H302, H315, H319. Low hazard.
- Must be removed after reduction — DTT competes with thiols for maleimide.
- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: aqueous, routine.

#### What you need

| Reagent | Amount | CAS |
|---|---|---|
| DTT | 5-10 mM final | 3483-12-3 |
| Reducing buffer: 0.1 M Na phosphate pH 7.5, 1 mM EDTA | 10 mL | — |

#### Procedure

1. Dissolve protein (≥ 1 mg/mL) in reducing buffer.
2. Add DTT from 100 mM stock to 5-10 mM final.
3. Incubate 30 min at RT (or 37 °C for tough disulfides).
4. **Immediately** buffer-exchange to remove DTT (§J.6.4). DTT competes
   with target thiols in coupling.
5. Verify free-thiol generation with Ellman's reagent.

#### Troubleshooting

- **Protein precipitates:** disulfides held the structure. Switch to TCEP
  (§J.6.2) or use a milder DTT concentration (0.5 mM).

---

### J.6.2  Protein Pretreatment — Disulfide reduction with TCEP

**Purpose:** Tris(2-carboxyethyl)phosphine. Reduces disulfides without
needing removal before thiol-maleimide coupling (TCEP doesn't react with
maleimide at normal concentrations). Stronger reductant than DTT. Useful
for disulfide-stable proteins.

**Evidence tier:** VALIDATED

#### Safety

- TCEP·HCl: H315, H319.
- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: aqueous, routine.

#### What you need

| Reagent | Amount | CAS |
|---|---|---|
| TCEP·HCl | 0.5-5 mM final | 51805-45-9 |
| Reducing buffer pH 7.5 | 10 mL | — |

#### Procedure

1. Dissolve protein in reducing buffer.
2. Add TCEP from 500 mM stock (prepared fresh or stored at –80 °C) to
   0.5-5 mM final.
3. Incubate 15-30 min RT.
4. Proceed DIRECTLY to thiol coupling. No buffer exchange needed (unlike DTT).

---

### J.6.3  Protein Pretreatment — Glycan oxidation with NaIO₄

**Purpose:** Convert diols on glycoprotein glycans to aldehydes. Required
for oriented IgG immobilisation (§J.3.5).

**Evidence tier:** VALIDATED

#### Safety

- NaIO₄: H272 (oxidiser), H315, H319. Keep dry, away from organics.
- PPE: nitrile gloves, safety glasses, lab coat.
- Work on ice, in the dark.
- Waste: aqueous; neutralise with ethylene glycol before disposal.

#### What you need (per 5 mg glycoprotein)

| Reagent | Amount | CAS |
|---|---|---|
| NaIO₄ | 10 mM final | 7790-28-5 |
| PBS | 5 mL | — |
| Ethylene glycol (quench) | 100 mM final | 107-21-1 |

#### Procedure

1. Dissolve protein at 1 mg/mL in PBS.
2. Add NaIO₄ from 100 mM stock to 10 mM final. On ice, in the dark, 30 min.
3. Quench with ethylene glycol to 100 mM, 15 min.
4. Buffer-exchange to coupling buffer (usually pH 6 phosphate for hydrazone,
   §J.6.4).

---

### J.6.4  Protein Pretreatment — Buffer exchange with PD-10 / Sephadex G-25

**Purpose:** Remove small molecules (amines, reducing agents, salt, dye)
before coupling. Fastest, cheapest method for 0.5-2.5 mL protein samples.

**Evidence tier:** VALIDATED

#### Safety

- No hazardous reagents.
- PPE: nitrile gloves, safety glasses, lab coat.
- Waste: aqueous, routine.

#### What you need

| Reagent | Amount |
|---|---|
| PD-10 column (or equivalent G-25 column) | 1 column |
| Target buffer | 25 mL |

#### Procedure

1. Equilibrate PD-10 with 25 mL target buffer (gravity flow).
2. Apply 2.5 mL protein sample, let enter bed.
3. Wash with 2.5 mL target buffer (discarded).
4. Elute with 3.5 mL target buffer (collected).
5. Protein is in a new buffer, diluted ~1.4-fold.

#### QC

- Verify buffer exchange by pH / conductivity measurement on the eluate.

---

### J.6.5  Protein Pretreatment — Dialysis

**Purpose:** Gentle, slow buffer exchange. Better for aggregation-prone
proteins where PD-10 handling causes denaturation.

**Evidence tier:** VALIDATED

#### Safety

- No hazardous reagents.
- PPE: nitrile gloves, lab coat.
- Waste: aqueous.

#### What you need

| Reagent | Amount |
|---|---|
| Dialysis tubing or cassette, appropriate MWCO (typically 10 kDa for proteins > 30 kDa) | 1 |
| Target buffer | 1 L per 10 mL sample |

#### Procedure

1. Pre-soak dialysis tubing 10 min in water (if dry tubing). Cassettes need
   no pre-soak.
2. Load protein. Place in 1 L target buffer with gentle stir bar at 4 °C.
3. Exchange 3 × 1 L buffer over 24 h (e.g., 4 h, then overnight, then 4 h
   again). Protein concentration unchanged.

---

### J.6.6  Protein Pretreatment — Concentration (Amicon / centrifugal)

**Purpose:** Raise protein concentration to ≥ 1 mg/mL for efficient coupling.

**Evidence tier:** VALIDATED

#### Safety / Waste

- Standard lab practice; no special hazards.

#### What you need

| Reagent | Amount |
|---|---|
| Centrifugal concentrator (appropriate MWCO, typically MWCO = 1/3 × protein MW) | 1 |
| Cold centrifuge | — |

#### Procedure

1. Pre-rinse concentrator with 1 mL buffer (removes glycerol residue).
2. Load protein, centrifuge at 4000 × g (or manufacturer spec), 4 °C.
3. Resuspend concentrate every 10 min to prevent membrane fouling.
4. Stop when target concentration reached. Measure A₂₈₀ or Bradford.

---

## J.7 Washing

Wash steps remove unbound reagents, non-covalently adsorbed contaminants,
and buffer residues. "Over-washing" is rare; "under-washing" is a common
cause of coupling-density overestimation and non-specific binding.

Universal rules:
- Wash at ≥ 3 CV (column volumes) per step.
- Use the same temperature as the step it follows.
- Never let the gel run dry on the funnel (air ingress damages matrix and
  can collapse pores irreversibly).

---

### J.7.1  Washing — Post-activation rinse

**Purpose:** Remove unreacted activator (CNBr, epichlorohydrin, DVS, CDI,
NHS-ester activator) before the gel sits for coupling.

**Buffer:** same chemistry as coupling buffer for the next step. E.g.,
ice-cold 0.1 M NaHCO₃ pH 8.3 for NHS-ester coupling; ice-cold 1 mM HCl for
CNBr-activated gels (preserves activation).

**Volume:** 3-5 CV.
**Temperature:** ice-cold where specified (CNBr, DVS), RT otherwise.

---

### J.7.2  Washing — Post-coupling rinse

**Purpose:** Remove unreacted ligand/protein, leaving only covalently bound.

**Buffer:** coupling buffer (e.g., 0.1 M NaHCO₃ pH 8.3).
**Volume:** 3 CV.

Save the first CV — measure ligand concentration in it to compute coupling
yield by difference (starting ligand minus wash-through ligand = incorporated).

---

### J.7.3  Washing — High-salt wash

**Purpose:** Disrupt ionic interactions between unreacted ligand/protein
and the matrix. Particularly important after CNBr-coupling (isourea is
positively charged) and for hydrophobic matrices.

**Buffer:** PBS + 0.5-1 M NaCl.
**Volume:** 3-5 CV.
**Temperature:** RT.

---

### J.7.4  Washing — Low-pH wash (for affinity media)

**Purpose:** Disrupt residual non-covalent binding on affinity media (e.g.,
Protein A/G columns post-DMP crosslinking, IgG-bound matrices).

**Buffer:** 0.1 M glycine pH 2.5 OR 0.1 M citric acid pH 3.0.
**Volume:** 3 CV.
**Temperature:** RT.
**Key note:** immediately neutralise the gel with 1 M Tris pH 9.0 (1/10
volume) in the collection tube so prolonged low pH doesn't damage the
ligand.

---

### J.7.5  Washing — Detergent wash

**Purpose:** Remove hydrophobically-adsorbed contaminants (lipids,
denatured protein). Use for hydrophobic ligands or when high background
is observed.

**Buffer:** PBS + 0.05 % Tween-20 OR PBS + 0.1 % Triton X-100.
**Volume:** 3 CV.
**Temperature:** RT.
**Follow-up:** always wash with 3 × 3 CV of plain PBS to remove detergent
residue (interferes with subsequent assays).

---

### J.7.6  Washing — Storage equilibration

**Purpose:** Place the finished gel into a stable storage buffer.

**Buffers:**
- For most protein-bearing affinity media: PBS + 0.05 % NaN₃ (antimicrobial)
  at 4 °C. Shelf life: 6-12 months.
- For organic-solvent-stable ligands: 20 % ethanol at 4 °C.
- For IMAC matrix: 20 % ethanol (stripped) or the charged form in binding
  buffer at 4 °C (only if used within days).

**Never freeze an agarose matrix.** Freezing destroys pore structure
irreversibly.

**Label** the storage tube with: date, matrix type, ligand, coupling
density (μmol/mL), buffer, storage T.

---

## J.8 Quenching

Quenching deactivates unreacted activated sites on the matrix after the
ligand has coupled. This is essential — unreacted sites continue to react
(slowly) with any amine, thiol, or other nucleophile that contacts the
gel, including your downstream samples. That produces non-specific
capture, reduced active-site accessibility, and irreproducible binding.

Quench by adding an excess of a harmless small-molecule nucleophile
that will saturate all remaining active sites.

---

### J.8.1  Quenching — NHS-ester / CDI / p-nitrophenyl carbonate

**Purpose:** Deactivate unreacted NHS-ester, CDI, or p-NP-carbonate sites.

**Reagent:** 1 M ethanolamine OR 1 M Tris, pH 8.0.

**Procedure:** drain post-coupling gel. Add 2 CV of 1 M Tris pH 8.0 or
1 M ethanolamine pH 8.0. Rotate 30-60 min at RT. Wash with 3 CV coupling
buffer to remove excess quench reagent.

**Safety:** Ethanolamine H302, H314. PPE: nitrile, face shield. Waste:
aqueous organic, neutralise.

---

### J.8.2  Quenching — Epoxide (BDDE, epichlorohydrin)

**Purpose:** Deactivate unreacted epoxide sites.

**Reagent:** 1 M ethanolamine pH 8.0 (for amine-quench) OR 0.5 M
β-mercaptoethanol (for thiol-quench).

**Procedure:** drain post-coupling gel. Add 2 CV of quench reagent. Rotate
4 h RT to overnight (epoxide is slow). Wash with 3 CV coupling buffer,
3 CV water.

**Safety:** β-mercaptoethanol H301, H310, H315, H317, H330, H410. Fume hood
mandatory. Ethanolamine as above. Waste: aqueous organic (neutralise and
dilute).

---

### J.8.3  Quenching — DVS (vinyl sulfone)

**Purpose:** Deactivate unreacted vinyl-sulfone sites.

**Reagent:** 1 M ethanolamine pH 9.0 OR 0.5 M β-mercaptoethanol in PBS.

**Procedure:** drain post-coupling gel. Add 2 CV quench. Rotate 2 h RT.
Wash with 3 CV PBS.

**Safety:** as §J.8.2. β-mercaptoethanol is the gold standard for vinyl
sulfone quenching (reacts fast); ethanolamine is slower but less hazardous.

---

### J.8.4  Quenching — Maleimide

**Purpose:** Deactivate unreacted maleimide sites.

**Reagent:** 10 mM L-cysteine OR 10 mM β-mercaptoethanol in PBS.

**Procedure:** drain post-coupling gel. Add 2 CV quench reagent. Rotate
30 min RT. Wash with 3 CV PBS.

**Safety:** cysteine low-hazard. β-mercaptoethanol as above.

---

### J.8.5  Quenching — Aldehyde (periodate-oxidised or glutaraldehyde-activated)

**Purpose:** Deactivate unreacted aldehyde sites AND reduce any existing
Schiff bases to stable secondary amines in one step.

**Reagent:** 50 mM NaBH₄ in 0.1 M Na phosphate pH 7.4 (fresh).
Alternative for less-damaging reduction: 50 mM NaCNBH₃ in 0.1 M Na
phosphate pH 7.4.

**Procedure:** drain post-coupling gel. Add 2 CV reducing agent + 100 mM
glycine or ethanolamine (to also cap remaining aldehydes). Rotate 30 min
at RT. Wash with 3 CV PBS.

**Safety:** NaBH₄ H260 (contact with water releases flammable gas), H314.
NaCNBH₃: H300+H311+H331. Fume hood mandatory for both. Waste: NaBH₄ aqueous
— let any residual gas evolve in the hood before capping the waste bottle.
NaCNBH₃ aqueous per §J.2.4 (alkaline + bleach neutralisation).

---

## J.bonus  Sim-to-Bench Decision Tree

A short index showing how to go from the simulator's output to the bench
protocol a first-timer should actually execute. If your M1 simulation
shipped a microsphere with target d32, pore size, and modulus, and you
want a surface functionality X:

```
Target functionality         →  Protocol chain

Affinity capture of IgG      →  §J.1.1 (CNBr)   +  §J.3.3 (Protein A + DMP) [oriented]
                             or  §J.1.5 (CDI)   +  §J.3.1 (NHS-ester)       [non-oriented]

Metal-binding protein (His-tag)→§J.1.3 (BDDE)   +  §J.2.2 (amine ligand: NTA)
                             +  §J.5.1 (Ni²⁺ charging)

Enzyme immobilisation         →  §J.1.2 (epichlorohydrin)  +  §J.3.2 (glutaraldehyde)
                             or  §J.1.4 (DVS)   +  §J.2.3 (thiol on cysteine)

Oriented IgG immobilisation  →  §J.6.3 (glycan oxidation) + §J.3.5 (Fc-hydrazone)

Click-functionalised matrix  →  §J.1.4 (DVS)   +  §J.2.3 (attach azide-thiol)
                             or vendor azide-agarose +  §J.2.5 (SPAAC)
```

Always end the chain with a §J.7 wash sequence and a §J.8 quench matched to
the last active chemistry. Store per §J.7.6.

---

## References consolidated

- Axén R, Porath J, Ernback S (1967) *Nature* 214:1302-1304.
- Porath J, Låås T, Janson J-C (1975) *J. Chromatogr.* 103:49-62.
- Porath J et al. (1975) *Nature* 258:598-599. (IMAC)
- Cuatrecasas P (1970) *J. Biol. Chem.* 245:3059-3065.
- Cuatrecasas P, Anfinsen CB (1971) *Annu. Rev. Biochem.* 40:259-278.
- Sundberg L, Porath J (1974) *J. Chromatogr.* 90:87-98.
- Matsumoto I, Mizuno Y, Seno N (1980) *J. Chromatogr.* 188:457-464.
- Nilsson K, Mosbach K (1984) *Meth. Enzymol.* 104:56-69.
- Kohn J, Wilchek M (1984) *Appl. Biochem. Biotechnol.* 9:285-305.
- O'Shannessy DJ et al. (1984) *J. Immunol. Methods* 75:11-17.
- Andersson L, Porath J (1986) *Anal. Biochem.* 154:250-254.
- Hearn MTW (1987) *Meth. Enzymol.* 135:102-117.
- Hochuli E, Döbeli H, Schacher A (1988) *Biotechnology* 6:1321-1325.
- Staros JV (1982) *Biochemistry* 21:3950-3955.
- Schneider C et al. (1982) *J. Biol. Chem.* 257:10766-10769.
- Monsan P (1978) *J. Mol. Catal.* 3:371-384.
- Kolb HC, Finn MG, Sharpless KB (2001) *Angew. Chem. Int. Ed.* 40:2004-2021.
- Meldal M, Tornøe CW (2008) *Chem. Rev.* 108:2952-3015.
- Mateo C et al. (2007) *Enzyme Microb. Technol.* 39:274-280.
- Morpurgo M et al. (1996) *Bioconjug. Chem.* 7:363-368.
- Mao H et al. (2004) *J. Am. Chem. Soc.* 126:2670-2671.
- Popp MW, Ploegh HL (2011) *Angew. Chem. Int. Ed.* 50:5024-5032.
- Bethell GS et al. (1979) *J. Biol. Chem.* 254:2572-2574.
- Migneault I et al. (2004) *BioTechniques* 37:790-802.
- Harris JM (ed.) (1992) *Poly(ethylene glycol) Chemistry*, Plenum Press.
- **Hermanson GT (2013) *Bioconjugate Techniques*, 3rd ed., Academic Press.**

---

## Disclaimer

This appendix is provided for informational, research, and training
purposes only. It does not constitute professional engineering advice,
medical advice, or formal peer review. Every protocol in this appendix
uses reagents some of which are classified as carcinogenic, mutagenic,
reproductively toxic, or strongly sensitising in one or more jurisdictions.
Before handling any reagent listed here, users must:

1. Consult the institution's safety office and the reagent's SDS.
2. Obtain approval for any restricted / CMR substance per local policy.
3. Receive hands-on training from a qualified supervisor.
4. Use the fume hood and PPE specified in each protocol without substitution.

The author is an AI assistant; all protocols should be validated on small
scale by a qualified researcher before routine use. Record every run in a
laboratory notebook. If your wet-lab outcome differs from the simulator's
prediction by more than a factor of 2, consult the Calibration panel of
EmulSim (§ Chapter 4 of the main manual) to update the model against your
actual bench data.

*End of Appendix J.*
