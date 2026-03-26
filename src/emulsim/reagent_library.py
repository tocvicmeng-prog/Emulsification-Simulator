"""Library of alternative crosslinkers and surfactants with simulation parameters.

Each entry contains literature-sourced or first-principles-estimated values
for use in the emulsification-gelation-crosslinking simulation pipeline.
Values marked [estimated] lack direct literature support and were derived
from analogous systems or theoretical scaling.

References
----------
Genipin:
    Butler et al. (2003) J. Polym. Sci. A 41:3941-3953
    Mi et al. (2005) Biomaterials 26:5983-5990
    Muzzarelli (2009) Carbohydr. Polym. 77:1-9
Glutaraldehyde:
    Migneault et al. (2004) BioTechniques 37:790-802
    Monteiro & Bhatt (2012) MRS Proc. 1498:mrss12-1498-l09-02
    Roberts & Taylor (1989) Makromol. Chem. 190:951-960
EDC/NHS:
    Hermanson (2013) Bioconjugate Techniques, 3rd ed., Academic Press
    Sehgal & Bhave (2020) J. Biomed. Mater. Res. A 108:1507
PEGDA:
    Lin & Anseth (2009) Pharm. Res. 26:631-643
    Nguyen & West (2002) Biomaterials 23:4307-4314
TPP:
    Calvo et al. (1997) J. Appl. Polym. Sci. 63:125-132
    Shu & Zhu (2002) Int. J. Pharm. 233:217-225
Epichlorohydrin:
    Zhao et al. (2020) Carbohydr. Polym. 227:115352
    Chen et al. (2016) Cellulose 23:631-639
Divinyl Sulfone (DVS):
    Porath et al. (1975) J. Chromatogr. 103:49-62
    Mateo et al. (2006) Enzyme Microb. Technol. 39:274-280
Citric Acid:
    Demitri et al. (2008) J. Appl. Polym. Sci. 110:2453-2460
    Reddy & Yang (2010) Food Chem. 118:702-711

Surfactant references:
    Santini et al. (2007) Colloids Surf. A 298:12-21
    Opawale & Burgess (1998) J. Colloid Interface Sci. 197:142-150
    Peltonen & Yliruusi (2000) J. Colloid Interface Sci. 227:1-6
    Wilson et al. (1998) J. Colloid Interface Sci. 205:94-103
    Marze (2009) Food Hydrocoll. 23:1972-1980
    Shrestha et al. (2010) Langmuir 26:7015-7024
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class CrosslinkerProfile:
    """Profile for a crosslinker candidate in the chitosan network.

    Parameters map directly to Level 3 solver inputs.
    """

    name: str
    cas: str
    mechanism: str  # "amine_bridge", "uv_radical", "ionic", "hydroxyl"
    k_xlink_0: float  # [m3/(mol*s)] Arrhenius rate prefactor
    E_a_xlink: float  # [J/mol] activation energy
    f_bridge: float  # [-] fraction of reactions producing elastically active crosslinks
    T_crosslink_default: float  # [K] recommended crosslinking temperature
    t_crosslink_default: float  # [s] recommended crosslinking time
    kinetics_model: str  # "second_order", "michaelis_menten", "uv_dose", "ionic_instant"
    suitability: int  # 1-10, overall suitability for agarose-chitosan chromatography beads
    notes: str


@dataclass
class SurfactantProfile:
    """Profile for a surfactant candidate for W/O emulsification.

    Szyszkowski-Langmuir parameters for use in the interfacial tension model:
        sigma = sigma_0 - R*T*gamma_inf*ln(1 + K_L*c_mol)
    """

    name: str
    cas: str
    hlb: float  # [-] hydrophilic-lipophilic balance
    mw: float  # [g/mol] molecular weight
    gamma_inf: float  # [mol/m2] maximum surface excess (Gibbs)
    K_L: float  # [m3/mol] Langmuir adsorption constant
    sigma_0_paraffin: float  # [N/m] bare IFT with paraffin oil (no surfactant, 20 degC)
    suitability: int  # 1-10
    notes: str


# ═══════════════════════════════════════════════════════════════════════════
#  CROSSLINKER LIBRARY
# ═══════════════════════════════════════════════════════════════════════════

CROSSLINKERS: dict[str, CrosslinkerProfile] = {
    # ── 1. Genipin (baseline) ─────────────────────────────────────────────
    "genipin": CrosslinkerProfile(
        name="Genipin",
        cas="6902-77-8",
        mechanism="amine_bridge",
        # Arrhenius prefactor calibrated so k(37 degC) ~ 5e-6 m3/(mol*s)
        # = 5e-3 L/(mol*s). From Butler et al. (2003): k ~ 0.005 L/(mol*s)
        # at 37 degC, Ea ~ 50-55 kJ/mol.
        k_xlink_0=2806.0,
        E_a_xlink=52000.0,
        # f_bridge 0.3-0.5 per Mi et al. (2005); genipin forms both
        # intra-chain ring opening products (monomeric) and inter-chain
        # bridges (oligomeric blue pigment).
        f_bridge=0.40,
        T_crosslink_default=310.15,  # 37 degC
        t_crosslink_default=86400.0,  # 24 h
        kinetics_model="second_order",
        suitability=8,
        notes=(
            "Natural crosslinker from gardenia fruit. Low cytotoxicity "
            "(~10000x less than glutaraldehyde). Forms stable amide and "
            "heterocyclic crosslinks with chitosan primary amines. Produces "
            "blue coloration (may affect optical monitoring). Slow reaction "
            "kinetics (hours to days). Does NOT interfere with agarose helix "
            "formation — genipin is amine-specific and does not react with "
            "agarose hydroxyls under mild conditions. Cost: ~USD 5-15/g."
        ),
    ),
    # ── 2. Glutaraldehyde ─────────────────────────────────────────────────
    "glutaraldehyde": CrosslinkerProfile(
        name="Glutaraldehyde",
        cas="111-30-8",
        mechanism="amine_bridge",
        # Glutaraldehyde reacts ~100-1000x faster than genipin with amines.
        # Roberts & Taylor (1989): k ~ 1 L/(mol*s) at 25 degC for amine-aldehyde
        # Schiff base formation. Ea ~ 40 kJ/mol (lower barrier, aldehyde
        # electrophilicity). k0 = k(25C)*exp(Ea/(R*298)) ~ 9.5e4 m3/(mol*s).
        k_xlink_0=9.5e4,
        E_a_xlink=40000.0,
        # Higher bridge efficiency — bifunctional aldehyde, nearly every
        # molecule forms inter-chain bridges (Schiff base at each end).
        f_bridge=0.65,
        T_crosslink_default=298.15,  # 25 degC — fast enough at RT
        t_crosslink_default=3600.0,  # 1 h sufficient
        kinetics_model="second_order",
        suitability=5,
        notes=(
            "OHC-(CH2)3-CHO, bifunctional aldehyde. Forms Schiff base "
            "(C=N) linkages with chitosan -NH2. Extremely fast, cheap "
            "(~USD 0.05/g), well-characterized. HIGH TOXICITY: LD50 oral "
            "~100 mg/kg, cytotoxic residuals require extensive washing. "
            "Schiff bases are hydrolytically unstable at low pH — may need "
            "NaBH4 reduction to stable secondary amines. Does not interfere "
            "with agarose helix formation. Produces yellow coloration. "
            "For chromatography: residual glutaraldehyde may leach and "
            "modify bound proteins — significant concern for bioprocessing."
        ),
    ),
    # ── 3. EDC/NHS (carbodiimide) ─────────────────────────────────────────
    "edc_nhs": CrosslinkerProfile(
        name="EDC/NHS (1-ethyl-3-(3-dimethylaminopropyl)carbodiimide / N-hydroxysuccinimide)",
        cas="25952-53-8 / 6066-82-6",
        mechanism="amine_bridge",
        # EDC activates carboxyl groups which then react with amines.
        # Zero-length crosslinker — no bridge molecule remains.
        # Requires COOH groups: chitosan has few unless partially reacetylated
        # or blended with carboxymethyl chitosan.
        # EDC hydrolysis competes: t1/2 ~ 10-20 min at pH 5 (Hermanson 2013).
        # Effective coupling k ~ 0.01 L/(mol*s) at pH 5, 25 degC [estimated
        # from typical EDC/NHS amide bond yields].
        # Ea ~ 45 kJ/mol [estimated from temperature studies in Sehgal 2020].
        k_xlink_0=5.0e4,
        E_a_xlink=45000.0,
        # Bridge efficiency is HIGH for available COOH/NH2 pairs, but
        # chitosan has very few COOH groups (~1-5% depending on source).
        # Effective f_bridge refers to fraction of EDC-activated sites
        # that form stable amide bonds (vs hydrolysis). Typically 10-40%.
        f_bridge=0.25,
        T_crosslink_default=277.15,  # 4 degC — cold to minimize EDC hydrolysis
        t_crosslink_default=14400.0,  # 4 h
        kinetics_model="michaelis_menten",
        suitability=3,
        notes=(
            "Zero-length crosslinker: forms direct amide bonds, no spacer "
            "molecule remains in the gel. REQUIRES carboxyl groups — pure "
            "chitosan has almost none. Would need carboxymethyl chitosan or "
            "addition of a dicarboxylic acid (e.g., adipic acid) as spacer. "
            "EDC hydrolyzes rapidly in water (t1/2 ~ 10-20 min at pH 5), "
            "limiting efficiency. NHS ester intermediate improves yield. "
            "Mild conditions (pH 4.5-6, 4-25 degC). Not directly applicable "
            "to standard chitosan without formulation changes. Does not "
            "interfere with agarose helices. Cost moderate (~USD 1-3/g EDC)."
        ),
    ),
    # ── 4. PEGDA + UV ─────────────────────────────────────────────────────
    "pegda_uv": CrosslinkerProfile(
        name="PEGDA (poly(ethylene glycol) diacrylate) + UV photoinitiation",
        cas="26570-48-9",
        mechanism="uv_radical",
        # UV-initiated radical polymerization of acrylate end groups.
        # Not a direct amine crosslinker — works by forming a PEG network
        # that interpenetrates the chitosan/agarose matrix.
        # k_p (propagation) ~ 10-100 L/(mol*s) for acrylates at 25 degC.
        # But overall crosslinking rate depends on UV dose and initiator.
        # Effective k_xlink_0 for simulation: use a dose-dependent model.
        # k_eff = k0 * I_UV^0.5 (radical generation scales as sqrt(intensity))
        # k0 ~ 1e3 m3/(mol*s) [estimated], Ea ~ 20 kJ/mol (radical, low barrier)
        k_xlink_0=1.0e3,
        E_a_xlink=20000.0,
        # PEG chains form a network — bridge efficiency depends on PEGDA MW
        # and concentration. For MW 700: ~0.5-0.7.
        f_bridge=0.55,
        T_crosslink_default=298.15,  # RT — UV driven
        t_crosslink_default=600.0,  # 10 min UV exposure
        kinetics_model="uv_dose",
        suitability=4,
        notes=(
            "Photo-crosslinkable: spatial and temporal control via UV lamp. "
            "Does NOT react with chitosan -NH2 directly — forms an "
            "independent PEG network (triple IPN: agarose + chitosan + PEG). "
            "Requires photoinitiator (Irgacure 2959 or LAP) which may be "
            "cytotoxic. UV penetration limited in opaque microspheres — "
            "core may be under-crosslinked for beads >50 um. Agarose helix "
            "formation not affected. PEG network adds hydrophilicity and "
            "reduces non-specific protein binding (advantage for chromato). "
            "Simulation needs UV dose model: k_eff = k0 * (I_UV)^0.5 * "
            "exp(-mu*r) where mu is optical attenuation. Cost: ~USD 0.5/g."
        ),
    ),
    # ── 5. Sodium tripolyphosphate (TPP) ──────────────────────────────────
    "tpp": CrosslinkerProfile(
        name="Sodium tripolyphosphate (TPP)",
        cas="7758-29-4",
        mechanism="ionic",
        # Ionic crosslinking: P3O10^5- bridges NH3+ groups on chitosan.
        # Essentially instantaneous electrostatic interaction — not a
        # covalent bond. "Rate" is diffusion-limited, not reaction-limited.
        # For simulation: k_xlink_0 is very large; gelation is nearly instant
        # upon contact. Use a high k0 with low Ea.
        # Calvo et al. (1997): TPP/chitosan nanoparticles form in seconds.
        k_xlink_0=1.0e8,
        E_a_xlink=10000.0,  # [estimated] minimal thermal barrier, mostly diffusion
        # Ionic crosslinks are multivalent but individually weak.
        # f_bridge is high (most TPP ions bridge chains) but bonds are
        # reversible and pH-sensitive.
        f_bridge=0.70,
        T_crosslink_default=298.15,  # RT
        t_crosslink_default=300.0,  # 5 min (effectively instant)
        kinetics_model="ionic_instant",
        suitability=4,
        notes=(
            "Ionic (non-covalent) crosslinker. P3O10^5- bridges protonated "
            "-NH3+ groups at pH < 6.5. Extremely fast gelation (seconds). "
            "Very cheap, food-grade, non-toxic. HOWEVER: ionic crosslinks "
            "are REVERSIBLE — sensitive to pH, ionic strength, and competing "
            "ions. At chromatography conditions (buffered salt solutions), "
            "TPP crosslinks will dissociate. Mechanical strength is low "
            "(G ~ 0.1-1 kPa vs 10-100 kPa for covalent). Not suitable as "
            "sole crosslinker for chromatography beads under pressure. "
            "Does not interfere with agarose. Simulation needs ionic "
            "equilibrium model rather than irreversible kinetics. "
            "Cost: ~USD 0.01/g."
        ),
    ),
    # ── 6. Epichlorohydrin ────────────────────────────────────────────────
    "epichlorohydrin": CrosslinkerProfile(
        name="Epichlorohydrin (ECH)",
        cas="106-89-8",
        mechanism="hydroxyl",
        # Reacts with -OH groups under alkaline conditions (pH > 10).
        # Used in Sepharose/Sephadex manufacturing (GE Healthcare).
        # Zhao et al. (2020): ECH crosslinks both agarose -OH and chitosan
        # -OH/-NH2 groups. Requires NaOH activation.
        # Rate: slower than glutaraldehyde, faster than genipin.
        # k ~ 0.01-0.1 L/(mol*s) at 50 degC, pH 12 (Chen et al. 2016).
        # Ea ~ 60 kJ/mol (epoxide ring-opening, higher barrier).
        # k0 = k(50C)*exp(Ea/(R*323)) ~ 1.5e7 m3/(mol*s)
        k_xlink_0=1.5e7,
        E_a_xlink=60000.0,
        # ECH crosslinks BOTH networks — reacts with agarose -OH and
        # chitosan -OH/-NH2. This is both advantage (reinforces both)
        # and complication (may interfere with agarose helix structure
        # if applied before gelation).
        f_bridge=0.50,
        T_crosslink_default=323.15,  # 50 degC
        t_crosslink_default=21600.0,  # 6 h
        kinetics_model="second_order",
        suitability=7,
        notes=(
            "Epoxide crosslinker used industrially in Sepharose CL "
            "production (Pharmacia/Cytiva). Reacts with -OH groups under "
            "strongly alkaline conditions (0.5-2 M NaOH, pH > 12). Also "
            "reacts with -NH2 at lower pH. CRITICAL: alkaline conditions "
            "required for -OH crosslinking WILL partially dissolve agarose "
            "helices — must crosslink AFTER complete gelation, or use a "
            "two-step protocol (gel first, then ECH at mild alkaline pH "
            "targeting -NH2 preferentially). Produces ether bridges "
            "(-O-CH2-CHOH-CH2-O-) that are extremely hydrolytically stable. "
            "TOXIC and CARCINOGENIC (IARC Group 1) — requires careful "
            "handling and thorough washing. Widely accepted in chromatography "
            "resin manufacturing with validated wash protocols. "
            "Cost: ~USD 0.02/g."
        ),
    ),
    "dvs": CrosslinkerProfile(
        name="Divinyl Sulfone (DVS)",
        cas="77-77-0",
        mechanism="michael_addition",
        # DVS reacts with -OH groups via Michael addition under alkaline
        # conditions (pH 11-12).  Porath et al. (1975) J. Chromatogr. 103:49.
        # Rate is faster than ECH at equivalent pH — typically 1-2 h at 25°C.
        # k0 estimated from published activation kinetics (Mateo et al. 2006,
        # Enzyme Microb. Technol. 39:274).
        k_xlink_0=5000.0,      # [m³/(mol·s)] [estimated]
        E_a_xlink=55000.0,     # [J/mol] [estimated from ECH analogy]
        # DVS produces highly efficient ether-sulfone bridges between
        # polysaccharide chains — virtually every reaction creates a crosslink.
        f_bridge=0.70,
        T_crosslink_default=298.15,  # 25 degC (room temperature, alkaline)
        t_crosslink_default=7200.0,   # 2 h
        kinetics_model="second_order",
        suitability=9,
        notes=(
            "Industry standard for high-performance agarose chromatography "
            "media (Sepharose HP, Sepharose CL). Reacts with -OH groups "
            "under alkaline conditions (pH 11-12, 0.5 M NaOH + 0.5 M "
            "Na2CO3). Produces vinyl sulfone ether bridges that are "
            "extremely stable (resistant to acid, base, urea, guanidine). "
            "Superior rigidity and flow rate vs uncrosslinked or ECH-"
            "crosslinked agarose. TOXIC and LACHRYMATORY — handle in fume "
            "hood with respiratory protection. Well-validated wash protocols "
            "in chromatography resin manufacturing. "
            "Ref: Porath et al. (1975) J. Chromatogr. 103:49-62. "
            "Cost: ~USD 0.50/g (Sigma-Aldrich V3700)."
        ),
    ),
    "citric_acid": CrosslinkerProfile(
        name="Citric Acid",
        cas="77-92-9",
        mechanism="ester_bond",
        # Citric acid crosslinks -NH2 and -OH groups via ester/amide bond
        # formation under heat (80-120°C).  Reddy & Yang (2010) Food Chem.
        # 118:702.  Reaction requires elevated temperature and is slow at RT.
        # k0 estimated from cure kinetics in Demitri et al. (2008).
        k_xlink_0=50.0,        # [m³/(mol·s)] [estimated, heat-activated]
        E_a_xlink=80000.0,     # [J/mol] [estimated — high Ea due to ester formation]
        # Trifunctional (3 COOH groups) but steric hindrance limits bridge
        # formation — many pendant carboxyls remain unreacted.
        f_bridge=0.25,
        T_crosslink_default=353.15,  # 80 degC (heat cure)
        t_crosslink_default=14400.0,  # 4 h
        kinetics_model="second_order",
        suitability=5,
        notes=(
            "Green chemistry crosslinker — cheap, non-toxic, FDA GRAS. "
            "Trifunctional (3 COOH groups) forms ester bonds with -OH and "
            "amide bonds with -NH2 under heat (80-120°C). CAUTION: heat "
            "cure at 80°C+ may partially melt agarose gel if performed "
            "before cooling below T_gel. Ester bonds are hydrolytically "
            "unstable at pH > 10 or < 3 — limits use in strong ion-exchange "
            "chromatography regeneration protocols. Suitable for "
            "single-use or mild-condition applications. "
            "Ref: Demitri et al. (2008) J. Appl. Polym. Sci. 110:2453. "
            "Cost: ~USD 0.01/g (commodity chemical)."
        ),
    ),
}


# ═══════════════════════════════════════════════════════════════════════════
#  SURFACTANT LIBRARY
# ═══════════════════════════════════════════════════════════════════════════

SURFACTANTS: dict[str, SurfactantProfile] = {
    # ── 1. Span-80 (baseline) ─────────────────────────────────────────────
    "span80": SurfactantProfile(
        name="Span-80 (sorbitan monooleate)",
        cas="1338-43-8",
        hlb=4.3,
        mw=428.6,
        # Santini et al. (2007): Gamma_inf = 3.5e-6 mol/m2 from pendant drop
        # tensiometry at paraffin oil/water interface.
        gamma_inf=3.5e-6,
        # K_L calibrated to give sigma ~ 5 mN/m at 2% w/v Span-80 at 90 degC
        # (Opawale & Burgess 1998).
        K_L=0.75,
        sigma_0_paraffin=0.050,  # bare water/paraffin oil IFT at 20 degC
        suitability=8,
        notes=(
            "Industry-standard W/O emulsifier. Liquid at RT (mp -25 degC), "
            "easy to handle. Oleate chain provides good oil solubility. "
            "Well-characterized IFT and CMC data. Compatible with agarose — "
            "does not interfere with helix formation. Produces stable "
            "emulsions in the 20-200 um range at 1-5% w/v. Biodegradable, "
            "food-grade (E494). Moderate IFT reduction. Cost: ~USD 0.10/g."
        ),
    ),
    # ── 2. Span-60 (sorbitan monostearate) ────────────────────────────────
    "span60": SurfactantProfile(
        name="Span-60 (sorbitan monostearate)",
        cas="1338-41-6",
        hlb=4.7,
        mw=430.6,
        # Similar headgroup to Span-80, slightly different packing due to
        # saturated C18 chain. Gamma_inf ~ 3.2e-6 mol/m2 [estimated from
        # molecular area ~52 A2, close to Span-80 at ~47 A2].
        gamma_inf=3.2e-6,  # [estimated]
        # K_L slightly lower due to reduced oil solubility of saturated chain.
        K_L=0.55,  # [estimated]
        sigma_0_paraffin=0.050,
        suitability=5,
        notes=(
            "Saturated C18 analog of Span-80. SOLID at RT (mp 53-57 degC) — "
            "must be dissolved in hot oil. At 90 degC emulsification "
            "temperature, fully dissolved and functional. Forms more rigid "
            "interfacial films than Span-80 (crystalline packing of "
            "saturated chains). May produce more stable initial emulsions "
            "but RISK of surfactant crystallization during cooling, which "
            "could disrupt microsphere surface or cause aggregation. "
            "Less IFT reduction than Span-80 at equal concentration. "
            "Cost: ~USD 0.08/g."
        ),
    ),
    # ── 3. Span-40 (sorbitan monopalmitate) ───────────────────────────────
    "span40": SurfactantProfile(
        name="Span-40 (sorbitan monopalmitate)",
        cas="26266-57-9",
        hlb=6.7,
        mw=402.6,
        # Higher HLB — less lipophilic. Gamma_inf lower due to shorter C16
        # chain providing less anchoring in oil phase.
        gamma_inf=2.8e-6,  # [estimated]
        K_L=0.35,  # [estimated] — weaker adsorption from oil phase
        sigma_0_paraffin=0.050,
        suitability=3,
        notes=(
            "C16 saturated sorbitan ester. HLB 6.7 is borderline for W/O "
            "emulsions — closer to the W/O-to-O/W transition (~7-8). At "
            "elevated temperature, effective HLB shifts further toward O/W, "
            "risking phase inversion during emulsification at 90 degC. "
            "Solid at RT (mp 45-47 degC). Less effective IFT reduction than "
            "Span-80. NOT RECOMMENDED as sole emulsifier; may work in "
            "blends with Span-85 or PGPR to tune HLB. Cost: ~USD 0.08/g."
        ),
    ),
    # ── 4. Span-85 (sorbitan trioleate) ───────────────────────────────────
    "span85": SurfactantProfile(
        name="Span-85 (sorbitan trioleate)",
        cas="26266-58-0",
        hlb=1.8,
        mw=957.5,
        # Three oleate chains: larger molecule, packs less tightly at interface.
        # Gamma_inf lower due to larger molecular footprint (~120 A2 vs 47 A2).
        gamma_inf=1.8e-6,  # [estimated from molecular area]
        # K_L high — very oil-soluble, strong partitioning to interface.
        K_L=1.5,  # [estimated]
        sigma_0_paraffin=0.050,
        suitability=6,
        notes=(
            "Tri-ester sorbitan: very lipophilic (HLB 1.8). Excellent oil "
            "solubility, liquid at RT. Produces highly stable W/O emulsions "
            "resistant to coalescence. Larger molecular footprint means fewer "
            "molecules per unit area — potentially larger droplets at equal "
            "mass concentration. Very viscous (300-400 mPa.s at 25 degC), "
            "which increases oil phase viscosity and may affect breakage "
            "kernels. Used as vaccine adjuvant emulsifier (MF59). "
            "Compatible with agarose gelation. Cost: ~USD 0.15/g."
        ),
    ),
    # ── 5. PGPR (polyglycerol polyricinoleate) ────────────────────────────
    "pgpr": SurfactantProfile(
        name="PGPR (polyglycerol polyricinoleate)",
        cas="29894-35-7",
        hlb=1.5,
        # MW is a distribution; typical ~3000 g/mol for commercial PGPR.
        mw=3000.0,
        # PGPR forms thick viscoelastic interfacial films.
        # Gamma_inf is high despite large MW because the molecule unfolds
        # at the interface. Wilson et al. (1998): Gamma ~ 2-3 umol/m2.
        gamma_inf=2.5e-6,  # [estimated from Wilson et al.]
        # K_L very high — PGPR adsorbs strongly and nearly irreversibly.
        K_L=3.0,  # [estimated]
        sigma_0_paraffin=0.050,
        suitability=7,
        notes=(
            "Polymeric emulsifier widely used in chocolate and food W/O "
            "emulsions. Excellent at very low HLB applications. Forms thick, "
            "viscoelastic interfacial films that resist coalescence even "
            "under shear. Can produce smaller droplets than Span-80 at "
            "equal mass loading. PGPR is polydisperse (MW 1000-5000) — "
            "batch variability is a concern for reproducible chromatography "
            "bead production. Food-grade (E476), low toxicity. Residual "
            "PGPR on beads is difficult to wash off due to strong adsorption "
            "— may affect chromatographic performance (non-specific binding). "
            "Szyszkowski model is approximate for polymeric surfactants; "
            "a Frumkin or multilayer model would be more accurate. "
            "Cost: ~USD 0.20/g."
        ),
    ),
    # ── 6. Lecithin (phosphatidylcholine) ─────────────────────────────────
    "lecithin": SurfactantProfile(
        name="Lecithin (soy phosphatidylcholine)",
        cas="8002-43-5",
        hlb=4.0,
        # Lecithin is a mixture; predominant species is dipalmitoyl-PC.
        # MW ~ 760 g/mol for DPPC.
        mw=760.0,
        # Phospholipid monolayer: well-characterized packing.
        # Gamma_inf ~ 4.0e-6 mol/m2 from Langmuir trough measurements
        # (molecular area ~40 A2 at collapse, Peltonen & Yliruusi 2000).
        gamma_inf=4.0e-6,
        # K_L moderate — lecithin partitions between oil bulk (reverse
        # micelles) and interface. Shrestha et al. (2010).
        K_L=0.60,  # [estimated]
        sigma_0_paraffin=0.050,
        suitability=6,
        notes=(
            "Natural phospholipid mixture from soy or egg. Zwitterionic "
            "headgroup provides good IFT reduction. Forms lamellar and "
            "reverse micellar structures in oil. HLB ~4 suitable for W/O. "
            "Biocompatible and biodegradable. CONCERNS: lecithin is a "
            "mixture (PC, PE, PI, PA) — batch variability. Can form "
            "liposomal/lamellar structures at the interface that may "
            "complicate droplet breakage kinetics. At pH < 4 or > 9, "
            "headgroup ionization changes and emulsion stability decreases. "
            "May interact with chitosan via electrostatic attraction "
            "(negative PA/PI components + positive chitosan NH3+), which "
            "could either stabilize or destabilize the interface depending "
            "on stoichiometry. Compatible with agarose gelation. "
            "Cost: ~USD 0.05/g."
        ),
    ),
}


# ═══════════════════════════════════════════════════════════════════════════
#  QUICK REFERENCE: SIMULATION PARAMETER MAPPING
# ═══════════════════════════════════════════════════════════════════════════
#
#  For each crosslinker, the Level 3 solver needs:
#    MaterialProperties.k_xlink_0   <-- CrosslinkerProfile.k_xlink_0
#    MaterialProperties.E_a_xlink   <-- CrosslinkerProfile.E_a_xlink
#    MaterialProperties.f_bridge    <-- CrosslinkerProfile.f_bridge
#    FormulationParameters.T_crosslink  <-- CrosslinkerProfile.T_crosslink_default
#    FormulationParameters.t_crosslink  <-- CrosslinkerProfile.t_crosslink_default
#
#  For kinetics_model != "second_order", the solver.py ODE RHS must be
#  extended:
#    - "michaelis_menten":  dX/dt = Vmax*[G]*[NH2] / (Km + [NH2])
#    - "uv_dose":           dX/dt = k0*sqrt(I_UV)*exp(-mu*r)*[PEGDA]*[radical]
#    - "ionic_instant":     X = min(n_TPP * valence, [NH3+]) (equilibrium, no ODE)
#
#  For each surfactant, the Level 1 interfacial tension model needs:
#    gamma_inf   <-- SurfactantProfile.gamma_inf
#    K_L         <-- SurfactantProfile.K_L
#    sigma_0     <-- SurfactantProfile.sigma_0_paraffin  (adjusted for T)
#    M_surfactant (MW for c_mass -> c_mol conversion) <-- SurfactantProfile.mw
# ═══════════════════════════════════════════════════════════════════════════
