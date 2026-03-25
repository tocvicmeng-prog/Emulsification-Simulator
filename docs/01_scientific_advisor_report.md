

# Comprehensive Scientific Decomposition: Emulsification and Double-Network Hydrogel Microsphere Formation

---

## 1. FUNDAMENTAL PHYSICS AND CHEMISTRY DECOMPOSITION

### 1A. Emulsification / Droplet Formation

The emulsification step is a W/O (water-in-oil) emulsion process in which the hot aqueous polysaccharide solution is dispersed into paraffin oil under intense mechanical agitation. The fundamental phenomena are:

#### 1A.1 Droplet Breakup Mechanisms

Two regimes govern droplet breakup depending on the ratio of droplet size to the Kolmogorov microscale $\eta_K = (\nu^3/\varepsilon)^{1/4}$:

**Turbulent Inertial Breakup (droplet >> $\eta_K$):**
The droplet is larger than the smallest eddies. Breakup occurs when the turbulent pressure fluctuations across the droplet overcome the restoring capillary pressure. The controlling dimensionless group is the **Weber number**:

$$We = \frac{\rho_c \, \overline{u'^2(d)} \, d}{\sigma}$$

where $\rho_c$ is the continuous phase (oil) density, $\overline{u'^2(d)}$ is the mean-square velocity fluctuation over scale $d$, and $\sigma$ is the interfacial tension. Breakup occurs when $We > We_{crit} \approx 1$ (order unity). Using Kolmogorov's second-order structure function scaling $\overline{u'^2(d)} \sim (\varepsilon d)^{2/3}$, this yields the classical Kolmogorov-Hinze prediction for maximum stable droplet size:

$$d_{max} = C_1 \left(\frac{\sigma}{\rho_c}\right)^{3/5} \varepsilon^{-2/5}$$

where $C_1 \approx 0.725$ (Hinze 1955), and $\varepsilon$ is the turbulent energy dissipation rate per unit mass ($\text{m}^2/\text{s}^3$).

**Turbulent Viscous Breakup (droplet << $\eta_K$):**
When the droplet is smaller than the Kolmogorov scale, it sits within the viscous subrange and experiences a locally laminar shear. The controlling group is the **Capillary number**:

$$Ca = \frac{\mu_c \dot{\gamma} \, d}{2\sigma}$$

where $\mu_c$ is the continuous phase viscosity and $\dot{\gamma}$ is the local shear rate. Critical capillary number depends on the viscosity ratio $\lambda = \mu_d/\mu_c$ (Taylor 1934, Grace 1982). For $\lambda \sim O(1)$, $Ca_{crit} \approx 0.5$. The maximum stable droplet size in this regime:

$$d_{max} = C_2 \, \sigma \, \mu_c^{-1} \, \varepsilon^{-1/2} \, \nu_c^{1/2}$$

**Critical assessment for target $d \sim 2\,\mu m$:** At 2 $\mu$m target size, the droplets will almost certainly be in the viscous sub-Kolmogorov regime for any reasonable stirring configuration (a homogenizing mixer at ~5000-20000 RPM in a typical vessel gives $\eta_K \sim 10$-$50\,\mu$m). This is a crucial design consideration: the Kolmogorov-Hinze equation in its classical inertial form is **not directly applicable**. One must use the viscous breakup correction or the modified Hinze correlation that accounts for dispersed-phase viscosity:

$$d_{max} = C_1 \left(\frac{\sigma}{\rho_c}\right)^{3/5} \varepsilon^{-2/5} \left(1 + C_3 \, \frac{\mu_d}{\sqrt{\rho_c \sigma \, d_{max}}}\right)$$

(Davies 1985 modification, or the Calabrese-Chang-Dang 1986 model). The dispersed phase viscosity is very high because the agarose-chitosan solution at 6% w/v is extremely viscous ($\mu_d \gg \mu_c$), which strongly resists deformation and breakup. This viscosity ratio effect is **the dominant challenge** for achieving 2 $\mu$m droplets.

#### 1A.2 Interfacial Tension

The interfacial tension $\sigma$ between the aqueous polysaccharide solution and paraffin oil, modified by Span-80, is the primary restoring force resisting droplet deformation.

- **Bare interface** (no surfactant): $\sigma_0$ for water/liquid paraffin $\approx 40$-$50$ mN/m at 20Â°C. At 90Â°C, expect reduction to $\sim 30$-$40$ mN/m due to thermal effects ($d\sigma/dT \approx -0.1$ mN/m/K for most hydrocarbon-water systems).
- **With dissolved polysaccharides**: Agarose is mildly surface-active; chitosan is more so (protonated amines act as weak surfactant). Expect modest reduction of $\sigma_0$ by 5-15%.
- **With Span-80**: At concentrations above CMC, Span-80 dramatically reduces $\sigma$ to $\sim 2$-$8$ mN/m. Span-80 (sorbitan monooleate, HLB = 4.3) is a low-HLB surfactant ideal for stabilizing W/O emulsions.

The time-dependent interfacial tension during emulsification follows the Ward-Tordai equation for diffusion-controlled adsorption:

$$\Gamma(t) = 2 C_0 \sqrt{\frac{D t}{\pi}} - 2\sqrt{\frac{D}{\pi}} \int_0^{\sqrt{t}} C_s(\tau) \, d(\sqrt{t - \tau})$$

where $\Gamma$ is surface excess, $C_0$ is bulk Span-80 concentration, $D$ is Span-80 diffusion coefficient in oil ($\sim 10^{-10}$ m$^2$/s), and $C_s$ is subsurface concentration. This is relevant because at very short times during rapid breakup, the interface may not be fully covered, giving a **dynamic interfacial tension** $\sigma_{dyn} > \sigma_{eq}$, which makes breakup harder than equilibrium predictions suggest.

#### 1A.3 Role of Span-80

- **HLB**: 4.3, placing it firmly in the W/O emulsifier category (Bancroft's rule: the phase in which the surfactant is more soluble becomes the continuous phase).
- **CMC in paraffin oil**: $\sim 0.5$-$2\%$ w/v depending on temperature and oil grade. Below CMC, interfacial tension drops roughly logarithmically with concentration (Gibbs adsorption isotherm: $\Gamma = -\frac{1}{RT}\frac{d\sigma}{d\ln C}$). Above CMC, $\sigma$ plateaus.
- **Adsorption kinetics**: Span-80 is oil-soluble and adsorbs from the continuous phase. Adsorption timescale $\tau_{ads} \sim \Gamma_{eq}^2 / (D C_0^2)$. For typical concentrations (2-5% v/v), $\tau_{ads} \sim O(1\,\text{ms})$, which is comparable to the droplet breakup timescale in a homogenizer. This means **Span-80 adsorption may be rate-limiting** for stabilization of freshly-formed droplets, leading to re-coalescence of insufficiently stabilized drops.
- **Steric/mechanical barrier**: Span-80 provides a mechanical film at the interface that resists coalescence via the Gibbs-Marangoni effect (local thinning of the film draws surfactant by surface tension gradients, restoring the film).

#### 1A.4 Population Balance Equation (PBE)

The droplet size distribution $n(v,t)$ (number density of droplets of volume $v$ at time $t$) evolves according to:

$$\frac{\partial n(v,t)}{\partial t} = \underbrace{\int_v^{\infty} \beta(v,v') \, g(v') \, n(v',t) \, dv'}_{\text{birth by breakage}} - \underbrace{g(v) \, n(v,t)}_{\text{death by breakage}} + \underbrace{\frac{1}{2}\int_0^v h(v-v',v') \, n(v-v',t) \, n(v',t) \, dv'}_{\text{birth by coalescence}} - \underbrace{n(v,t)\int_0^{\infty} h(v,v') \, n(v',t) \, dv'}_{\text{death by coalescence}}$$

(Ramkrishna 2000)

where:
- $g(v)$ = breakage frequency (rate at which a droplet of volume $v$ breaks), $[\text{s}^{-1}]$
- $\beta(v,v')$ = daughter-size distribution (probability density that breakage of $v'$ produces a fragment of size $v$)
- $h(v,v')$ = coalescence frequency between droplets of volume $v$ and $v'$, $[\text{m}^3/\text{s}]$

**Breakage frequency** â€” Coulaloglou & Tavlarides (1977):

$$g(d) = \frac{C_1 \, \varepsilon^{1/3}}{d^{2/3}} \exp\left(-\frac{C_2 \, \sigma}{\rho_c \, \varepsilon^{2/3} \, d^{5/3}}\right)$$

This model assumes inertial-regime breakup. For the viscous regime relevant to 2 $\mu$m targets, the Luo & Svendsen (1996) model is more rigorous:

$$g(d) = 0.9238 \, \varepsilon^{1/3} \, d^{-2/3} \int_{\xi_{min}}^{1} \frac{(1+\xi)^2}{\xi^{11/3}} \exp\left(-\frac{12 c_f \sigma}{\beta \rho_c \varepsilon^{2/3} d^{5/3} \xi^{11/3}}\right) d\xi$$

where $\xi = \lambda/d$ is the ratio of eddy size to droplet size, $c_f$ is the increase in surface energy, and $\beta = 2.045$.

For high-viscosity dispersed phase, the Alopaeus et al. (2002) modification is recommended, which adds a viscous resistance term:

$$g(d) = C_1 \varepsilon^{1/3} d^{-2/3} \exp\left(-C_2 \frac{\sigma}{\rho_c \varepsilon^{2/3} d^{5/3}} - C_3 \frac{\mu_d}{\rho_c^{1/2} \varepsilon^{1/3} d^{4/3} \sigma^{1/2}}\right)$$

**Coalescence frequency** â€” Coulaloglou & Tavlarides (1977):

$$h(d_i, d_j) = C_4 \, \varepsilon^{1/3} (d_i^2 + d_j^2)(d_i^{2/3} + d_j^{2/3})^{1/2} \exp\left(-\frac{C_5 \, \mu_c \, \rho_c \, \varepsilon}{\sigma^2} \left(\frac{d_i d_j}{d_i + d_j}\right)^4\right)$$

The exponential term represents the film drainage efficiency â€” high surfactant concentration increases the drainage resistance and effectively suppresses coalescence.

---

### 1B. Double-Network Hydrogel Formation

#### 1B.1 Agarose Gelation

Agarose gelation is a **physical gelation** driven by a thermally-reversible coil-to-helix conformational transition followed by helix aggregation.

**Molecular mechanism** (Arnott et al. 1974; Djabourov et al. 1989):
1. At $T > T_{sol} \approx 85$-$90Â°$C: random coil conformation in solution
2. Upon cooling below $T_{gel} \approx 35$-$45Â°$C (varies with MW, concentration): chains undergo coil â†’ double-helix transition (left-handed 3-fold helix)
3. Helices laterally aggregate into bundles (junction zones), forming a fibrillar network
4. The large thermal hysteresis ($T_{sol} - T_{gel} \approx 40$-$50Â°$C) is due to the cooperative nature of helix bundle formation

**Gelation temperature** depends on agarose type and concentration:

$$T_{gel} \approx T_{gel,0} + k_c \cdot c$$

where $c$ is polymer concentration and $k_c \approx 0.5$-$2$ Â°C/(% w/v).

**Gelation kinetics** â€” Aymard et al. (2001) isothermal model:

$$\alpha(t) = 1 - \exp\left(-(t/\tau_{gel})^n\right)$$

where $\alpha$ is the fraction of gel formed (0 = sol, 1 = fully gelled), $\tau_{gel}$ is a characteristic gelation time that depends strongly on temperature (Arrhenius-like), and $n$ is an Avrami exponent ($n \approx 2$-$3$, indicating nucleation-and-growth kinetics).

The temperature dependence of the gelation rate:

$$\frac{1}{\tau_{gel}(T)} = A \exp\left(-\frac{E_a}{R}\left(\frac{1}{T} - \frac{1}{T_{gel}}\right)\right) \quad \text{for } T < T_{gel}$$

with $E_a \sim 200$-$400$ kJ/mol for agarose (high cooperativity).

**Rheological signature**: Storage modulus $G'$ during gelation follows:

$$G'(t) = G'_{\infty} \left[1 - \exp\left(-(t/\tau)^n\right)\right]$$

with $G'_{\infty}$ for 6% agarose $\approx 10$-$50$ kPa depending on MW.

#### 1B.2 Chitosan Behavior

Chitosan (degree of deacetylation > 90%) is a cationic polysaccharide soluble in dilute acid (pH < 6.2) due to protonation of the primary amine ($-NH_2 + H^+ \to -NH_3^+$, $pK_a \approx 6.3$-$6.5$).

**Critical consideration**: At neutral/alkaline pH, chitosan precipitates. In the aqueous phase preparation, if the chitosan is dissolved in acidic solution and mixed with hot agarose (prepared in water), the pH of the mixture determines whether chitosan remains soluble during emulsification. If the mixture is near-neutral pH, chitosan may begin to phase-separate from the agarose solution even before gelation, creating a **spinodal or nucleation-and-growth microstructure** that templates the pore structure.

**Chitosan "gelation"**: Unlike agarose, chitosan does not undergo a thermoreversible physical gelation. Instead, it forms a gel by:
1. **pH-induced gelation**: Raising pH above $pK_a$ causes deprotonation â†’ hydrophobic interactions â†’ physical gel
2. **Chemical crosslinking** (genipin, glutaraldehyde, etc.)
3. **Ionic crosslinking** (tripolyphosphate, etc.)

In this system, the relevant mechanism is primarily pH-induced chain association during or after the emulsification step, followed by genipin chemical crosslinking.

#### 1B.3 IPN vs Semi-IPN Architecture

A true **interpenetrating polymer network (IPN)** requires both networks to be independently crosslinked and mutually interlocked. A **semi-IPN** has only one network crosslinked, with the second polymer physically entangled.

In this system:
- **First network**: Agarose physical gel (helix-bundle junction zones) â€” formed upon cooling
- **Second network**: Chitosan crosslinked by genipin â€” formed during post-gelation crosslinking step

This is technically a **sequential IPN** (the agarose network forms first during cooling, then the chitosan network is locked in by genipin crosslinking). However, if the chitosan is not independently crosslinked into a continuous percolating network but merely trapped within the agarose matrix, it is a semi-IPN.

The **double-network (DN)** concept (Gong et al. 2003) typically refers to a stiff, brittle first network and a soft, ductile second network, where the first network sacrificially breaks to dissipate energy. In the agarose-chitosan system, agarose is the rigid first network and chitosan is the flexible second network â€” consistent with DN mechanics.

#### 1B.4 Phase Separation During Cooling

This is a **key phenomenon** for pore formation. During cooling, as the agarose undergoes coil-to-helix transition, two types of phase separation can occur:

1. **Polymer-solvent phase separation**: As agarose chains aggregate, they exclude water, creating polymer-rich (fibrillar) and polymer-poor (aqueous pore) domains. This is the primary mechanism of macropore formation in agarose gels.

2. **Polymer-polymer phase separation**: Agarose and chitosan are thermodynamically incompatible (different backbone chemistry, charge states). The Flory-Huggins interaction parameter $\chi_{12}$ between agarose and chitosan determines the miscibility:

$$\Delta G_{mix} = RT\left[\frac{\phi_1}{N_1}\ln\phi_1 + \frac{\phi_2}{N_2}\ln\phi_2 + \chi_{12}\phi_1\phi_2\right]$$

where $\phi_i$ are volume fractions and $N_i$ are degrees of polymerization. When $\chi_{12}$ exceeds a critical value, the mixture undergoes spinodal decomposition, creating a bicontinuous or droplet morphology that is frozen in place by the agarose gelation.

**This polymer-polymer phase separation is likely the dominant mechanism for achieving the 60-100 nm pore structure** in the double-network system â€” the chitosan-rich domains template the pores within the agarose matrix, or vice versa.

---

### 1C. Pore Structure Formation

#### 1C.1 Spinodal Decomposition vs Nucleation-and-Growth

The morphology depends on the cooling trajectory through the phase diagram:

- **Spinodal decomposition** (quench into the unstable region): produces a bicontinuous, interconnected pore network with a characteristic wavelength $\lambda_c$ given by:

$$\lambda_c = 2\pi \sqrt{\frac{-2\kappa}{\partial^2 f/\partial \phi^2}}$$

where $\kappa$ is the gradient energy coefficient and $f$ is the free energy density. This yields **interconnected, wormlike** pore morphology â€” desirable for chromatographic applications.

- **Nucleation-and-growth** (quench into the metastable region): produces discrete, spherical pores â€” less interconnected.

For agarose at 6% w/v, the gel typically forms with a fibrillar network structure where fiber diameter $\sim 10$-$30$ nm and pore size depends on concentration. The agarose fiber network can be described by the Ogston model:

$$K_{av} = \exp\left(-\pi r_s^2 / (\pi r_f^2 \cdot n_f)\right)$$

or more directly, the average pore radius from the Ogston-Morris-Rodbard-Chrambach theory:

$$\langle r_{pore} \rangle \approx r_f \sqrt{\frac{\pi}{4\phi_p}} - r_f$$

where $r_f$ is the fiber radius and $\phi_p$ is the polymer volume fraction.

#### 1C.2 Role of Cooling Rate

**Fast cooling** (quench): 
- Deeper penetration into the unstable region â†’ spinodal decomposition dominates
- Many nucleation events â†’ smaller domain size â†’ smaller pores
- Less time for coarsening â†’ finer structure preserved

**Slow cooling** (gradual):
- Passes through metastable region â†’ nucleation-and-growth may dominate initially
- More time for Ostwald ripening / coarsening â†’ larger pores
- Coarsening follows $\lambda(t) \sim t^{1/3}$ (Lifshitz-Slyozov-Wagner)

**For the target 60-100 nm pores**: This is relatively fine compared to typical agarose gels (which naturally have pores of $\sim 100$-$500$ nm at 2-6% depending on type). Achieving 60-100 nm consistently requires controlled rapid cooling, high polymer concentration, and/or the presence of chitosan to modify the phase separation length scale.

#### 1C.3 Mesh Size Models

**Pore model from the paper** (Zhao et al. 2020):

$$\xi = 0.071 \cdot \nu_{2,s}^{-1/3} \cdot M_c^{1/2}$$

where:
- $\xi$ = mesh size (m)
- $\nu_{2,s}$ = polymer volume fraction in the swollen state
- $M_c$ = molecular weight between crosslinks (g/mol)

This is derived from **rubber elasticity theory combined with Flory-Rehner**. The general form (Canal & Peppas 1989) is:

$$\xi = \nu_{2,s}^{-1/3} \cdot \left(\bar{l}^2 \cdot C_n \cdot \frac{2M_c}{M_r}\right)^{1/2}$$

where $\bar{l}$ is the bond length along the backbone, $C_n$ is the characteristic ratio (backbone stiffness), and $M_r$ is the repeat unit molecular weight. The coefficient 0.071 in the Zhao equation bundles $\bar{l}$, $C_n$, and $M_r$ for agarose.

**Flory-Rehner theory** for swelling equilibrium:

$$\ln(1-\nu_{2,s}) + \nu_{2,s} + \chi \nu_{2,s}^2 = -V_1 \cdot n_e \left(\nu_{2,s}^{1/3} - \frac{\nu_{2,s}}{2}\right)$$

where $V_1$ is the molar volume of solvent, $n_e$ is the effective crosslink density (mol/m$^3$), and $\chi$ is the Flory-Huggins interaction parameter.

---

### 1D. Crosslinking (Genipin)

#### 1D.1 Reaction Mechanism

Genipin reacts with **primary amines** on chitosan via a two-step mechanism (Butler et al. 2003; Muzzarelli 2009):

1. **Nucleophilic attack**: The C-3 atom of the genipin dihydropyran ring undergoes nucleophilic attack by the chitosan $-NH_2$ group, opening the ring to form a secondary amine linkage.
2. **Bridge formation**: The ester group of genipin undergoes amide bond formation with a second amine, creating a **crosslink bridge** between two chitosan chains.

The crosslinking creates covalent bonds ($C-N$ bonds, bond energy $\sim 305$ kJ/mol), which are much stronger and more stable than the physical crosslinks in the agarose network.

#### 1D.2 Crosslinking Kinetics

The genipin-chitosan crosslinking follows second-order kinetics (Dimida et al. 2017):

$$\frac{d[\text{crosslinks}]}{dt} = k_{xlink}(T) \cdot [\text{NH}_2]_{free} \cdot [\text{genipin}]_{free}$$

where:

$$k_{xlink}(T) = A_{xlink} \exp\left(-\frac{E_{a,xlink}}{RT}\right)$$

Literature values: $E_{a,xlink} \approx 40$-$60$ kJ/mol; at 37Â°C, $k_{xlink} \sim 10^{-3}$-$10^{-2}$ L/(molÂ·s) depending on pH, degree of deacetylation, and solvent.

The **degree of crosslinking** $\alpha_{xlink}$ evolves as:

$$\alpha_{xlink}(t) = \frac{[\text{crosslinks}](t)}{[\text{NH}_2]_{0}/2}$$

where the factor of 2 accounts for the fact that each crosslink consumes two amine groups.

#### 1D.3 Effect on Mesh Size

As crosslinking proceeds, $M_c$ decreases (shorter chains between crosslinks):

$$M_c = \frac{M_0}{\alpha_{xlink} \cdot f/2}$$

where $M_0$ is the molecular weight of the uncrosslinked polymer and $f$ is the functionality of the crosslinker ($f = 2$ for genipin bridge). This causes the mesh size $\xi$ to decrease (per the Zhao equation). However, the Zhao paper reports that **increased crosslinking leads to larger pores** â€” this apparent contradiction is resolved by recognizing that in the agarose system, chemical crosslinking by epichlorohydrin can cause local chain contraction that opens up macropores (different from the mesh size of the gel phase itself). There are two length scales: the **macropore** (visible in SEM, 60-100+ nm) and the **gel-phase mesh** ($\xi$, typically 1-10 nm).

---

### 1E. Mechanical Properties

#### 1E.1 Rubber Elasticity

For a chemically crosslinked hydrogel, the shear modulus is:

$$G = n_e k_B T = \frac{\rho_{polymer} R T}{M_c} \cdot \nu_{2,s}^{1/3}$$

(Flory 1953, affine network model)

where $n_e$ is the effective crosslink density (mol/m$^3$).

The phantom network model (James & Guth 1947) gives:

$$G = \left(1 - \frac{2}{f}\right) n_e k_B T$$

For a tetrafunctional network ($f=4$), $G_{phantom} = 0.5 \, G_{affine}$.

#### 1E.2 Double-Network Mechanics

The double-network concept (Gong 2010) gives enhanced toughness through sacrificial bond breaking. The simplest composite model (assuming parallel coupling of the two networks):

$$G_{DN} = G_{agarose} + G_{chitosan-genipin}$$

More realistically, for a sequential IPN where the networks are not perfectly independent, the Lake-Thomas theory gives the fracture energy:

$$\Gamma = \frac{1}{2} n_{chains} \cdot N_{monomers} \cdot U_{bond}$$

where $n_{chains}$ is the areal density of first-network chains crossing the fracture plane, $N_{monomers}$ is the number of monomers per strand, and $U_{bond}$ is the bond dissociation energy. This is relevant for predicting the mechanical robustness of the microspheres under chromatographic column packing pressures.

#### 1E.3 Compression of Microspheres

For a microsphere under uniaxial compression between two plates, the Hertz contact model gives:

$$F = \frac{4}{3} E^* R^{1/2} \delta^{3/2}$$

where $E^* = E/(1-\nu^2)$ is the reduced modulus, $R$ is the microsphere radius, and $\delta$ is the indentation depth. For hydrogels, $\nu \approx 0.45$-$0.5$ (nearly incompressible). The Young's modulus relates to shear modulus as $E = 2G(1+\nu) \approx 3G$.

---

## 2. GOVERNING EQUATIONS AND MATHEMATICAL FRAMEWORKS

### 2.1 Emulsification â€” Population Balance Equation

**Equation:**

$$\frac{\partial n(v,t)}{\partial t} + \nabla \cdot [\mathbf{u}(\mathbf{x},t) \, n(v,\mathbf{x},t)] = B_{break} - D_{break} + B_{coal} - D_{coal}$$

The left-hand side includes a convective transport term when coupled to CFD (spatially inhomogeneous turbulence).

**Variables:**
| Symbol | Meaning | Units |
|--------|---------|-------|
| $n(v,t)$ | Number density of droplets of volume $v$ | m$^{-6}$ |
| $g(v)$ | Breakage frequency | s$^{-1}$ |
| $\beta(v,v')$ | Daughter size distribution | m$^{-3}$ |
| $h(v,v')$ | Coalescence frequency | m$^3$ s$^{-1}$ |
| $\varepsilon$ | Turbulent energy dissipation rate | m$^2$ s$^{-3}$ |
| $\sigma$ | Interfacial tension | N/m |

**Source:** Ramkrishna, D. (2000) *Population Balances: Theory and Applications to Particulate Systems in Engineering*, Academic Press.

**Validity/Assumptions:** Assumes well-mixed system (0D) or resolved spatial distribution (CFD-coupled). Droplets are spherical. No mass transfer between phases during breakup. Kernels assume isotropic turbulence.

**Parameter availability:** $C_1$-$C_5$ in C&T model must be fitted to experimental data for this specific system (oil type, surfactant, polymer solution); no direct literature values exist for agarose-chitosan/paraffin/Span-80.

### 2.2 Kolmogorov-Hinze Equation (with viscosity correction)

**Equation (Calabrese et al. 1986):**

$$d_{max} = C_1 \left(\frac{\sigma}{\rho_c}\right)^{3/5} \varepsilon^{-2/5} \left[1 + C_3 \cdot Vi\right]^{3/5}$$

where the viscosity group $Vi = \mu_d / \sqrt{\rho_c \sigma d_{max}}$.

**Source:** Calabrese, R.V., Chang, T.P.K., Dang, P.T. (1986) Drop breakup in turbulent stirred-tank contactors. Part I. *AIChE J.* 32:657-666.

**Validity:** Applicable for $d > \eta_K$ (inertial subrange). For $d < \eta_K$, must switch to viscous-subrange model. The iterative nature (implicit in $d_{max}$) requires numerical solution.

**Parameter availability:**
- $\rho_c$ (paraffin oil): $\sim 830$-$870$ kg/m$^3$ at 20Â°C â€” **direct** from manufacturer data
- $\mu_c$ (paraffin oil at 90Â°C): $\sim 3$-$10$ mPaÂ·s â€” **direct** from manufacturer data
- $\mu_d$ (6% agarose + chitosan at 90Â°C): $\sim 10$-$500$ mPaÂ·s â€” **indirect**, must be measured or estimated from Mark-Houwink
- $\sigma$ (with Span-80): $\sim 2$-$8$ mN/m â€” **indirect**, must be measured or estimated from Span-80 isotherm
- $\varepsilon$: **indirect**, depends on mixer geometry and RPM (see Section 2.9)

### 2.3 Agarose Gelation Kinetics

**Equation (modified Avrami/Aymard):**

$$\frac{d\alpha}{dt} = n \cdot k_{gel}(T) \cdot (1-\alpha) \left[-\ln(1-\alpha)\right]^{(n-1)/n}$$

This is the differential form of the Avrami equation $\alpha(t) = 1 - \exp(-(k_{gel} t)^n)$.

**Temperature dependence:**

$$k_{gel}(T) = k_0 \exp\left(-\frac{E_a}{R T}\right) \cdot \Theta(T_{gel} - T)$$

where $\Theta$ is a Heaviside function ensuring gelation only occurs below $T_{gel}$.

**Source:** Aymard, P., Martin, D.R., Plucknett, K., Foster, T.J., Clark, A.H., Norton, I.T. (2001) Influence of thermal history on the structural and mechanical properties of agarose gels. *Biopolymers* 59:131-144.

**Parameters:**
- $n \approx 2$-$3$: **direct** from Aymard et al.
- $k_0, E_a$: **indirect**, must be fitted from oscillatory rheology cooling ramps specific to this agarose grade + chitosan mixture
- $T_{gel}$: $\approx 36$-$42Â°$C for standard agarose â€” **direct** from DSC or rheology; shifts with chitosan addition â€” **indirect**

### 2.4 Mark-Houwink Equation for Agarose

**Equation (from Zhao et al.):**

$$[\eta] = 0.07 \cdot M^{0.72}$$

where $[\eta]$ is intrinsic viscosity (dL/g) and $M$ is molecular weight (g/mol).

**Practical viscosity at concentration $c$** â€” Huggins equation:

$$\eta_{sp}/c = [\eta] + k_H [\eta]^2 c$$

or for concentrated solutions, the Martin equation:

$$\ln(\eta_{sp}/c) = \ln[\eta] + k_M [\eta] c$$

For 6% w/v agarose solution at 90Â°C, the zero-shear viscosity can be estimated from these, but the agarose-chitosan blend viscosity requires additional mixing rules (e.g., logarithmic mixing: $\ln \eta_{mix} = \phi_1 \ln \eta_1 + \phi_2 \ln \eta_2$).

**Source:** Zhao et al. (2020), *Eng. Life Sci.* 20:504-513.

### 2.5 Flory-Rehner Swelling Theory

**Equation:**

$$-[\ln(1-\nu_{2,s}) + \nu_{2,s} + \chi_1 \nu_{2,s}^2] = V_1 n_e \left[\nu_{2,s}^{1/3} - \frac{\nu_{2,s}}{2}\right]$$

Rearranged for crosslink density:

$$n_e = \frac{-[\ln(1-\nu_{2,s}) + \nu_{2,s} + \chi_1 \nu_{2,s}^2]}{V_1 [\nu_{2,s}^{1/3} - \nu_{2,s}/2]}$$

And $M_c = \rho_{polymer} / n_e$.

**Variables:**
| Symbol | Meaning | Typical Value |
|--------|---------|---------------|
| $\nu_{2,s}$ | Polymer volume fraction in swollen gel | 0.02-0.10 for 2-10% gels |
| $\chi_1$ | Flory-Huggins interaction parameter (polymer-solvent) | 0.47-0.50 for agarose-water |
| $V_1$ | Molar volume of water | $1.8 \times 10^{-5}$ m$^3$/mol |

**Source:** Flory, P.J., Rehner, J. (1943) Statistical mechanics of cross-linked polymer networks. *J. Chem. Phys.* 11:521-526.

**Parameter availability:** $\chi_1$ for agarose-water is $\approx 0.497$ (**direct**, Normand et al. 2000). For chitosan-water, $\chi_1 \approx 0.4$-$0.5$ depending on DD and pH. For the blend, $\chi_{eff}$ is **indirect** and must be estimated.

### 2.6 Genipin Crosslinking Kinetics

**Equation:**

$$\frac{d[X]}{dt} = k(T) \cdot [NH_2]_{free}(t) \cdot [Gen]_{free}(t)$$

$$[NH_2]_{free}(t) = [NH_2]_0 - 2[X](t)$$

$$[Gen]_{free}(t) = [Gen]_0 - [X](t)$$

where $[X]$ is crosslink concentration. This is a standard bimolecular kinetics ODE.

**Rate constant:** $k(37Â°C, pH 7.4) \approx 5 \times 10^{-3}$ L/(molÂ·s) for chitosan (DD > 85%) â€” **indirect**, compiled from Dimida et al. (2017) and Mi et al. (2005), but sensitive to pH, DD, and solvent.

**Source:** Dimida, S., et al. (2017) Genipin-cross-linked chitosan-based hydrogels: Reaction kinetics and structure-related characteristics. *J. Appl. Polym. Sci.* 134:45006.

### 2.7 Heat Transfer During Cooling

**For a single microsphere** (Biot number check first):

$$Bi = \frac{h_{conv} \cdot R}{k_{gel}} \approx \frac{(100\text{-}500)(1 \times 10^{-6})}{0.5} \sim 10^{-4}\text{-}10^{-3}$$

Since $Bi \ll 0.1$, the **lumped capacitance model** is valid â€” temperature is uniform within each microsphere:

$$\rho_{drop} c_p V \frac{dT}{dt} = -h_{conv} A (T - T_{oil}(t))$$

$$T(t) = T_{oil}(t) + [T_0 - T_{oil,0}] \exp\left(-\frac{t}{\tau_{cool}}\right)$$

where $\tau_{cool} = \frac{\rho c_p R}{3 h_{conv}}$. For $R = 1\,\mu$m, $\tau_{cool} \sim 10^{-6}$ s â€” effectively instantaneous equilibration with the oil phase. **Therefore, the cooling rate of the microsphere is entirely controlled by the bulk oil cooling rate**, and the microsphere tracks the oil temperature with negligible lag.

**For the bulk emulsion** (oil bath cooling):

$$\rho_{eff} c_{p,eff} V_{vessel} \frac{dT_{oil}}{dt} = -U A_{vessel} (T_{oil} - T_{ambient}) + Q_{stirrer}$$

where $U$ is the overall heat transfer coefficient of the vessel wall, and $Q_{stirrer}$ is the viscous dissipation from the mixer (typically negligible).

**Source:** Standard heat transfer; Incropera & DeWitt (2002), *Fundamentals of Heat and Mass Transfer*, Wiley.

### 2.8 Viscosity of the Agarose-Chitosan Solution at Emulsification Temperature

At 90Â°C, the solution is well above $T_{gel}$, so the viscosity is purely hydrodynamic (no gel contribution). Using the Mark-Houwink relation for agarose:

$$[\eta]_{agarose} = 0.07 \cdot M_w^{0.72}$$

For $M_w \approx 120{,}000$ g/mol (typical agarose): $[\eta] \approx 0.07 \times 120000^{0.72} \approx 0.07 \times 5600 \approx 390$ mL/g. At $c = 0.06$ g/mL ($6\%$ w/v), $c[\eta] \approx 23$, which is well above the entanglement concentration $c^* \approx 1/[\eta] \approx 0.0026$ g/mL. The solution is in the **concentrated entangled regime**, and zero-shear viscosity scales as:

$$\eta_0 \sim c^{3.4/(3\nu-1)} \cdot M_w^{3.4}$$

in the reptation regime. Practically, expect $\eta_0 \sim 0.1$-$10$ PaÂ·s at 90Â°C for 6% agarose, depending on MW and measurement conditions. This is **orders of magnitude higher** than the oil phase viscosity, giving a viscosity ratio $\lambda = \mu_d/\mu_c \gg 1$, which strongly disfavors breakup to small sizes.

### 2.9 Energy Dissipation Rate from the Mixer

For a homogenizing mixer (rotor-stator type):

$$\varepsilon = \frac{P}{\rho_c V_{gap}}$$

where $P$ is the power input and $V_{gap}$ is the high-shear gap volume. For a rotor-stator:

$$P = N_p \rho_c N^3 D^5$$

where $N_p$ is the power number ($\sim 1$-$5$ depending on geometry), $N$ is rotational speed (rps), and $D$ is the rotor diameter.

Alternatively, for a high-shear homogenizer, the energy dissipation in the gap:

$$\varepsilon_{gap} = \frac{N_p N^3 D^5}{V_{gap}} = \frac{N_p N^3 D^2}{\delta}$$

where $\delta$ is the gap width. Typical values for a lab rotor-stator at 10,000 RPM: $\varepsilon \sim 10^4$-$10^6$ m$^2$/s$^3$ in the gap region.

**Source:** Padron, G.A. (2005) *Measurement and Comparison of Power Draw in Batch Rotor-Stator Mixers*, PhD thesis. Atiemo-Obeng & Calabrese (2004) in *Handbook of Industrial Mixing*.

---

## 3. OPTIMAL PARAMETERS â€” DIRECT vs. INDIRECT SOURCES

### 3.1 Stirring Speed (RPM) for Target Droplet Size ~2 Âµm

**Classification: INDIRECT**

No direct literature value exists for this specific system. Must be derived by solving the Kolmogorov-Hinze (with viscosity correction) or the viscous-subrange breakup model iteratively:

1. Estimate $\varepsilon$ from mixer geometry and RPM: $\varepsilon = N_p N^3 D^5 / V_{gap}$
2. Estimate $\sigma$ from Span-80 concentration
3. Estimate $\mu_d$ from Mark-Houwink + mixing rule
4. Solve $d_{max}(\varepsilon, \sigma, \mu_d, \rho_c) = 2\,\mu$m for $N$

**Preliminary estimate**: For $\sigma \approx 5$ mN/m, $\rho_c \approx 850$ kg/m$^3$, $\mu_d \approx 1$ PaÂ·s, and using the Davies (1985) correlation for high-viscosity drops:

$$d_{max} \approx C \sigma^{0.6} \rho_c^{-0.6} \varepsilon^{-0.4} (1 + f(\mu_d))$$

For $d_{max} = 2\,\mu$m, back-calculation gives $\varepsilon \sim 10^6$-$10^8$ m$^2$/s$^3$, corresponding to rotor-stator speeds of **15,000-25,000 RPM** with a tight-gap homogenizer, or the use of an **ultrasonic homogenizer/microfluidizer**. A conventional magnetic stirrer or even an overhead impeller is **insufficient** â€” a rotor-stator homogenizer (e.g., IKA Ultra-Turrax) or ultrasonic processor is required.

**Key uncertainty**: The extremely high $\mu_d$ of the agarose solution at 6% may make 2 $\mu$m impossible with conventional homogenization. This is a critical feasibility question that should be addressed first.

### 3.2 Surfactant Concentration (Span-80)

**Classification: SEMI-DIRECT**

- **Minimum**: Above CMC to ensure equilibrium coverage. CMC of Span-80 in mineral oil $\approx 0.5$-$1\%$ w/v (**direct**, Opawale & Burgess 1998, *J. Colloid Interface Sci.* 197:142).
- **Practical range**: 2-10% v/v relative to oil phase (**direct** from emulsion literature for similar W/O systems).
- **Optimal**: Must be determined by experiment or PBE simulation balancing breakup facilitation (lower $\sigma$) against coalescence suppression. **Indirect** for precise optimization.
- **Note**: Excess Span-80 forms reverse micelles that can extract water from droplets, potentially affecting polymer concentration inside the droplet.

### 3.3 Polymer Concentration (Total and Ratio)

**Classification: DIRECT (partially)**

- **Agarose**: 6% w/v following Zhao et al. (2020) â€” **direct**.
- **Chitosan**: The 7:3 or 8:2 mass ratio gives chitosan at $\sim 1.5$-$2.6\%$ w/v â€” **specified by the user's design**.
- **Effect on emulsification**: Higher total polymer concentration â†’ higher $\mu_d$ â†’ harder to break into small droplets. The agarose-chitosan ratio affects phase separation morphology during gelation â†’ affects pore structure â€” **indirect** relationship requiring phase diagram mapping.

### 3.4 Oil Phase Temperature

**Classification: DIRECT**

Must be above $T_{gel}$ of agarose ($\sim 36$-$42Â°$C) to prevent premature gelation during emulsification. The user's 90Â°C is well above this and is chosen to match the agarose dissolution temperature. Practical range: 60-90Â°C.
- Lower temperature â†’ higher oil viscosity â†’ lower $\lambda$ â†’ easier breakup, BUT risk of premature gelation.
- Higher temperature â†’ lower oil viscosity â†’ higher $\lambda$ â†’ harder breakup, but guaranteed sol state.
- **Optimal**: 70-80Â°C (balance between sol state and manageable viscosity ratio). **Indirect** â€” requires coupled modeling of $\mu_d(T)$, $\mu_c(T)$, and $T_{gel}$.

### 3.5 Cooling Rate

**Classification: INDIRECT**

No single optimal value in literature for this system. Depends on desired pore structure:
- **Fast cooling** (>10Â°C/min): finer pores, more nucleation events
- **Slow cooling** (<2Â°C/min): coarser pores, more coarsening
- For target 60-100 nm: likely requires **moderate** cooling rate (2-10Â°C/min) â€” must be determined from phase-field or kinetic Monte Carlo simulation of the gelation/phase-separation process.

### 3.6 Genipin Concentration and Crosslinking Time

**Classification: INDIRECT (but guided by literature)**

- **Genipin concentration**: 0.1-1.0% w/v relative to chitosan mass is typical (**direct**, Sung et al. 1999; Muzzarelli 2009). Molar ratio of genipin to free $-NH_2$ groups is the meaningful parameter: typically 0.05-0.5 mol/mol.
- **Crosslinking time**: Depends on temperature and desired degree of crosslinking. At 37Â°C, significant crosslinking in 1-6 hours; plateau at 24 hours (**direct** from Dimida et al. 2017).
- **Optimal for target mechanical properties**: Must solve the kinetics ODE coupled with rubber elasticity to find the crosslinking degree that gives target $G$. **Indirect**.

### 3.7 Oil:Water Volume Ratio

**Classification: SEMI-DIRECT**

- Typical W/O emulsion for microsphere preparation: oil:water = 5:1 to 20:1 (v/v) (**direct** from general microsphere literature).
- Higher oil fraction â†’ more dilute droplets â†’ less coalescence â†’ narrower size distribution.
- For 2 $\mu$m microspheres: recommend **10:1 to 20:1** to minimize coalescence at high shear. **Direct** guidance from emulsion literature (Bibette et al. 1999; Mason & Bibette 1997).

---

## 4. MODELING AND SIMULATION METHODOLOGY

### 4.1 Level 1 â€” Emulsification (Droplet Formation)

**Recommended approach: Tiered strategy**

**Tier 1A â€” Analytical/0D PBE (start here):**
- Solve the 0D (spatially homogeneous) PBE using the Method of Moments (MOM) or Quadrature Method of Moments (QMOM) (McGraw 1997; Marchisio et al. 2003).
- Use the Alopaeus et al. (2002) breakage kernel (accounts for dispersed phase viscosity).
- Use the Coulaloglou & Tavlarides coalescence kernel.
- Estimate $\varepsilon$ from the power input correlation for the specific mixer.
- **Advantage**: Fast (seconds on a laptop), gives the steady-state droplet size distribution as a function of RPM, $\sigma$, $\mu_d$.
- **Disadvantage**: Assumes spatially uniform $\varepsilon$ (unrealistic for a rotor-stator where $\varepsilon$ varies by orders of magnitude between the gap and the bulk).
- **Software**: Custom Python/MATLAB script or OpenQBMM library.

**Tier 1B â€” CFD-PBE (if 0D is insufficient):**
- Resolve the spatial variation of $\varepsilon(\mathbf{x})$ using RANS (Reynolds-Averaged Navier-Stokes) with a $k$-$\varepsilon$ turbulence model.
- Couple with the PBE using QMOM or class methods (CM) at each CFD cell.
- Use OpenFOAM (multiphaseEulerFoam + PBE extension) or ANSYS Fluent (built-in PBE module).
- **Advantage**: Captures the inhomogeneous $\varepsilon$ field; critical for rotor-stator mixers where breakup is localized.
- **Disadvantage**: Computationally expensive (days on a workstation for 3D transient); requires good mesh resolution in the gap region.
- **Geometry**: Needs CAD model of the actual mixer (IKA T25 or similar). If not available, use a simplified 2D axisymmetric model.

**Tier 1C â€” Lattice Boltzmann or VOF (NOT recommended for production):**
- Could resolve individual droplet breakup events at 2 $\mu$m, but requires resolution $\sim 0.1\,\mu$m per cell over a domain of $\sim 1$-$10$ mm â†’ $10^{12}$+ cells. **Computationally prohibitive** for the full emulsification process. Only useful for validating breakup models on a single-droplet scale.

**Recommendation**: Start with **Tier 1A** for rapid parameter screening and design-space exploration. Move to **Tier 1B** only if the 0D model gives unrealistic distributions or if spatial optimization of the mixer geometry is needed.

### 4.2 Level 2 â€” Gelation and Pore Formation

**Recommended approach: Phase-field model coupled with thermal transport**

The **Cahn-Hilliard** equation for spinodal decomposition:

$$\frac{\partial \phi}{\partial t} = \nabla \cdot \left[M(\phi) \nabla \left(\frac{\delta F}{\delta \phi}\right)\right]$$

where $\phi$ is the local polymer composition (e.g., agarose volume fraction or the composition difference between agarose-rich and chitosan-rich phases), $M$ is the mobility, and $F$ is the Flory-Huggins free energy functional:

$$F[\phi] = \int \left[f_{FH}(\phi, T) + \frac{\kappa}{2}|\nabla\phi|^2\right] dV$$

The temperature $T(t)$ enters through the $T$-dependence of $\chi(T)$ and through the gelation kinetics (which arrest the phase separation when the agarose gels).

**Coupling with gelation**: When $\alpha(T,t)$ exceeds a critical value ($\alpha_{gel} \approx 0.05$-$0.1$, percolation threshold), the mobility $M$ drops to zero (gelation arrests coarsening):

$$M(\phi, \alpha) = M_0(\phi) \cdot (1 - \alpha/\alpha_{arrest})^+ $$

**Implementation**:
- **FiPy** (Python, NIST) or **MOOSE** (Idaho National Lab) for phase-field PDEs on a 2D/3D domain.
- Domain: a single microsphere (sphere of radius 1 $\mu$m), exploit spherical symmetry â†’ reduce to 1D radial problem or 2D slice.
- Grid: $\Delta x \sim 1$-$5$ nm to resolve the 60-100 nm pore features â†’ $\sim 400$-$2000$ cells in 1D, or $\sim 10^5$-$10^6$ cells in 2D.
- Time stepping: implicit scheme to handle the stiff kinetics; total time $\sim 10$-$100$ s (cooling duration).

**Alternative â€” Kinetic Monte Carlo (KMC)**:
- Model gelation as a bond-percolation process on a lattice.
- Each lattice site represents a polymer segment; bonds form with temperature-dependent probability.
- Can naturally capture the stochastic nature of gelation and the resulting pore-size distribution.
- **Advantage**: Naturally produces a distribution of pore sizes (not just the mean).
- **Disadvantage**: Harder to parametrize from continuum thermodynamic data; results depend on lattice choice.

**Recommendation**: **Phase-field** as the primary method (more direct connection to measurable thermodynamic quantities), with KMC as a validation/sensitivity check.

### 4.3 Level 3 â€” Crosslinking Optimization

**Recommended approach: ODE system with reaction-diffusion PDE if needed**

**Step 1 â€” ODE kinetics** (spatially homogeneous, valid if Thiele modulus $\Phi = R\sqrt{k[NH_2]/D_{gen}} \ll 1$):

For microsphere radius $R = 1\,\mu$m, $D_{genipin} \sim 10^{-10}$ m$^2$/s (in hydrogel), $k \sim 10^{-3}$ L/(molÂ·s), $[NH_2] \sim 10$ mol/L:

$$\Phi = R\sqrt{k[NH_2]/D_{gen}} \approx 10^{-6} \sqrt{10^{-3} \times 10 / 10^{-10}} = 10^{-6} \times 10^{3.5} \approx 3 \times 10^{-3}$$

$\Phi \ll 1$, so **diffusion is NOT rate-limiting** at 2 $\mu$m scale. The ODE approach is fully justified. No need for a reaction-diffusion PDE.

The ODE system:

$$\frac{d[X]}{dt} = k(T) [NH_2]_0(1-2\alpha_x) \cdot ([Gen]_0 - [X])$$
$$\alpha_x = [X]/([NH_2]_0/2)$$
$$G(\alpha_x) = (n_{e,0} + \Delta n_e(\alpha_x)) k_B T$$
$$\xi(\alpha_x) = 0.071 \cdot \nu_{2,s}(\alpha_x)^{-1/3} \cdot M_c(\alpha_x)^{1/2}$$

This gives $G(t)$ and $\xi(t)$ as functions of crosslinking time, genipin concentration, and temperature.

**Software**: SciPy `odeint` / MATLAB `ode45`.

### 4.4 Level 4 â€” Mechanical Property Prediction

**Recommended approach: Analytical rubber elasticity + FEM validation**

**Step 1 â€” Analytical**:
For the double network:

$$G_{DN} = G_{agarose}(\nu_{2,agar}, T_{gel}) + G_{chitosan}(\alpha_x, n_{e,genipin})$$

$$G_{agarose} \approx \frac{\rho_{agar} R T}{M_{c,agar}} \nu_{2,agar}^{1/3}$$

$$G_{chitosan} = \frac{n_{e,genipin} k_B T}{\nu_{2,chit}^{-1/3}}$$

$$E_{DN} \approx 3 G_{DN} \quad (\nu \approx 0.5)$$

**Step 2 â€” FEM** (if detailed stress distribution needed):
- Model the microsphere as a hyperelastic sphere (Neo-Hookean or Ogden) under compression.
- Use COMSOL, Abaqus, or FEniCS.
- Include biphasic (poroelastic) behavior if drained compression is relevant.

**Recommendation**: Analytical models are sufficient for optimization; FEM only for validation against AFM nanoindentation or microcompression data.

### 4.5 Integration Strategy

The four levels couple as follows:

```
Level 1: Emulsification PBE
    Input: RPM, Span-80 conc., T_oil, Ï†_d, Î¼_d(T), Ïƒ(Span-80)
    Output: Droplet size distribution n(d)
        â†“
Level 2: Phase-field gelation
    Input: d from Level 1, cooling rate dT/dt, polymer concentrations, Ï‡(T)
    Output: Pore structure (pore size distribution, connectivity)
        â†“
Level 3: Crosslinking kinetics
    Input: Network from Level 2, genipin concentration, time, T
    Output: Crosslink density n_e(t), mesh size Î¾(t)
        â†“
Level 4: Mechanical properties
    Input: n_e from Level 3, network topology from Level 2
    Output: G, E, compression resistance
```

**The coupling is sequential (one-way)**: each level's output feeds the next. No feedback loops are needed because:
- Droplet formation is complete before gelation begins (the emulsion is formed at 90Â°C, gelation occurs below 42Â°C)
- Gelation and phase separation are complete before crosslinking (genipin is added after microsphere recovery)
- Crosslinking determines the final mechanical properties

This sequential coupling is a **major simplification** that makes the problem tractable.

**Overall Optimization Framework:**

**Recommended: Bayesian Optimization (BO) with Gaussian Process surrogate**

Rationale:
- The parameter space has ~7 dimensions (RPM, Span-80 conc., polymer ratio, oil temperature, cooling rate, genipin conc., crosslinking time).
- Each simulation run (Levels 1-4) takes minutes to hours (dominated by the phase-field in Level 2).
- BO is sample-efficient (finds optima in 50-200 evaluations), ideal for expensive simulations.
- Use a multi-objective acquisition function (Expected Hypervolume Improvement) to simultaneously optimize: (a) droplet size â†’ 2 $\mu$m, (b) pore size â†’ 60-100 nm, (c) mechanical modulus â†’ target value.

**Software**: BoTorch (PyTorch-based BO library) or GPyOpt.

**Alternative**: If computational budget allows >1000 evaluations (e.g., using only Tier 1A + analytical models for Levels 2-4), use a **Design of Experiments (DOE)** approach with a **Response Surface Method (RSM)**: Central Composite Design or Box-Behnken, fit quadratic response surfaces, optimize analytically.

**Validation Strategy:**

1. **Level 1**: Compare predicted droplet size distributions with optical microscopy / laser diffraction measurements at 3-5 different RPM values.
2. **Level 2**: Compare predicted pore sizes with SEM / cryo-SEM measurements at 3 different cooling rates.
3. **Level 3**: Compare predicted crosslink density with swelling ratio measurements ($Q = 1/\nu_{2,s}$) at 3 different genipin concentrations.
4. **Level 4**: Compare predicted modulus with AFM nanoindentation or rheometry on bulk gel samples.
5. **Chromatographic validation**: Measure $K_{av}$ for a series of protein standards and compare with the predicted pore-size distribution using the Ogston-Laurent-Killander model:

$$K_{av} = \exp\left[-\pi (r_s + r_f)^2 n_f L\right]$$

where $r_s$ is the solute hydrodynamic radius, $r_f$ is the fiber radius, $n_f$ is the fiber number density, and $L$ is the column length.

---

## 5. KEY LITERATURE REFERENCES

### Emulsification and Droplet Breakup

1. **Hinze, J.O.** (1955) Fundamentals of the hydrodynamic mechanism of splitting in dispersion processes. *AIChE J.* 1:289-295. â€” Original Kolmogorov-Hinze theory.

2. **Kolmogorov, A.N.** (1949) On the disintegration of drops in a turbulent flow. *Doklady Akad. Nauk SSSR* 66:825-828. â€” Foundational turbulence-breakup scaling.

3. **Coulaloglou, C.A. & Tavlarides, L.L.** (1977) Description of interaction processes in agitated liquid-liquid dispersions. *Chem. Eng. Sci.* 32:1289-1297. â€” Breakage and coalescence kernel models.

4. **Luo, H. & Svendsen, H.F.** (1996) Theoretical model for drop and bubble breakup in turbulent dispersions. *AIChE J.* 42:1225-1233. â€” Improved breakage model.

5. **Alopaeus, V., Koskinen, J., Keskinen, K.I.** (2002) Simulation of the population balances for liquid-liquid systems in a nonideal stirred tank. Part 2. *Chem. Eng. Sci.* 57:1815-1825. â€” Viscosity-corrected breakage kernel.

6. **Calabrese, R.V., Chang, T.P.K., Dang, P.T.** (1986) Drop breakup in turbulent stirred-tank contactors. Part I: Effect of dispersed-phase viscosity. *AIChE J.* 32:657-666. â€” Modified Hinze for viscous drops.

7. **Davies, J.T.** (1985) Drop sizes of emulsions related to turbulent energy dissipation rates. *Chem. Eng. Sci.* 40:839-842. â€” Practical correlation for emulsions.

8. **Ramkrishna, D.** (2000) *Population Balances: Theory and Applications to Particulate Systems in Engineering*. Academic Press. â€” Definitive PBE textbook.

9. **Marchisio, D.L., Vigil, R.D., Fox, R.O.** (2003) Quadrature method of moments for aggregation-breakage processes. *J. Colloid Interface Sci.* 258:322-334. â€” QMOM solution method.

### Agarose Gelation and Pore Structure

10. **Arnott, S., Fulmer, A., Scott, W.E., et al.** (1974) The agarose double helix and its function in agarose gel structure. *J. Mol. Biol.* 90:269-284. â€” Agarose helix structure.

11. **Aymard, P., Martin, D.R., Plucknett, K., Foster, T.J., Clark, A.H., Norton, I.T.** (2001) Influence of thermal history on the structural and mechanical properties of agarose gels. *Biopolymers* 59:131-144. â€” Gelation kinetics and cooling rate effects.

12. **Djabourov, M., Clark, A.H., Rowlands, D.W., Ross-Murphy, S.B.** (1989) Small-angle X-ray scattering characterization of agarose sols and gels. *Macromolecules* 22:180-188. â€” Agarose network structure.

13. **Normand, V., Lootens, D.L., Amici, E., Plucknett, K.P., Aymard, P.** (2000) New insight into agarose gel mechanical properties. *Biomacromolecules* 1:730-738. â€” Mechanical properties and Flory-Huggins parameter.

### Chitosan, Genipin, and Crosslinking

14. **Muzzarelli, R.A.A.** (2009) Genipin-crosslinked chitosan hydrogels as biomedical and pharmaceutical aids. *Carbohydr. Polym.* 77:1-9. â€” Comprehensive review of genipin-chitosan chemistry.

15. **Butler, M.F., Ng, Y.-F., Pudney, P.D.A.** (2003) Mechanism and kinetics of the crosslinking reaction between biopolymers containing primary amine groups and genipin. *J. Polym. Sci. Part A: Polym. Chem.* 41:3941-3953. â€” Reaction mechanism.

16. **Dimida, S., Demitri, C., De Benedictis, V.M., Scalera, F., Gervaso, F., Sannino, A.** (2017) Genipin-cross-linked chitosan-based hydrogels: Reaction kinetics and structure-related characteristics. *J. Appl. Polym. Sci.* 134:45006. â€” Crosslinking kinetics.

17. **Mi, F.-L., Sung, H.-W., Shyu, S.-S.** (2000) Synthesis and characterization of a novel chitosan-based network prepared using naturally occurring crosslinker. *J. Polym. Sci. Part A: Polym. Chem.* 38:2804-2814. â€” Genipin crosslinking of chitosan.

### Double-Network Hydrogels

18. **Gong, J.P., Katsuyama, Y., Kurokawa, T., Osada, Y.** (2003) Double-network hydrogels with extremely high mechanical strength. *Adv. Mater.* 15:1155-1158. â€” Original DN concept.

19. **Gong, J.P.** (2010) Why are double network hydrogels so tough? *Soft Matter* 6:2583-2590. â€” DN mechanics theory.

### Polymer Physics and Swelling

20. **Flory, P.J.** (1953) *Principles of Polymer Chemistry*. Cornell University Press. â€” Flory-Rehner theory, rubber elasticity.

21. **Canal, T. & Peppas, N.A.** (1989) Correlation between mesh size and equilibrium degree of swelling of polymeric networks. *J. Biomed. Mater. Res.* 23:1183-1193. â€” Mesh size equation.

22. **Zhao, J., Li, Q., Huayan, S., et al.** (2020) Preparation and chromatographic property evaluation of macroporous double-network agarose-based microspheres. *Eng. Life Sci.* 20:504-513. â€” The reference paper for this system.

### Surfactant and Interfacial Science

23. **Opawale, F.O. & Burgess, D.J.** (1998) Influence of interfacial properties of lipophilic surfactants on water-in-oil emulsion stability. *J. Colloid Interface Sci.* 197:142-150. â€” Span-80 CMC and interfacial properties.

24. **Grace, H.P.** (1982) Dispersion phenomena in high viscosity immiscible fluid systems and application of static mixers as dispersion devices in such systems. *Chem. Eng. Commun.* 14:225-277. â€” Critical capillary number vs. viscosity ratio.

### Phase-Field Modeling

25. **Cahn, J.W. & Hilliard, J.E.** (1958) Free energy of a nonuniform system. I. Interfacial free energy. *J. Chem. Phys.* 28:258-267. â€” Cahn-Hilliard equation.

### Heat Transfer

26. **Incropera, F.P. & DeWitt, D.P.** (2002) *Fundamentals of Heat and Mass Transfer*, 5th ed. Wiley. â€” Lumped capacitance model, convective heat transfer.

### Optimization

27. **Snoek, J., Larochelle, H., Adams, R.P.** (2012) Practical Bayesian optimization of machine learning algorithms. *NeurIPS* 25:2951-2959. â€” BO methodology.

---

## CRITICAL FEASIBILITY ASSESSMENT

Before proceeding with the full simulation pipeline, there is one **critical feasibility concern** that must be addressed:

**Achieving 2 Âµm droplets with 6% agarose + chitosan solution may be extremely difficult or impossible via conventional rotor-stator homogenization.** The dispersed-phase viscosity is likely $\sim 0.1$-$10$ PaÂ·s at 90Â°C (highly entangled polymer solution), giving viscosity ratios $\lambda \gg 1$. The Grace curve shows that for $\lambda > 4$, droplet breakup in simple shear becomes impossible regardless of capillary number. In extensional flow (which rotor-stator mixers partially provide), breakup remains possible but requires very high $Ca$.

**Recommended pre-study**: Before committing to the full simulation campaign, perform a **single-point experimental test** or a **0D PBE calculation** (Tier 1A, ~1 day of work) to verify that 2 Âµm is achievable. If not, consider:
1. Reducing polymer concentration (but this affects pore structure)
2. Using an ultrasonic homogenizer or microfluidizer instead of rotor-stator
3. Using membrane emulsification (Shirasu Porous Glass) for precise size control
4. Accepting larger microspheres (~10-50 Âµm) and adjusting the target application accordingly

This feasibility check should be **Step 0** of the computational program.