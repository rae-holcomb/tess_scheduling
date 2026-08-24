# Literature review: TESS long-period window function project

Compiled 2026-08-21 in response to `LITERATURE_BRIEF.md`. Everything here comes
from live searches (ADS/arXiv/journal pages/MIT/HEASARC), not recall. Confidence
is flagged per item. Fields I could not confirm from a primary page are marked
`[unverified]`.

---

## 0. Headline answer to §5.1's last question

**No — nobody has published the thing you are building. But the gap is narrower
than you probably think, and one paper is much closer than the rest.**

There is no published "TESS long-period window function" in the sense of an
all-sky, per-orbit-window, phase-marginalised map of *what period regimes are
recoverable where*. But three groups have built real TESS window functions and
used them for adjacent purposes:

| Paper | Real per-orbit windows? | Alias ladder? | Scores future pointings? |
|---|---|---|---|
| Cooke et al. 2019 (A&A 631, A83) | **Yes** (SPOC 2-min timestamps, S1–11) | No | One fixed scenario (Y1→Y4) |
| **Cooke et al. 2021 (MNRAS 500, 5088)** | **Yes** | **Yes** | No (scores *follow-up*, not pointings) |
| Rodel et al. 2024 (MNRAS 529, 715) | **Yes** (FITS TIME headers, gaps >0.5 d) | No | No |
| Bouma et al. 2017; Kunimoto et al. 2022 | Idealised | No | **Yes**, but by total planet yield |

**Cooke et al. 2021 is your closest prior art** and you must engage with it
directly in the introduction. It does the `P = dt/n` enumeration, it rules
aliases out using TESS coverage itself, and it reports a resolution fraction.
What it does *not* do is the thing your paper does — treat the pointing schedule
as the free variable. See §3 (Synthesis) for how to position against it.

Nobody has combined (a) real windows + (b) alias-state classification + (c)
future-schedule scoring. That intersection is your paper.

---

## 1. Annotated bibliography

### 1.1 Prior mono-/duo-transit studies (brief §5.1)

**Cooke, B. F., Pollacco, D., West, R., McCormac, J., & Wheatley, P. J. 2018,
A&A, 619, A175** — "Single site observations of TESS single transit detections".
*Verified: exists, A&A 2018. Volume/page `[unverified]`.*
The origin of the "~10% of TESS discoveries will be monotransits" figure your
brief's §6 half-remembers. Simulates the two-year primary mission with emphasis
on single transits, and assesses single-site ground recovery.
**Reusable:** their monotransit yield normalisation; the framing that TESS's
sector strategy *increases* monotransit probability rather than merely failing.

**Cooke, B. F., Pollacco, D., & Bayliss, D. 2019, A&A, 631, A83** — "An
examination of the effect of the TESS extended mission on southern hemisphere
monotransits". *Verified from journal page.*
Directly analogous to your strategy-scoring goal, for one fixed scenario. Uses
**realistic window functions derived from actual SPOC 2-minute light curves,
sectors 1–11**, extrapolated to Year 4. Headline numbers:
- 339 monotransits predicted from Year 1 (southern hemisphere)
- **~80% (266/339) transit again during Year 4 re-observation**
- 149 remain monotransits in the combined Y1+Y4 dataset
- of the 266 recovered, **75% (189) show only a *single* transit in Year 4** —
  i.e. they become duotransits with an alias ladder, not solved systems.

**This is the single most directly comparable result to your 6.6% upgrade
figure**, and it is the number a referee will ask you to reconcile with. Note
the structural agreement with your finding: re-observation mostly moves systems
*mono → ambiguous*, not *ambiguous → solved*. Cooke et al. observed the same
effect and treated it as a problem; you are treating it as a measurable quantity.
That is a defensible reframing, but say so explicitly.

**Cooke, B. F., Pollacco, D., Anderson, D. R., Bayliss, D., Bouchy, F., Gill, S.,
Grieves, N., Lendl, M., Nielsen, L. D., Udry, S., & Wheatley, P. J. 2021, MNRAS,
500, 5088–5097**, doi:10.1093/mnras/staa3569 — "Resolving period aliases for TESS
monotransits recovered during the extended mission". *Verified from journal page,
full author list confirmed. Cited variously as 2020 or 2021 — ADS bibcode is
`2021MNRAS.500.5088C`; use 2021.*

**Read this one first.** Direct quotes from the methodology:

> "Based on only two transits we cannot infer the true orbital period... The
> maximum period is given by the separation of the two observed transits which is
> on the order of ~2 years. The minimum period is given by the constraint that
> only one transit is seen in each TESS observing run, approximately 10 d for a
> single sector of observations."

> "For each alias, we create a list of all times at which the system would transit
> (including times before the Year 1 transit and after the Year 3 transit) and
> check if there are any TESS observations within half a transit duration of any
> of these."

> "If no transit is seen it shows that the true period is not one of those that
> predicted a transit this night and we therefore remove all periods that
> predicted a transit that night from our list of aliases."

That second quote is **your rule-out criterion (a), already published**, and the
third is criterion (b). Your logic does not duplicate `MonoTools` — it duplicates
*this*. Cite it as the origin of the method and claim the extension, not the idea.

Key numbers for direct comparison:
- ~400 systems with one transit in each of primary and extended missions
- **average 38 period aliases per duotransit system**
- 207/395 (~52%) fully resolved by 50-day NGTS photometry + CORALIE spectroscopy

Note the **13 d vs 10 d floor discrepancy**: they justify ~10 d as "only one
transit per sector"; you use 13 d as "an ordinary BLS search would have found it".
These are different arguments for the same quantity. Yours is defensible but you
must state the justification, because the two conventions give different alias
counts and therefore non-comparable "ambiguous" fractions.

**Hawthorn, F., et al. 2024, MNRAS, 528, 1841–1862**, doi:10.1093/mnras/stad3783
— "TESS duotransit candidates from the Southern Ecliptic Hemisphere". *Verified.*
85 duotransit candidates (25 known, 60 new) from Cycles 1 and 3; 8 < Tmag < 14,
depths 0.1–1.8%, durations 2–10 h, periods typically >20 d with the lower bound
set by TESS observation duration.
**Reusable:** this is the *observed* duo population your simulated AMBIGUOUS class
should be validated against. If your ambiguous fraction predicts far more or
fewer than 85 for that footprint, you have a calibration problem worth explaining.

**Rodel, T., Bayliss, D., Gill, S., & Hawthorn, F. 2024, MNRAS, 529, 715–731**,
doi:10.1093/mnras/stae474 — "TIaRA TESS 1: estimating exoplanet yields from Years
1 and 3 SPOC light curves". *Verified; author list and pages confirmed.*
The nearest thing to a published TESS long-period completeness function. Their
window function construction, in their words:

> "To use the real window function for each star in the SPOC FFI list, we
> determine blocks of continuous TESS observations for each target using the
> timestamps from the FITS file TIME header for each TESS sector."

They identify gaps >0.5 d (perigee downlink gaps and anomalies); one CVZ example
is 79.3% duty cycle. Transit durations from Winn (2014), trapezoidal SNR.
Detection uses an incomplete-gamma detection-efficiency function from Kepler,
**linearly extrapolated down to 1–2 transits** — this is exactly the soft spot
your pure-geometry upper bound sidesteps.
Headline: **403(+64/−38) detections at P > 25 d**, 215(+37/−23) monotransits,
113(+22/−13) biennial duotransits, 110 at P > 100 d. Giant-planet completeness
falls from >80% at P < 6.25 d to **<10% for 200–400 d**.
**How it relates to you:** TIaRA is the SNR-aware version of your calculation for
a subset of the sky. Your numbers should *bound* theirs from above. If they do
not, something is wrong. This is your best quantitative validation check and I
would build a figure around it.

**Villanueva, S., Jr., Dragomir, D., & Gaudi, B. S. 2019, AJ, 157, 84** — "An
Estimate of the Yield of Single-transit Planetary Events from TESS". *Verified as
AJ, Jan 2019, doi:10.3847/1538-3881/aaf85e. Volume/page `[unverified]`.*
First TESS single-transit yield estimate; predicts ~100× more single transits than
Kepler. **Methodologically the closest to yours in spirit**: uses an *analytic
integration* for the probability of exactly one vs two-or-more transits over a
finite baseline, rather than Monte Carlo injection. Cited elsewhere as predicting
~1100 single-transit events in the primary mission.
**Reusable:** their analytic P(1 transit), P(≥2 transits) expressions are the
closed-form limit of your simulation for an idealised single window. Reproducing
them analytically for one sector would be a cheap, convincing correctness check to
put in an appendix — and it is a *different* check from your "true period never
ruled out" invariant, which only tests the alias logic, not the window logic.

**Osborn, H. P., et al. 2016, MNRAS, 457, 2273** — "Single transit candidates from
K2: detection and period estimation". *Verified.*
Introduces **NAMASTE** ("An MCMC Analysis of Single Transit Exoplanets"),
estimating period from transit duration + stellar density. 7 candidates from K2
C1–3; EPIC 203311200 at 540(+410/−230) d.
**Reusable:** the duration→period inversion is the *complementary* constraint to
your combinatorial rule-out. Worth one paragraph — a referee will ask why you do
not use it, and the honest answer is that it is a prior on `n`, not a rule-out,
and it is degenerate with eccentricity.

**Eisner, N. L., et al. 2020/2021, Planet Hunters TESS I & II** — *Verified that
the series exists; PHT I is TOI-813, an 84-d planet on a subgiant (arXiv:1909.09094).
Volumes `[unverified]`.* PHT has surfaced **90 candidates, 73 of them single
transits**. Follow-ups combine citizen-science labels with CNNs for non-phase-folded
single-event vetting (see also *Accelerating Long-period Exoplanet Discovery by
Combining Deep Learning and Citizen Science*, AJ 2025, doi:10.3847/1538-3881/add46d).
**Relevance:** evidence that the BLS-based detection assumption underlying your
13 d floor is *not* how many real long-period detections happen. Worth a caveat
sentence.

**Kepler/K2 long-period occurrence from few transits:**
- **Foreman-Mackey, D., Morton, T. D., Hogg, D. W., Agol, E., & Schölkopf, B.
  2016, AJ, 152, 206** — "The Population of Long-period Transiting Exoplanets".
  *Verified, exactly as recalled in §6.* First fully automated search for
  long-period planets with only one or two transits; 39,036 bright GK dwarfs.
  **This is the template for how to turn a few-transit detection into an
  occurrence rate, and therefore the paper your completeness function is
  ultimately feeding.** Their Section on the search-completeness simulation is
  the structural analogue of what you are doing for TESS.
- **Herman, M. K., Zhu, W., & Wu, Y. 2019, AJ, 157, 248** — "Revisiting the
  Long-period Transiting Planets from Kepler". *Verified, exactly as recalled.*
  Re-runs the Foreman-Mackey pipeline on a larger sample with Gaia DR2 stellar
  parameters; several previous candidates become false positives. 15 candidates,
  0.3–1 R_J, most with 2–10 yr periods.
- **Wang, J., et al. 2015, ApJ, 815, 127** — "Planet Hunters VIII: Characterization
  of 41 Long-Period Exoplanet Candidates from Kepler Archival Data". *Verified.*
  41 candidates in 38 systems; **17 with one transit, 14 with two transits**.
  Note the explicit mono/duo split — useful precedent for reporting your outcome
  classes as counts rather than only fractions.

### 1.2 Resolving ambiguous periods (brief §5.2)

**`MonoTools` — VERIFIED, IT EXISTS.** Hugh Osborn, `github.com/hposborn/MonoTools`,
ASCL entry `2022ascl.soft04020O`, cited as **Osborn 2022**. Description: "A package
for detecting, vetting, and modelling transiting exoplanets on uncertain periods."

What it actually does, and how it relates to your rule-out logic:
> it "calculates a marginalised probability distribution across all allowed aliases
> for a given transit model by combining priors for each alias"

with two stated assumptions: (i) short periods are favoured over long ones "due to
a combination of geometric probability **and window function**", and (ii) planets
in multiplanet systems have low eccentricities. It fits transits agnostic of
period, then computes a PDF from the orbital velocity implied by impact parameter,
radius ratio, duration and stellar density.

**Your logic complements it; it does not duplicate it.** `MonoTools` puts a
*probability* on each surviving alias (soft, physical, per-system). You determine
*which aliases survive at all* (hard, combinatorial, per-schedule). Note the
important detail for your framing: **`MonoTools` already uses a window function as
part of its alias prior** — so a referee could reasonably ask whether your
population-level window function is the correct object to feed that prior. That is
an argument *for* your paper, and I would make it explicitly: you are computing,
at population level and as a function of schedule, the quantity `MonoTools`
currently approximates per-system.

**Follow-up in practice — the three channels, all verified:**
- **Ground, NGTS.** The dedicated NGTS duotransit programme monitors candidates
  with depths ≳1000 ppm and schedules observations "when one of the period aliases
  is predicted to transit". NGTS-11 b (TOI-1847 b) is the canonical success:
  **Gill, S., et al. 2020, ApJL**, doi:10.3847/2041-8213/ab9eb9 — *the first
  exoplanet discovered from a TESS single-transit event*, recovered after
  **79 nights** of NGTS monitoring; the second transit restricted the period to
  **13 discrete aliases**, resolved to 35.46 d with RVs. Also **Gill, S., et al.
  2020, MNRAS, 491, 1548** — "NGTS and WASP photometric recovery of a single-transit
  candidate from TESS" (*archival* recovery, a distinct and cheaper channel).
- **Space, CHEOPS.** **Osborn, H. P., et al. 2020, MNRAS, 494, 736** — "CHEOPS
  observations of TESS primary mission monotransits". Plus a dedicated CHEOPS GTO
  programme for TESS duotransits, targeting depths ~2–4 ppt and small planets that
  cannot be done from the ground; worked examples in **Osborn et al. 2023, MNRAS,
  523, 3069** (HIP 9618) and the HD 15906 paper, **MNRAS, 523, 3090** (*first
  author `[unverified]`*). Also **Garai, Z., et al. 2023, A&A** on HD 22946 d.
- **RVs / long-term spectroscopy.** Ulmer-Moll et al. 2022, 2023 (A&A); the
  CORALIE arm of Cooke et al. 2021; TOI-5678 b (48 d, CHEOPS+HARPS).

**"Which alias first?" — there is no published, general prioritisation scheme.**
This is a real finding, and I searched for it specifically. Practice is:
1. `MonoTools` alias probabilities (the closest thing to a principled ranking);
2. RV pre-screening — e.g. NGTS-38 b, where the 180.5 d alias was targeted as the
   best match to existing RVs;
3. observability/schedulability — whichever alias transits next from a given site.
No paper I found frames this as an optimisation problem with an objective function.
**This is an exploitable gap adjacent to your paper** — see §3.

**Has anyone framed TESS re-observation itself as the alias-resolution strategy?**
**Partially — Cooke et al. 2021, and you must not overclaim here.** They use TESS
coverage as one of three rule-out mechanisms (the first quote in §1.1). What is
*not* done anywhere I could find: treating the *future pointing schedule* as the
decision variable and scoring candidate schedules by alias-state upgrades. Phrase
your novelty claim at that level of precision and it will survive review.

**Also relevant, and recent:** *Orbital Periods and Equilibrium Temperatures from
Single TESS Transits with a Physics-Informed Neural Network* (arXiv:2608.01101,
Aug 2026; *authors `[unverified]`*). Recovers period from duration via Kepler's
third law with a PINN marginalising over unobserved transit geometry; median
absolute error 40.5% vs BLS 79.5%. Explicitly notes "with one transit, every trial
period longer than the observing baseline fits the data identically." Two weeks
old — worth a citation to show currency, and it makes the point that the
central-transit assumption (which your `T = 13 hr·(P/365)^(1/3)` also makes)
biases periods low.

### 1.3 The mathematics of transit aliases (brief §5.3)

**Short version: your instinct is right that there is a formal literature, and
your instinct is also right that nobody has connected it to transit scheduling.**

**Formal treatment of the `P = dt/n` ladder:** the closest is
**Becker, J. C., Vanderburg, A., Rodriguez, J. E., et al. 2019, AJ** (arXiv:1809.10688)
— "A Discrete Set of Possible Transit Ephemerides for Two Long Period Gas Giants
Orbiting HIP 41378". *Verified it exists and is about exactly the discrete-alias
problem; volume/page `[unverified]`, and I could not confirm from the abstract page
whether they use the integer-divisor enumeration explicitly.* They narrow the
discrete set using transit durations, stellar properties, system dynamics, and
**archival HATNet/KELT/WASP coverage** — i.e. your rule-out criterion (a) applied
to ground surveys — then recommend targeted follow-up at the most probable times.
**This is the single-system version of your calculation.** Worth reading in full.

**RV aliasing — mature and transferable:**
- **Dawson, R. I., & Fabrycky, D. C. 2010, ApJ, 722, 937** — "Radial Velocity
  Planets De-aliased: A New, Short Period for Super-Earth 55 Cnc e". *Verified;
  your §6 recollection is correct.* Compares amplitude *and phase* of predicted
  aliases against peaks in the data using the spectral window function; revised
  55 Cnc e from 2.8 d to 0.74 d. Presented as a "cookbook".
  **What transfers:** the conceptual move of treating the *window function as the
  object that generates the aliases*, which is precisely your paper's thesis. What
  does *not* transfer: their aliases are additive in frequency (`f_alias = f_true ± k·f_sample`),
  yours are integer divisors of a gap. Different algebra, same epistemics. Say this
  explicitly — it is the kind of distinction referees like.
- **Baluev, R. V. 2012, MNRAS, 422, 2372–2385** — "Distinguishing between a true
  period and its alias, and other tasks of model discrimination". *Verified.*
  Applies the **Vuong closeness test** (Kullback–Leibler-based) to decide whether
  two rival period models are yet observationally distinguishable, or still
  equivalent. Asymptotically normal, works under misspecification.
  **Directly reusable:** this is a principled, published answer to "are these two
  aliases distinguishable given the data I have?" — the statistical counterpart of
  your hard rule-out. If you want a probabilistic version of your AMBIGUOUS class,
  this is the tool.

**Optimal experimental design / adaptive scheduling — the real find:**
- **Dzigan, Y., & Zucker, S. 2011/2012, MNRAS** (arXiv:1105.5393) — "Directed
  follow-up strategy of low-cadence photometric surveys in search of transiting
  exoplanets — I. Bayesian approach for adaptive scheduling". *Verified it exists.*
  **This is the closest published work to your §5.3 minimum-observation question**
  and it was not in your §6 list. Bayesian adaptive scheduling of follow-up for
  transit candidates from sparse photometry.
- **Loredo, T. J., Berger, J. O., Chernoff, D. F., Clyde, M. A., & Liu, B. 2011**
  (arXiv:1108.0020) — "Bayesian Methods for Analysis and Adaptive Scheduling of
  Exoplanet Observations". Adaptive Bayesian exploration maximising information per
  observation; finds most information is gained within the next orbit or two, with
  preferred later times of decreasing expected gain.
  **This gives you the information-theoretic ranking you asked for in §5.3's last
  bullet, essentially off the shelf**: expected information gain = reduction in
  entropy over the alias set. For your discrete alias set, the entropy is just
  `log(#surviving aliases)` under a flat prior, and a candidate observing window's
  expected gain is analytic. That is a very cheap upgrade to your scoring metric
  and I would strongly consider adding it — it converts "count of upgrades" into a
  quantity with units of bits, which is far more defensible in referee terms.

**Sparse rulers / covering systems / group testing — negative result.** I searched
specifically for Golomb rulers, covering and separating systems, and group testing
applied to transit or observation scheduling, and **found nothing**. The Golomb
ruler literature is real and well developed (phased arrays, X-ray crystallography,
error-correcting codes) but I found no application to transit period disambiguation.
Structurally your problem *is* a group-testing / separating-system problem: each
candidate observing window is a test that partitions the alias set into
"predicted to transit" and "not", and you want a minimal set of tests that
separates all aliases. **This connection appears to be genuinely unmade in the
literature.** That is either a nice contribution or a sign that it is harder than
it looks (the tests are not freely choosable — they are constrained by the
spacecraft), and I would not build the paper on it. But one paragraph framing your
rule-out procedure as a separating system, with a reference to the group-testing
literature, would be distinctive and cheap.

### 1.4 TESS observing strategy, Year 8/9 and beyond (brief §5.4)

**Sectors 97/98, the four-orbit pointings — likely deliberate and likely
long-period-motivated, but I could not fully close this.** *Read the flag carefully.*

Confirmed from MIT: Year 8 = Sectors 97–107, 2025-09-15 → 2026-09-07. Sectors 97
(orbits 201–204, Sep 15 – Nov 9 2025) and 98 (orbits 205–208, Nov 9 2025 – Jan 5
2026) are **four consecutive orbits** instead of the standard two. Sectors 99–107
have the FOV "rotated about its axis by 40° and shifted to provide overlapping
observations of the middle latitudes", all southern.

**The MIT Year 8 page states no scientific rationale** — I fetched it and checked.
The mechanism and motivation appear in the *Guide for TESS Extended Mission
Planning* (HEASARC): the anti-solar pointing requirement (nominally within 15°,
no more than 30° in ecliptic longitude) is what limits a sector to two orbits, and
**relaxing it toward ±30° permits four-orbit, ~54 d sectors**. The stated science
argument is the long-period one: with periods approaching a sector length, few
transits are recorded, and when the star is next observed the gap is large, so
"the orbital period will not be well defined, making follow-up difficult."

> **FLAG — verify before you rely on this.** That last connection comes from a
> search engine's summary of the HEASARC planning-guide PDF. My PDF fetch returned
> undecodable binary and PDF tooling was unavailable in this environment, so **I
> did not read the primary source directly.** Given your brief says this "matters a
> lot", download
> `https://heasarc.gsfc.nasa.gov/docs/tess/docs/Guide-to-TESS-for-EM-planning.pdf`
> and confirm the ±30° / four-orbit / ~54 d passage yourself. If it says what the
> summary claims, **your paper is contributing to an active decision process and
> you should say so in the first paragraph of the introduction.** Also check TESS
> Data Release Notes for Sector 97 (DRN132) and the 2025 Astrophysics Senior Review
> TESS submission, neither of which I could read directly.

**Pointing constraints — what is and is not flyable.** From the TESS Observatory
Guide and mission documentation (verified via search summary of primary docs):
- **Anti-solar constraint.** Combined FOV centre ideally within **15°** of the
  antisolar direction in ecliptic longitude, **no more than 30°**. This is a
  *power* constraint (solar panel illumination), and it is the binding one.
- **Thermal/sunshade.** The sunshade geometry makes it hard to hold a field for
  **more than two spacecraft orbits (~42 d)** — this is precisely what the
  four-orbit sectors are relaxing.
- **Earth and Moon.** Passages through the FOV cause scattered light; affected
  regions are ~120° in extent in ecliptic latitude. Avoidable by *timing*, not by
  pointing alone.
- **Nominal cadence.** ~28° eastward in ecliptic longitude per lunar month.
- **Sector geometry.** 13 partially overlapping sectors per hemisphere, 24°×96°,
  from ecliptic latitude 6° to the pole.
**Implication for your strategy scoring:** any proposed window must sit at an
ecliptic longitude within 30° of antisolar *at that date*. This is a hard filter
you can apply mechanically, and applying it will make your proposed strategies
credible. **I would add this as an explicit feasibility check in the paper** — it
is the difference between "here are some windows" and "here is a flyable schedule",
and it is the first thing a TESS-team referee will look for.

**Prior work scoring candidate TESS strategies by simulated yield:**
- **Bouma, L. G., et al. 2017** (arXiv:1705.08891) — "Planet Detection Simulations
  for Several Possible TESS Extended Missions". *Verified. Your §6 instinct that
  this is closest to your strategy-comparison work is correct — for the machinery,
  though not for the metric.* Monte Carlo over **six** one-year extended-mission
  scenarios. Named strategies now standard vocabulary: **`pole`** (stare at an
  ecliptic pole for ≥1 yr), **`hemi`** (repeat the primary mission on one
  hemisphere), **`all-sky`** (both hemispheres in one year at ~14 d/sector).
  *`[unverified]`: whether this was ever refereed or remains arXiv-only.*
  **Use their scenario names for your hypothetical strategies.** Free comparability.
- **Kunimoto, M., Winn, J., Ricker, G., & Vanderspek, R. 2022, AJ, 163, 290** —
  "Predicting the Exoplanet Yield of the TESS Prime and Extended Missions Through
  Years 1–7". *Verified.* 9.4M AFGKM stars from CTL v8.01; 4719±334 (Prime),
  3707±209 (EM1), 4093±180 (EM2), total 12519±678. **>1200 with P > 20 d in EM2.**
  Explicitly models "photometric performance, temporal window functions, and transit
  detection probability".
  **This is your closest methodological sibling on the yield side, and it uses the
  same catalogue lineage as your 8,158 CTLv8 targets.** Their P > 20 d subsample is
  the natural denominator for your outcome fractions.
- **Sullivan, P. W., et al. 2015, ApJ, 809, 77** — original TESS yield simulation.
  Superseded but always cited.
- **Barclay, T., Pepper, J., & Quintana, E. V. 2018, ApJS, 239, 2** — "A Revised
  Exoplanet Yield from TESS". *Verified as existing with that title.* The standard
  correction to Sullivan et al.
- **Ricker, G. R., et al. 2015, JATIS, 1, 014003** — mission paper. Certain.

**Year 9 / Cycle 9.** Confirmed: **Sectors 108–121, 2026-09-07 → 2027-09-19**, third
extended mission. TESS-point has been updated with these sectors — **so the Cycle 9
pointing table already exists and you can score it directly rather than
hypothetically.** That is a strong, timely addition to the paper. Data products:
200 s FFIs, ~8000 targets at 120 s, ~2000 at 20 s per sector. Primary reference:
`TESS_Cycle9_D3CS.pdf` (ROSES-2025 D.3, Phase-1 proposals due 2026-03-10).

---

## 2. Terminology note (brief §7.3)

Adopt these. They are what the TESS/NGTS/CHEOPS community actually uses.

| Use | Not | Why |
|---|---|---|
| **monotransit** | "single-transit candidate", "mono" | Dominant in TESS-era literature (Cooke, Osborn, Ulmer-Moll, NGTS papers). "Single transit" is used by US-based groups (Villanueva et al.; Foreman-Mackey et al.) and in Kepler-era work — acceptable but weaker. Define once as equivalent. |
| **duotransit** | "duo", "double-transit" | Established by Cooke et al. 2019/2021 and Hawthorn et al. 2024. |
| **biennial duotransit** | — | Rodel et al. 2024's term for the specific Year 1 + Year 3 case. Use it when you mean that; it is more precise than "duotransit". |
| **period alias** / **alias** | "candidate period" | Universal. |
| **solved** / "period recovered" | — | Community phrasing is "resolving period aliases" (Cooke et al. 2021 title). Your `SOLVED` label is fine; introduce it as "period-resolved". |

**On your four-way scheme.** Nobody uses your exact ordered four-way classification
— the literature uses mono/duo as *categories of candidate*, not as *outcome states
of a survey*. That difference is a genuine contribution, not a problem, but it means
you must:
1. define the four states explicitly and early;
2. give the mapping: your `MONO` = the literature's monotransit, your `AMBIGUOUS` =
   duotransit (or multi-transit with surviving aliases), your `SOLVED` = a
   period-resolved system, your `NO_TRANSIT` has **no counterpart in the
   literature** because existing studies condition on detection;
3. **be careful with `AMBIGUOUS`.** In your scheme it can include systems with ≥3
   transits that still have surviving aliases. That is broader than "duotransit".
   Either rename it or state the superset relation, otherwise your 20.3% is not
   comparable to Hawthorn et al.'s 85 duos.

Two conventions to declare because they break comparability if silent:
- **Your 13 d alias floor** vs Cooke et al.'s ~10 d (§1.1).
- **Your central-transit duration** `T = 13 hr·(P/365)^(1/3)`, vs the Winn (2014)
  formulation with impact parameter and eccentricity used by Rodel et al. Yours is
  the `b=0`, `e=0` limit, so it is an upper bound on duration and therefore on
  recoverability — consistent with your "upper bound" framing, but state it.

---

## 3. Synthesis: where the gap is, and how to position

**The gap, stated precisely.** Three literatures exist and do not intersect:

1. **Real-window TESS completeness** (Rodel et al. 2024; Cooke et al. 2019) —
   builds true window functions, but asks "how many planets?", not "which period
   states?"
2. **Alias resolution** (Cooke et al. 2021; Osborn's `MonoTools`; Becker et al. 2019)
   — enumerates and kills aliases, but takes the observing schedule as *given* and
   optimises *follow-up*.
3. **Pointing-strategy comparison** (Bouma et al. 2017; Kunimoto et al. 2022) —
   treats the schedule as the free variable, but scores it by *total planet yield*,
   which is dominated by short periods and is therefore nearly blind to exactly the
   effect you care about.

**Your paper is the intersection: scoring candidate TESS schedules by their effect
on the alias state of long-period systems.** That is a real and defensible gap.

**How to position, concretely:**

- **Lead with Cooke et al. 2021, generously.** They published the rule-out logic.
  Your contribution is not the method, it is (i) applying it across the *whole*
  real 8-year baseline rather than a Y1/Y3 pair, (ii) the four-state outcome
  algebra, and (iii) making the schedule the decision variable. Claiming the
  alias-ladder method as novel would be the fastest route to a hostile referee
  report.
- **Make the metric change your headline.** Your finding that upgrades out of MONO
  and NO_TRANSIT together outweigh the classic duo-resolution channel is the paper's
  most interesting claim, *and it is corroborated by Cooke et al. 2019's 189/266
  result* (most re-observed monotransits become duotransits, not solved systems).
  Present it as: the community optimises for resolving duos, but the schedule's
  leverage is mostly upstream of that. That is a clean, quotable thesis.
- **Validate against TIaRA and Hawthorn.** Your pure-geometry numbers should bound
  Rodel et al. 2024 from above and should be consistent with Hawthorn et al. 2024's
  85 observed duos. Two external checks cost you one figure each and buy a lot of
  credibility, because your primary correctness check (the true period is never
  ruled out) is internal and self-referential — it validates the alias logic but
  says nothing about whether the window function is right.
- **Score Cycle 9, not just hypotheticals.** Sectors 108–121 are public. Scoring the
  *actual* proposed schedule alongside your hypotheticals converts the paper from a
  methods note into an input to a live decision.
- **Apply the anti-solar feasibility filter** (§1.4) to every proposed strategy, and
  say you did.
- **Consider adding expected information gain in bits** (Loredo et al. 2011;
  Dzigan & Zucker 2011) as a second metric alongside upgrade counts. Cheap, and it
  answers "how do you rank strategies?" in a way that "6.6% of systems upgraded"
  does not.

**Two things I would push back on in the brief:**
- Your detection model has no SNR *and* no eccentricity. Rodel et al. 2024 show
  completeness for long-period *giants* is <10% at 200–400 d — so the gap between
  your geometric upper bound and reality is large and period-dependent. Frame the
  paper as a window function (a necessary condition) and be explicit that
  multiplying by a detection efficiency is future work, otherwise the outcome
  fractions read as yields and will be misread.
- Your BLS-based justification for the 13 d floor sits awkwardly with the fact that
  a large fraction of real long-period TESS detections came from **visual inspection
  and ML on non-phase-folded data** (Planet Hunters TESS: 73 of 90 candidates were
  single transits). The floor is still a reasonable modelling choice; just justify
  it as a definitional convention rather than as a statement about what searches
  can find.

---

## 4. Verification table for brief §6

| Lead as recalled | Status | Correction |
|---|---|---|
| Cooke, Pollacco, Bayliss — TESS/NGTS single-transit yields | **Verified** | Three papers: 2018 A&A (with West, McCormac, Wheatley); 2019 A&A 631 A83; 2021 MNRAS 500 5088. |
| Hugh Osborn — K2/TESS single transits, `MonoTools` | **Verified** | Osborn et al. 2016 MNRAS 457 2273 (NAMASTE); `MonoTools` = Osborn 2022, ASCL 2204.020. Both real. |
| Gill et al. — NGTS recovery of a TESS single transit, "NGTS-11b" | **Verified, and you half-remembered *two* papers** | Gill et al. 2020 ApJL = NGTS-11 b/TOI-1847 b. Separately Gill et al. 2020 MNRAS 491 1548 = NGTS+WASP archival recovery. |
| Ulmer-Moll et al. — NGTS/CHEOPS long-period recovery | **Verified** | 2022 and 2023 A&A. Also lead author of the 2026 review (below). |
| Villanueva, Dragomir & Gaudi — single-transit occurrence method | **Verified** | AJ 2019, doi:10.3847/1538-3881/aaf85e. Method is analytic, not Monte Carlo. |
| Battley et al. | **Weakly verified** | The author is real and works on TESS single transits, but I could not pin a specific paper matching your description. Search ADS directly. |
| Eisner et al. — Planet Hunters TESS | **Verified** | PHT I (TOI-813, 84 d) and PHT II. 90 candidates, 73 single transits. |
| Foreman-Mackey, Morton, Hogg, Agol & Schölkopf ~2016 | **Verified exactly** | AJ 152, 206 (2016). |
| Wang et al. | **Verified** | Planet Hunters VIII, ApJ 815, 127 (2015). 41 candidates, 17 mono / 14 duo. |
| Herman, Zhu & Wu | **Verified exactly** | AJ 157, 248 (2019). |
| Ricker et al. 2015 | **Verified** | JATIS 1, 014003. |
| Sullivan et al. 2015 | **Verified** | ApJ 809, 77. |
| Barclay, Pepper & Quintana 2018 | **Verified** | "A Revised Exoplanet Yield from TESS", ApJS 239, 2. |
| Bouma et al. ~2017 — "closest to my work, check first" | **Verified; instinct correct** | arXiv:1705.08891. Closest on *machinery*; but scores yield, not alias states. Gave us `pole`/`hemi`/`all-sky`. |
| Kunimoto et al. — extended-mission yields | **Verified** | Kunimoto, Winn, Ricker & Vanderspek 2022, AJ 163, 290. |
| Dawson & Fabrycky ~2010 — RV aliases | **Verified exactly** | ApJ 722, 937 (2010). Best analogue, with the caveat in §1.3. |

**Not in your §6, and you should have them:**
- **Rodel et al. 2024 (TIaRA TESS 1)** — real TESS window functions, long-period
  completeness. Your closest completeness comparison.
- **Hawthorn et al. 2024** — 85 observed duotransit candidates. Your validation set.
- **Baluev 2012, MNRAS 422, 2372** — formal test for alias distinguishability.
- **Dzigan & Zucker 2011** and **Loredo et al. 2011** — Bayesian adaptive scheduling
  and expected information gain. The answer to your §5.3 information-theory bullet.
- **Becker et al. 2019** (HIP 41378) — the single-system version of your discrete
  alias enumeration with archival rule-outs.
- **Ulmer-Moll, S., Akinsanmi, B., & Müller, S. 2026, arXiv:2604.09254** — "Long-period
  transiting exoplanets: advances in detection and characterization", chapter in the
  NCCR PlanetS Legacy Book (Springer). **Four months old and it is a review of
  exactly your subfield.** Read it cover to cover before writing your introduction;
  it will give you the current citation graph for free. Useful stat from it: 37 new
  transiting warm giants with P > 20 d published from TESS; the PlanetS Monotransit
  Initiative alone accounts for 13 published + 13 solved-and-forthcoming planets.

---

## 5. Things I could not resolve — worth your own follow-up

1. **The Sector 97/98 rationale, from a primary source.** See the FLAG in §1.4.
   This is the highest-value open item given how much it matters to your framing.
2. **The 2025 Astrophysics Senior Review TESS submission.** Not findable via search;
   likely the authoritative statement of the Year 8/9 strategy rationale. Try the
   NASA Astrophysics Senior Review page directly, or ask the TESS Science Office.
3. **Whether Bouma et al. 2017 was ever refereed** or remains an arXiv preprint.
   Affects how you cite it.
4. **Exact volume/page for Villanueva et al. 2019, Gill et al. 2020 (ApJL),
   Becker et al. 2019, and the Eisner PHT papers.** I have DOIs or arXiv IDs for all
   but did not confirm the bibliographic fields; pull them from ADS before submission.
5. **First author of the HD 15906 paper** (MNRAS 523, 3090).
6. **`MonoTools`' actual window-function implementation.** Since it already folds a
   window function into its alias prior, it is worth reading
   `github.com/hposborn/MonoTools` to see how coarse that treatment is. If it is
   idealised 27 d sectors, that is a concrete, citable motivation for your work and
   arguably the strongest single argument for the paper's utility.
7. **Authors of arXiv:2608.01101** (the PINN single-transit period paper).

---

## 6. BibTeX starter

Entries below are those whose bibliographic fields I could confirm. Fields I could
not confirm are omitted rather than guessed — fill from ADS.

```bibtex
@article{Cooke2019,
  author  = {Cooke, B. F. and Pollacco, D. and Bayliss, D.},
  title   = {An examination of the effect of the {TESS} extended mission on
             southern hemisphere monotransits},
  journal = {Astronomy \& Astrophysics}, year = 2019, volume = 631, pages = {A83}
}

@article{Cooke2021,
  author  = {Cooke, B. F. and Pollacco, D. and Anderson, D. R. and Bayliss, D. and
             Bouchy, F. and Gill, S. and Grieves, N. and Lendl, M. and
             Nielsen, L. D. and Udry, S. and Wheatley, P. J.},
  title   = {Resolving period aliases for {TESS} monotransits recovered during the
             extended mission},
  journal = {MNRAS}, year = 2021, volume = 500, number = 4, pages = {5088--5097},
  doi     = {10.1093/mnras/staa3569}
}

@article{Rodel2024,
  author  = {Rodel, T. and Bayliss, D. and Gill, S. and Hawthorn, F.},
  title   = {{TIaRA TESS} 1: estimating exoplanet yields from Years 1 and 3
             {SPOC} light curves},
  journal = {MNRAS}, year = 2024, volume = 529, number = 1, pages = {715--731},
  doi     = {10.1093/mnras/stae474}
}

@article{Hawthorn2024,
  author  = {Hawthorn, F. and others},
  title   = {{TESS} duotransit candidates from the Southern Ecliptic Hemisphere},
  journal = {MNRAS}, year = 2024, volume = 528, number = 2, pages = {1841--1862},
  doi     = {10.1093/mnras/stad3783}
}

@article{ForemanMackey2016,
  author  = {Foreman-Mackey, Daniel and Morton, Timothy D. and Hogg, David W. and
             Agol, Eric and Sch{\"o}lkopf, Bernhard},
  title   = {The Population of Long-period Transiting Exoplanets},
  journal = {AJ}, year = 2016, volume = 152, pages = {206}
}

@article{Herman2019,
  author  = {Herman, Miranda K. and Zhu, Wei and Wu, Yanqin},
  title   = {Revisiting the Long-period Transiting Planets from {Kepler}},
  journal = {AJ}, year = 2019, volume = 157, pages = {248}
}

@article{Wang2015,
  author  = {Wang, Ji and others},
  title   = {Planet Hunters. {VIII}. Characterization of 41 Long-period Exoplanet
             Candidates from {Kepler} Archival Data},
  journal = {ApJ}, year = 2015, volume = 815, pages = {127}
}

@article{Kunimoto2022,
  author  = {Kunimoto, Michelle and Winn, Joshua and Ricker, George R. and
             Vanderspek, Roland K.},
  title   = {Predicting the Exoplanet Yield of the {TESS} Prime and Extended
             Missions through Years 1--7},
  journal = {AJ}, year = 2022, volume = 163, pages = {290},
  doi     = {10.3847/1538-3881/ac68e3}
}

@article{Osborn2016,
  author  = {Osborn, H. P. and others},
  title   = {Single transit candidates from {K2}: detection and period estimation},
  journal = {MNRAS}, year = 2016, volume = 457, number = 3, pages = {2273}
}

@article{Osborn2020,
  author  = {Osborn, H. P. and others},
  title   = {{CHEOPS} observations of {TESS} primary mission monotransits},
  journal = {MNRAS}, year = 2020, volume = 494, number = 1, pages = {736}
}

@misc{MonoTools,
  author = {Osborn, Hugh P.},
  title  = {{MonoTools}: Planets of uncertain periods detector and modeler},
  year   = 2022, note = {Astrophysics Source Code Library, ascl:2204.020},
  howpublished = {\url{https://github.com/hposborn/MonoTools}}
}

@article{Dawson2010,
  author  = {Dawson, Rebekah I. and Fabrycky, Daniel C.},
  title   = {Radial Velocity Planets De-aliased: A New, Short Period for
             Super-Earth 55 {Cnc} e},
  journal = {ApJ}, year = 2010, volume = 722, pages = {937}
}

@article{Baluev2012,
  author  = {Baluev, Roman V.},
  title   = {Distinguishing between a true period and its alias, and other tasks
             of model discrimination},
  journal = {MNRAS}, year = 2012, volume = 422, number = 3, pages = {2372--2385}
}

@article{Ricker2015,
  author  = {Ricker, George R. and others},
  title   = {Transiting Exoplanet Survey Satellite},
  journal = {Journal of Astronomical Telescopes, Instruments, and Systems},
  year    = 2015, volume = 1, pages = {014003}
}

@article{Sullivan2015,
  author  = {Sullivan, Peter W. and others},
  title   = {The Transiting Exoplanet Survey Satellite: Simulations of Planet
             Detections and Astrophysical False Positives},
  journal = {ApJ}, year = 2015, volume = 809, pages = {77}
}

@article{Barclay2018,
  author  = {Barclay, Thomas and Pepper, Joshua and Quintana, Elisa V.},
  title   = {A Revised Exoplanet Yield from the Transiting Exoplanet Survey
             Satellite ({TESS})},
  journal = {ApJS}, year = 2018, volume = 239, pages = {2}
}

@misc{Bouma2017,
  author = {Bouma, L. G. and others},
  title  = {Planet Detection Simulations for Several Possible {TESS} Extended
            Missions},
  year   = 2017, eprint = {1705.08891}, archivePrefix = {arXiv}
}

@misc{UlmerMoll2026,
  author = {Ulmer-Moll, Sol{\`e}ne and Akinsanmi, Babatunde and M{\"u}ller, Simon},
  title  = {Long-period transiting exoplanets: advances in detection and
            characterization},
  year   = 2026, eprint = {2604.09254}, archivePrefix = {arXiv},
  note   = {Chapter in NCCR PlanetS Legacy Book, Springer}
}

@misc{DziganZucker2011,
  author = {Dzigan, Yifat and Zucker, Shay},
  title  = {Directed follow-up strategy of low-cadence photometric surveys in
            search of transiting exoplanets -- {I}. {B}ayesian approach for
            adaptive scheduling},
  year   = 2011, eprint = {1105.5393}, archivePrefix = {arXiv}
}

@misc{Loredo2011,
  author = {Loredo, Thomas J. and Berger, James O. and Chernoff, David F. and
            Clyde, Merlise A. and Liu, Bin},
  title  = {Bayesian Methods for Analysis and Adaptive Scheduling of Exoplanet
            Observations},
  year   = 2011, eprint = {1108.0020}, archivePrefix = {arXiv}
}
```
