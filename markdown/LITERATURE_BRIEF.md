# Context brief: literature search for the TESS long-period window-function project

**Paste this whole file into a fresh conversation as the opening message.** It
is written to be read cold by an assistant with no prior context.

---

## What I need from you

A literature review supporting a paper in preparation. I have the code and the
simulation working; what I lack is grounding in prior work. I need to know what
has already been done, what conventions I should match, whose results I should
compare against, and where the genuine gap is.

**Please search rather than recall.** Use ADS, arXiv and the web. Section 6
below lists leads from memory — treat every one of them as a *guess that needs
verifying*, including author names and years. Some are probably wrong or
misattributed. Tell me which ones don't check out; that is useful information,
not a failure.

---

## 1. Who I am

Observational astronomer working on TESS detection of long-period transiting
planets. Comfortable with the domain — you do not need to explain what a
transit or an ephemeris is. I prototype in notebooks and care about methodology
being defensible in referee report terms.

## 2. The project in one paragraph

TESS observes most of the sky in ~27 day blocks separated by months or years.
A planet whose period exceeds a sector's span cannot produce the three-plus
evenly spaced transits a box-least-squares search needs. I inject planets at
arbitrary period and conjunction phase, determine which transits land inside
the *real* TESS observing windows (built from MIT's per-orbit table, not an
idealised 27 d block), and classify each outcome as **missed / mono-transit /
ambiguous duo / solved**. Two goals: (1) score proposed *future* TESS pointing
strategies by how many otherwise-unusable planets each converts into usable
ones, and (2) measure TESS completeness to long-period planets.

## 3. Technical definitions you need

**The alias ladder.** Two transits separated by a gap `dt` are consistent with
*any* period `dt/n` for integer `n`. I enumerate `dt/n` from the smallest
observed gap down to a floor of 13 d (below which an ordinary search would have
found the planet).

**Rule-out criteria.** An alias dies if it (a) predicts a detectable transit
inside a window that came up empty, or (b) fails to reproduce a transit that was
observed. Whatever survives is the candidate period set. **Invariant: the true
period must never be ruled out** — this is my primary correctness check.

**Outcomes**, ordered: `NO_TRANSIT(0) < MONO(1) < AMBIGUOUS(2) < SOLVED(3)`.

**Strategy scoring.** A proposed future schedule is expressed as extra
observing windows, each tagged with which past sector's sky footprint it
re-observes. I re-run the whole injection-recovery against baseline + new
windows and count outcome *upgrades*.

**Detection model.** A transit counts if the entire event `tc ± T/2` fits
inside a single window, with `T = 13 hr · (P/365 d)^(1/3)` (central transit,
solar density). **No SNR, depth, or systematics modelling** — this is a pure
window function, an upper bound on recovery.

## 4. Numbers already established, for comparison against the literature

- Baseline: TESS sectors 1–106, 320 orbit windows, 2018-07-25 → 2026-08-09
  (2,937 d).
- Catalog: 8,158 CTLv8 targets with coverage in that span; 7,446 seen in ≥2
  sectors, collapsing to 2,200 distinct sector combinations.
- Outcome fractions over a log-spaced 20–500 d grid with uniform conjunction
  phase: **solved 28.9%, missed 27.9%, mono 22.9%, ambiguous 20.3%.**
- A 13-sector extended mission re-observing Year-1 fields upgraded ~6.6% of all
  injected systems. `mono → ambiguous` and `ambiguous → solved` are comparable
  in size and trade places between grids; `no transit → mono` is close behind.
  The point is that **upgrades out of mono and missed states together outweigh
  the classic "resolve a duo" channel**, which is why I changed what the study
  measures.
- Sector geometry is not uniform: spans 24.3–57.4 d, 2–8 orbit segments per
  sector. Sectors 97–98 are ~55 d double-length pointings. Sectors 99–107 roll
  the field of view 40° and shift toward middle ecliptic latitudes (median
  |ecliptic latitude| of observed stars drops 53.8° → 41.2°).

## 5. The four questions

### 5.1 Prior duo-transit and mono-transit studies

Who has done this, for TESS, Kepler and K2? I specifically want:

- **Yield predictions** for single- and double-transit candidates, and how they
  defined the categories. Do others use my four-way classification, or a
  different one? What vocabulary dominates — "monotransit", "duotransit",
  "single-transit candidate"? I want to match established terminology rather
  than invent it.
- **Real recovered systems** — planets whose periods were pinned down after a
  mono or duo detection, and by what means (follow-up photometry, TESS
  re-observation, RVs). Rough count of how many exist.
- **Occurrence-rate work** that uses long-period transits with few events, since
  that is what a completeness function ultimately feeds.
- Whether anyone has published a TESS window function specifically for the
  long-period regime. **If this exists, I need to know immediately** — it is
  the closest prior art and determines how I frame the paper.

### 5.2 Strategies for resolving ambiguous periods

- What methods are used in practice: targeted photometric follow-up at
  predicted alias times, ground-based networks (NGTS, and whichever others are
  active), space-based follow-up (CHEOPS is the obvious one), RVs, archival
  data mining.
- How do people **choose which alias to observe first**? Is there a published
  prioritisation scheme, or is it ad hoc?
- Is there existing software for this? I have a vague recollection of a
  `MonoTools` package that does alias probability weighting — verify whether it
  exists, who wrote it, what it does, and whether my rule-out logic duplicates
  or complements it.
- Has anyone framed *TESS re-observation itself* as an alias-resolution
  strategy, which is exactly what my paper argues?

### 5.3 The mathematics of transit aliases

This is where I most suspect there is relevant work outside the exoplanet
literature that I am unaware of.

- Is there a formal treatment of the `P = dt/n` ladder — its statistics, the
  expected number of surviving aliases, the prior on `n`?
- **Minimum-observation problem:** given transits at known times and a set of
  candidate periods, what is the minimum number (or total duration) of
  additional observations needed to guarantee disambiguation? This smells like a
  known combinatorial problem. Possible framings worth checking: **sparse /
  Golomb rulers**, **covering and separating systems**, **group testing**,
  **optimal experimental design** (D-optimality etc.). Has anyone applied any of
  these to transit scheduling?
- **Radial-velocity aliasing** is a mature analogous literature — sampling
  aliases in periodograms, and formal methods for distinguishing true periods
  from aliases. What transfers?
- Any information-theoretic treatment: expected information gain per unit
  observing time, which is the natural way to rank my candidate strategies.

### 5.4 TESS observing-strategy decisions, Year 8/9 and beyond

- What drove the Year 8 change — the 40° field-of-view roll and the shift to
  middle ecliptic latitudes? Find the actual justification, not my inference
  from the pointing table.
- Why are sectors 97 and 98 double-length (~55 d, eight orbit segments)? Was
  that a deliberate long-period-sensitivity experiment, a scheduling
  constraint, or something else? **This matters a lot to me** — if it was
  deliberate, my paper is contributing to an active decision process rather
  than commenting from outside.
- What is publicly known or proposed for Year 9 and later? Senior review
  material, extended-mission proposals, community white papers, GI program
  documentation.
- What *constrains* the pointing? Spacecraft thermal and pointing limits,
  downlink geometry, Moon and Earth avoidance, camera roll constraints. I need
  to know which of my hypothetical strategies are physically impossible, so I
  don't propose something that cannot be flown.
- Prior work evaluating *candidate* TESS extended-mission strategies by
  simulated planet yield. If someone has already built this comparison
  machinery, I want their metrics.

## 6. Leads from memory — **all unverified, verify or discard**

Author/year attributions below are recalled, not checked. Expect errors.

*Long-period and single-transit detection:*
Cooke, Pollacco, Bayliss and collaborators (TESS/NGTS single-transit yields);
Hugh Osborn (K2 and TESS single transits, and the possible `MonoTools`
package); Gill et al. (NGTS recovery of a TESS single transit — NGTS-11b is the
example I half-remember); Ulmer-Moll et al. (NGTS/CHEOPS long-period recovery);
Villanueva, Dragomir & Gaudi (single-transit occurrence-rate method);
Battley et al.; Eisner et al. (Planet Hunters TESS, which surfaces long-period
single transits by visual inspection).

*Kepler/K2 long-period occurrence with few transits:*
Foreman-Mackey, Morton, Hogg, Agol & Schölkopf (~2016); Wang et al.;
Herman, Zhu & Wu.

*TESS yield and extended-mission planning:*
Ricker et al. 2015 (mission paper, certain); Sullivan et al. 2015;
Barclay, Pepper & Quintana 2018; Bouma et al. ~2017 on extended-mission
simulations — **this one sounds closest to my strategy-comparison work and
should be checked first**; Kunimoto et al. on extended-mission yields.

*Aliasing methodology:*
Dawson & Fabrycky (~2010) on RV aliases — I am fairly confident this exists and
is the best analogue for the alias-discrimination problem.

## 7. What I want back

1. **An annotated bibliography**, grouped by the four questions above.
   For each entry: full citation, what they did, what is directly reusable,
   and how it relates to my method specifically. BibTeX-ready if convenient.
2. **A short synthesis** — where the genuine gap is, and how I should position
   this paper relative to the closest prior work.
3. **A terminology note** — the conventional vocabulary and category
   definitions I should adopt so my results are comparable to others'.
4. **Explicit flags** on anything in §6 that you could not verify, and anything
   you found that contradicts an assumption stated in this brief.
5. If you find that §5.1's last question is a *yes* — someone has already
   published a long-period TESS window function — say so up front rather than
   burying it.

## 8. Scope

Literature and synthesis only. No code, and no need to open the repository —
everything you need about the method is in this document. If a question turns
on a methodological detail I have not specified, ask rather than assume. 
Prioritize peer-reviewed work and work published by the official TESS, HEASARC, and MIT teams.
