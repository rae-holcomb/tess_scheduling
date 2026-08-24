"""Main batch run: TESS window function for long-period planets.

Injects planets on a log-spaced period grid across a star-weighted sample of
sky positions, classifies each as missed / mono / ambiguous / solved, then
scores candidate extended-mission strategies by how many ambiguous systems
each one resolves.

    python run_batch.py --n-sets 400 --n-phases 12 --out results/
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from tess_window import aliases as al
from tess_window import experiment, mission, targets
from tess_window.transits import transit_duration

# Candidate extended-mission strategies. Each maps new sector slots onto the
# pointing of an earlier sector, i.e. which patch of sky gets re-observed.
#
# Sector -> sky is not obvious from the sector number, so the median ecliptic
# latitude each block actually covers is noted alongside it.
STRATEGIES = {
    # southern hemisphere, median ec_lat -53 deg
    "repeat_y1_south": np.arange(1, 14),
    # northern hemisphere, median ec_lat +60 deg
    "repeat_y2_north": np.arange(14, 27),
    # interleaved south/north, one slot each
    "alternate_ns": np.ravel(np.column_stack([np.arange(1, 14), np.arange(14, 27)]))[:13],
    # Year 5: five northern sectors then eight southern ones. Despite the old
    # name "ecliptic_y5", sectors 56-68 sit at |ec_lat| ~ 55 and never approach
    # the plane -- this is a lopsided north/south mix, not an ecliptic survey.
    "y5_mixed_ns": np.arange(56, 69),
    # the sectors that really do sit on the ecliptic (median |ec_lat| = 6 deg).
    # Only ten exist, so the 13-slot budget cycles back over the first three.
    "ecliptic": np.resize(np.array([42, 43, 44, 45, 46, 70, 71, 72, 91, 92]), 13),
}


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--p-min", type=float, default=20.0)
    p.add_argument("--p-max", type=float, default=500.0)
    p.add_argument("--n-periods", type=int, default=60)
    p.add_argument("--n-phases", type=int, default=12)
    p.add_argument("--n-sets", type=int, default=400,
                   help="number of sky positions, sampled by star count")
    p.add_argument("--max-sector", type=int, default=None,
                   help="last sector to use; defaults to the last one that has "
                        "finished observing in the orbit table")
    p.add_argument("--n-future", type=int, default=13,
                   help="sectors of extended mission to simulate per strategy")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", type=Path, default=Path("results"))
    return p


def sample_sector_sets(sets, n, seed):
    """Draw sky positions with probability proportional to star count."""
    if n >= len(sets):
        return sets
    rng = np.random.default_rng(seed)
    probs = sets["weight"] / sets["weight"].sum()
    idx = rng.choice(sets.index, size=n, replace=False, p=probs)
    return sets.loc[np.sort(idx)]


def main():
    args = build_parser().parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    windows = mission.load_orbit_windows(
        max_sector="auto" if args.max_sector is None else args.max_sector)
    # Resolve "auto" to a number so the target filter and the future-sector
    # slot numbering below agree with the windows we actually loaded.
    args.max_sector = int(windows["sector"].max())
    mission.validate_sector_alignment(windows)

    all_sets = targets.unique_sector_sets(
        targets.load_targets(max_sector=args.max_sector), min_sectors=2)
    sets = sample_sector_sets(all_sets, args.n_sets, args.seed)

    periods = np.geomspace(args.p_min, args.p_max, args.n_periods)
    phases = np.linspace(0, 1, args.n_phases, endpoint=False)

    n_rows = len(periods) * len(phases) * len(sets)
    print(f"windows:  sectors 1-{args.max_sector}, {len(windows)} rows")
    print(f"sky:      {len(sets)} positions "
          f"({sets['weight'].sum()} of {all_sets['weight'].sum()} stars)")
    print(f"grid:     {len(periods)} periods x {len(phases)} phases "
          f"x {len(sets)} positions = {n_rows:,} systems\n")

    results = experiment.run_grid(periods, phases, sets, windows,
                                  duration=transit_duration, progress=True)
    results.to_pickle(args.out / "results_base.pkl")

    summary = experiment.summarize(results)
    summary.to_csv(args.out / "completion_base.csv", index=False)

    print("\nbaseline completion (star-weighted):")
    print(summary[["period", "solved", "ambiguous", "mono", "missed"]]
          .round(3).to_string(index=False))

    counts = results["flag"].value_counts()
    n_ambiguous = int(counts.get(al.AMBIGUOUS, 0))
    n_mono = int(counts.get(al.MONO, 0))
    n_missed = int(counts.get(al.NO_TRANSIT, 0))
    n_improvable = n_ambiguous + n_mono + n_missed
    print(f"\nimprovable systems: {n_improvable:,} "
          f"({n_missed:,} missed, {n_mono:,} mono, {n_ambiguous:,} ambiguous)")

    scores = []
    for name, pointings in STRATEGIES.items():
        pointings = np.asarray(pointings)[:args.n_future]
        slots = np.arange(args.max_sector + 1,
                          args.max_sector + 1 + len(pointings))
        future = mission.synthesize_windows(windows, slots, pointings)

        print(f"\nscoring {name} ...")
        scored = experiment.score_strategy(results, windows, future, sets,
                                           progress=True)
        scored.to_pickle(args.out / f"results_{name}.pkl")

        base, new = scored["flag"], scored["flag_new"]
        up = scored["upgraded"]
        row = {
            "strategy": name,
            "n_upgraded": int(up.sum()),
            "stars_upgraded": int(scored.loc[up, "weight"].sum()),
            # the original headline: duos whose period got pinned down
            "ambiguous_solved": int(((base == al.AMBIGUOUS) &
                                     (new == al.SOLVED)).sum()),
            # newly usable planets that previously had 0 or 1 transit
            "mono_to_multi": int(((base == al.MONO) &
                                  (new >= al.AMBIGUOUS)).sum()),
            "missed_to_detected": int(((base == al.NO_TRANSIT) &
                                       (new >= al.MONO)).sum()),
            "n_new_transits": int(scored["new_transits"].sum()),
        }
        row["frac_improvable_upgraded"] = row["n_upgraded"] / n_improvable
        scores.append(row)

        print(f"  upgraded {row['n_upgraded']:,} systems "
              f"({100 * row['n_upgraded'] / n_improvable:.2f}% of improvable), "
              f"{row['stars_upgraded']:,} stars")
        print(experiment.transition_matrix(scored, weighted=False)
              .to_string())

    score_df = pd.DataFrame(scores).sort_values("stars_upgraded", ascending=False)
    score_df.to_csv(args.out / "strategy_scores.csv", index=False)
    print("\nstrategy comparison:")
    print(score_df.to_string(index=False))

    print(f"\nwrote outputs to {args.out}/")


if __name__ == "__main__":
    main()
