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
STRATEGIES = {
    "repeat_y1_south": np.arange(1, 14),
    "repeat_y2_north": np.arange(14, 27),
    "alternate_ns": np.ravel(np.column_stack([np.arange(1, 14), np.arange(14, 27)]))[:13],
    "ecliptic_y5": np.arange(56, 69),
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
    p.add_argument("--max-sector", type=int, default=mission.LAST_COMPLETE_SECTOR)
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

    windows = mission.load_orbit_windows(max_sector=args.max_sector)
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

    n_ambiguous = int((results["flag"] == al.AMBIGUOUS).sum())
    print(f"\nambiguous systems available to resolve: {n_ambiguous:,}")

    scores = []
    for name, pointings in STRATEGIES.items():
        pointings = np.asarray(pointings)[:args.n_future]
        slots = np.arange(args.max_sector + 1,
                          args.max_sector + 1 + len(pointings))
        future = mission.synthesize_windows(windows, slots, pointings)

        scored = experiment.score_strategy(results, future, sets)
        scored.to_pickle(args.out / f"results_{name}.pkl")

        resolved = scored["resolved"]
        scores.append({
            "strategy": name,
            "n_resolved": int(resolved.sum()),
            "frac_ambiguous_resolved": float(resolved.sum() / n_ambiguous),
            "stars_resolved": int(scored.loc[resolved, "weight"].sum()),
        })
        print(f"  {name:18s} resolved {resolved.sum():6,d} "
              f"({100 * resolved.sum() / n_ambiguous:5.2f}%)  "
              f"{scored.loc[resolved, 'weight'].sum():7,d} stars")

    score_df = pd.DataFrame(scores).sort_values("stars_resolved", ascending=False)
    score_df.to_csv(args.out / "strategy_scores.csv", index=False)

    print(f"\nwrote outputs to {args.out}/")


if __name__ == "__main__":
    main()
