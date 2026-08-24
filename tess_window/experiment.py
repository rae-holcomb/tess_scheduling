"""Batch injection-and-recovery over grids of period, phase and sky position."""

import numpy as np
import pandas as pd

from . import aliases as al
from .mission import mission_span
from .transits import (check_observability, resolve_duration,
                       transit_times_from_phase)

RESULT_COLUMNS = ["period", "tc_phase", "tc", "set_id", "num_s", "weight",
                  "ec_lat", "ec_long", "duration", "n_observed", "flag",
                  "n_aliases", "n_surviving", "surviving"]


def run_grid(periods, tc_phases, sector_sets, windows, duration=0.0,
             min_period=13.0, start_time=None, stop_time=None, progress=False):
    """Inject planets across a grid of period, phase and sky position.

    Args:
        periods: injected periods in days.
        tc_phases: conjunction phases in [0, 1).
        sector_sets: DataFrame from ``targets.unique_sector_sets``. Its index
            is recorded as ``set_id`` so results can be traced back to a sky
            position.
        windows: window table defining the observing strategy.
        duration: transit duration in days, or a callable of period such as
            ``transits.transit_duration``.

    Returns a tidy DataFrame with one row per (period, phase, sector set).
    """
    span_start, span_stop = mission_span(windows)
    if start_time is None:
        start_time = span_start
    if stop_time is None:
        stop_time = span_stop

    periods = np.atleast_1d(periods)
    tc_phases = np.atleast_1d(tc_phases)
    rows = []

    for i, period in enumerate(periods):
        if progress:
            print(f"period {i + 1}/{len(periods)}: {period:.2f} d", flush=True)

        dur = resolve_duration(duration, period)
        all_times = transit_times_from_phase(period, tc_phases, start_time, stop_time)

        for set_id, target in sector_sets.iterrows():
            transit_obs, window_obs, visible = check_observability(
                all_times, target["sectors"], windows, duration=dur)

            for j, phase in enumerate(tc_phases):
                observed = all_times[j][transit_obs[j]]
                empty = visible[~window_obs[j]]
                result = al.analyze(observed, empty, min_period=min_period,
                                    duration=dur)

                rows.append((
                    period, phase, result.tc, set_id, target["num_s"],
                    target["weight"], target["ec_lat"], target["ec_long"],
                    dur, len(observed), result.flag, len(result.alias_periods),
                    result.n_surviving, result.surviving,
                ))

    return pd.DataFrame(rows, columns=RESULT_COLUMNS)


GRID_KEYS = ["period", "tc_phase", "set_id"]

# Columns carrying the post-strategy outcome, and the baseline column each one
# falls back to for systems the strategy never revisits.
SCORED_COLUMNS = ["flag", "n_observed", "n_aliases", "n_surviving", "surviving"]


def _grid_from_results(results):
    """Recover the injection grid, and each period's duration, from a run."""
    periods = np.sort(results["period"].unique())
    phases = np.sort(results["tc_phase"].unique())
    durations = results.groupby("period")["duration"].first().to_dict()
    return periods, phases, durations


def touched_sets(sector_sets, new_windows):
    """Index labels of the sky positions a strategy actually re-observes.

    A new window only informs a target whose sky position was already covered
    by the sector whose pointing that window reuses. Everything else is
    untouched, and its outcome cannot change -- skipping those is what keeps
    scoring affordable.
    """
    pointings = set(int(p) for p in np.unique(new_windows["pointing"]))
    covered = sector_sets["sectors"].apply(
        lambda s: not pointings.isdisjoint(int(x) for x in s))
    return sector_sets.index[covered.to_numpy()]


def score_strategy(results, windows, new_windows, sector_sets, min_period=13.0,
                   progress=False):
    """Re-run the full analysis with a candidate strategy's windows appended.

    Every system is re-injected against ``windows + new_windows`` on the same
    (period, phase, sky position) grid and with the same period-dependent
    durations, so *any* outcome can improve: a missed planet can become a
    mono-transit, a mono can gain a second transit and become a duo, and a
    duo's surviving aliases can be pruned to one.

    This deliberately re-analyses rather than incrementally pruning the stored
    alias ladder. Pruning can only ever remove candidate periods, so it cannot
    promote a mono or a miss, and it compares aliases at window granularity --
    two periods that both place a transit somewhere in the same window look
    identical to it, even though their predicted transit *times* differ by days
    and the real data would separate them.

    Args:
        results: output of ``run_grid`` for the baseline mission.
        windows: the baseline window table the results were computed from.
        new_windows: extra windows from ``mission.synthesize_windows``.
        sector_sets: the same table passed to ``run_grid``.

    Returns a copy of ``results`` with the post-strategy outcome in ``*_new``
    columns, plus ``new_transits``, ``upgraded``, ``resolved`` and a readable
    ``transition`` label.
    """
    periods, phases, durations = _grid_from_results(results)
    combined = (pd.concat([windows, new_windows], ignore_index=True)
                .sort_values("t_start").reset_index(drop=True))

    keep = touched_sets(sector_sets, new_windows)
    start_time, _ = mission_span(windows)
    _, stop_time = mission_span(combined)

    if len(keep):
        extended = run_grid(periods, phases, sector_sets.loc[keep], combined,
                            duration=lambda p: durations[p],
                            min_period=min_period, start_time=start_time,
                            stop_time=stop_time, progress=progress)
    else:
        extended = pd.DataFrame(columns=RESULT_COLUMNS)

    scored = results.merge(extended[GRID_KEYS + SCORED_COLUMNS], on=GRID_KEYS,
                           how="left", suffixes=("", "_new"), validate="1:1")

    # Untouched systems keep their baseline outcome.
    for col in SCORED_COLUMNS:
        new = f"{col}_new"
        missing = scored[new].isna().to_numpy()
        if col == "surviving":       # object column of arrays; fillna can't
            scored[new] = [b if m else n for b, n, m in
                           zip(scored[col], scored[new], missing)]
        else:
            scored[new] = scored[new].where(~missing, scored[col]).astype(int)

    scored["new_transits"] = scored["n_observed_new"] - scored["n_observed"]
    scored["upgraded"] = scored["flag_new"] > scored["flag"]
    scored["resolved"] = (scored["flag_new"] == al.SOLVED) & \
                         (scored["flag"] != al.SOLVED)
    scored["transition"] = (scored["flag"].map(al.FLAG_NAMES) + " -> "
                            + scored["flag_new"].map(al.FLAG_NAMES))
    return scored


def transition_matrix(scored, weighted=True):
    """Baseline outcome (rows) against post-strategy outcome (columns).

    With ``weighted`` the cells count catalog stars rather than injected
    systems, matching how ``summarize`` reports completeness.
    """
    order = [al.NO_TRANSIT, al.MONO, al.AMBIGUOUS, al.SOLVED, al.ERROR]
    names = [al.FLAG_NAMES[f] for f in order]

    w = scored["weight"] if weighted else pd.Series(1, index=scored.index)
    table = (pd.crosstab(scored["flag"], scored["flag_new"],
                         values=w, aggfunc="sum")
             .reindex(index=order, columns=order, fill_value=0)
             .fillna(0).astype(int))
    table.index = pd.Index(names, name="before")
    table.columns = pd.Index(names, name="after")
    return table.loc[(table.sum(axis=1) > 0), :]


def summarize(results, by="period", weighted=True, flag_col="flag"):
    """Fraction of systems in each outcome class, grouped by ``by``.

    Pass ``flag_col="flag_new"`` to summarize a scored strategy instead of the
    baseline, which makes the two directly comparable.
    """
    df = results.copy()
    df["_w"] = df["weight"] if weighted else 1

    names = [(al.SOLVED, "solved"), (al.AMBIGUOUS, "ambiguous"),
             (al.MONO, "mono"), (al.NO_TRANSIT, "missed"), (al.ERROR, "error")]
    for flag, name in names:
        df[name] = (df[flag_col] == flag) * df["_w"]

    grouped = df.groupby(by, as_index=False).agg(
        solved=("solved", "sum"), ambiguous=("ambiguous", "sum"),
        mono=("mono", "sum"), missed=("missed", "sum"),
        error=("error", "sum"), total=("_w", "sum"))

    for _, name in names:
        grouped[name] = grouped[name] / grouped["total"]

    return grouped
