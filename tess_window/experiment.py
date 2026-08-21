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


def score_strategy(results, new_windows, sector_sets):
    """Re-evaluate ambiguous systems after adding new observing windows.

    Each system is re-tested with the same transit duration ``run_grid`` used
    for it, so period-dependent durations stay consistent.

    Returns a copy of ``results`` with ``flag``, ``n_surviving`` and
    ``surviving`` updated, plus a ``resolved`` column marking systems that the
    new windows newly solved.
    """
    updated = results.copy()
    updated["resolved"] = False

    for idx in results.index[results["flag"] == al.AMBIGUOUS]:
        row = results.loc[idx]
        target_sectors = sector_sets.loc[row["set_id"], "sectors"]
        visible_new = new_windows[new_windows["pointing"].isin(target_sectors)]
        if len(visible_new) == 0:
            continue

        surviving = np.asarray(row["surviving"], dtype=float)
        stub = al.AliasResult(
            flag=al.AMBIGUOUS,
            observed_times=np.array([row["tc"]]),
            alias_periods=surviving,
            ruled_out=np.zeros(len(surviving), dtype=bool),
        )
        new = al.update_with_windows(stub, row["period"], visible_new,
                                     duration=row["duration"])

        updated.at[idx, "flag"] = new.flag
        updated.at[idx, "n_surviving"] = new.n_surviving
        updated.at[idx, "surviving"] = new.surviving
        updated.at[idx, "resolved"] = new.flag == al.SOLVED

    return updated


def summarize(results, by="period", weighted=True):
    """Fraction of systems in each outcome class, grouped by ``by``."""
    df = results.copy()
    df["_w"] = df["weight"] if weighted else 1

    names = [(al.SOLVED, "solved"), (al.AMBIGUOUS, "ambiguous"),
             (al.MONO, "mono"), (al.NO_TRANSIT, "missed"), (al.ERROR, "error")]
    for flag, name in names:
        df[name] = (df["flag"] == flag) * df["_w"]

    grouped = df.groupby(by, as_index=False).agg(
        solved=("solved", "sum"), ambiguous=("ambiguous", "sum"),
        mono=("mono", "sum"), missed=("missed", "sum"),
        error=("error", "sum"), total=("_w", "sum"))

    for _, name in names:
        grouped[name] = grouped[name] / grouped["total"]

    return grouped
