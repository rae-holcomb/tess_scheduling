"""Invariants the window-function calculation must satisfy.

Run with: pytest tests/
"""

import numpy as np
import pandas as pd
import pytest

from tess_window import aliases as al
from tess_window import mission, targets
from tess_window.transits import (check_observability, transit_times_from_phase,
                                  transit_times_from_tc)

DUR = 4 / 24.0


def make_windows(spans, pointings=None):
    """Build a minimal window table from (t_start, t_end) pairs."""
    spans = np.asarray(spans, dtype=float)
    n = len(spans)
    pointings = np.arange(1, n + 1) if pointings is None else pointings
    return pd.DataFrame({
        "sector": np.arange(1, n + 1), "pointing": pointings,
        "t_start": spans[:, 0], "t_end": spans[:, 1],
        "predicted": False,
    })


# --- ephemerides ---------------------------------------------------------

def test_phase_grid_spans_requested_window():
    """The transit grid must cover the window it was asked for, not a global one."""
    start, stop = 1000.0, 3000.0
    for period in [13.0, 20.8, 137.0, 400.0]:
        tt = transit_times_from_phase(period, [0.0], start, stop)[0]
        assert tt[0] == pytest.approx(start)
        assert tt[-1] >= stop - period


def test_tc_grid_stays_inside_window():
    tt = transit_times_from_tc(50.0, 1234.0, 1000.0, 2000.0)
    assert tt.min() >= 1000.0 and tt.max() < 2000.0
    assert np.allclose(np.diff(tt), 50.0)


def test_tc_grid_anchors_on_the_given_epoch():
    """tc need not be the first transit, or even inside the window."""
    tt = transit_times_from_tc(10.0, 5.0, 100.0, 200.0)
    assert np.allclose((tt - 5.0) % 10.0, 0.0)


# --- observability -------------------------------------------------------

def test_transit_must_fit_entirely_inside_a_window():
    w = make_windows([(0.0, 10.0)])
    # centred: caught. straddling the edge: not.
    caught, _, _ = check_observability(np.array([[5.0]]), [1], w, duration=DUR)
    edge, _, _ = check_observability(np.array([[9.99]]), [1], w, duration=DUR)
    assert caught[0, 0]
    assert not edge[0, 0]


def test_transit_in_a_downlink_gap_is_missed():
    w = make_windows([(0.0, 5.0), (6.0, 11.0)])
    obs, _, _ = check_observability(np.array([[5.5]]), [1, 2], w)
    assert not obs[0, 0]


def test_only_windows_covering_the_target_are_used():
    w = make_windows([(0.0, 10.0), (20.0, 30.0)])
    _, _, visible = check_observability(np.array([[5.0]]), [1], w)
    assert list(visible["sector"]) == [1]


# --- aliases -------------------------------------------------------------

def test_alias_ladder_divides_the_smallest_gap():
    _, periods = al.find_all_aliases(np.array([0.0, 100.0, 250.0]), min_period=13.0)
    assert periods[0] == pytest.approx(100.0)
    assert np.allclose(periods, 100.0 / np.arange(1, len(periods) + 1))
    assert periods.min() >= 13.0


def test_fewer_than_two_transits_yields_no_aliases():
    for times in [np.array([]), np.array([500.0])]:
        n, p = al.find_all_aliases(times)
        assert len(n) == 0 and len(p) == 0


def test_alias_that_reproduces_all_transits_survives():
    """The true period must never be ruled out by its own data."""
    empty = make_windows([(1000.0, 1010.0)])
    observed = np.array([0.0, 200.0, 400.0])
    ruled = al.rule_out_aliases(observed, np.array([200.0]), empty)
    assert not ruled[0]


def test_alias_predicting_a_transit_in_an_empty_window_is_ruled_out():
    empty = make_windows([(95.0, 105.0)])
    observed = np.array([0.0, 200.0])
    ruled = al.rule_out_aliases(observed, np.array([100.0]), empty)
    assert ruled[0]


def test_alias_that_misses_an_observed_transit_is_ruled_out():
    empty = make_windows([(1000.0, 1010.0)])
    observed = np.array([0.0, 150.0, 300.0])
    # 300 d cannot produce the transit at 150 d
    ruled = al.rule_out_aliases(observed, np.array([300.0]), empty)
    assert ruled[0]


def test_flags_track_the_number_of_observed_transits():
    empty = make_windows([(1000.0, 1010.0)])
    assert al.analyze(np.array([]), empty).flag == al.NO_TRANSIT
    assert al.analyze(np.array([5.0]), empty).flag == al.MONO


# --- end to end ----------------------------------------------------------

@pytest.fixture(scope="module")
def real_setup():
    w = mission.load_orbit_windows()
    sets = targets.unique_sector_sets(
        targets.load_targets(max_sector=98), min_sectors=2)
    return w, sets


def test_true_period_always_survives(real_setup):
    """The injected period must appear among the surviving aliases, always."""
    w, sets = real_setup
    start, stop = mission.mission_span(w)
    rng = np.random.default_rng(7)
    tested = 0

    for _ in range(600):
        period = rng.uniform(15.0, 400.0)
        target = sets.iloc[rng.integers(len(sets))]
        tt = transit_times_from_phase(period, rng.random(), start, stop)
        obs, wobs, visible = check_observability(tt, target["sectors"], w, duration=DUR)
        result = al.analyze(tt[0][obs[0]], visible[~wobs[0]], duration=DUR)

        if result.flag in (al.NO_TRANSIT, al.MONO):
            continue
        tested += 1
        assert result.flag != al.ERROR
        assert np.min(np.abs(result.surviving - period)) < 0.1

    assert tested > 100, "not enough multi-transit systems to be meaningful"


def test_synthesized_windows_inherit_recent_cadence(real_setup):
    w, _ = real_setup
    future = mission.synthesize_windows(w, np.arange(99, 105))
    bounds = mission.sector_bounds(future)

    assert future["predicted"].all()
    assert future["t_start"].min() > w["t_end"].max()
    assert np.median(np.diff(bounds["t_start"].values)) == pytest.approx(26.8, abs=1.0)
    assert (bounds["t_end"] - bounds["t_start"]).median() == pytest.approx(26.5, abs=2.0)


def test_single_sector_completion_is_normalised():
    periods = np.arange(14.0, 300.0, 1.0)
    solved, mono, missed = targets.completion_for_single_sectors(periods)
    assert np.allclose(solved + mono + missed, 1.0)
    assert (missed[periods > 100] > 0.5).all()
