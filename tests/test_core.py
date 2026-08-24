"""Invariants the window-function calculation must satisfy.

Run with: pytest tests/
"""

import warnings

import numpy as np
import pandas as pd
import pytest
from astropy.time import Time

from tess_window import aliases as al
from tess_window import experiment, mission, targets
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


def make_sets(*sector_lists):
    """A minimal ``unique_sector_sets`` table covering the given sectors."""
    return pd.DataFrame([
        {"sectors": np.array(s), "num_s": len(s), "weight": 1,
         "ec_lat": 0.0, "ec_long": 0.0} for s in sector_lists
    ])


def score_one(windows, new_windows, sectors, period, phase):
    """Run one injected planet through baseline scoring and return the row."""
    sets = make_sets(sectors)
    base = experiment.run_grid([period], [phase], sets, windows)
    scored = experiment.score_strategy(base, windows, new_windows, sets)
    return scored.iloc[0]


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
    last = int(w["sector"].max())
    sets = targets.unique_sector_sets(
        targets.load_targets(max_sector=last), min_sectors=2)
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
    first = int(w["sector"].max()) + 1
    future = mission.synthesize_windows(w, np.arange(first, first + 6))
    bounds = mission.sector_bounds(future)

    assert future["predicted"].all()
    assert future["t_start"].min() > w["t_end"].max()
    assert np.median(np.diff(bounds["t_start"].values)) == pytest.approx(26.8, abs=1.0)
    assert (bounds["t_end"] - bounds["t_start"]).median() == pytest.approx(26.5, abs=2.0)


# --- orbit table / pointing agreement ------------------------------------

def test_loader_excludes_sectors_still_observing():
    """A sector in progress lists only some of its segments; using it would
    look like real coverage with a fake gap at the end."""
    everything = mission.load_orbit_windows(max_sector=None)
    default = mission.load_orbit_windows()

    last = int(default["sector"].max())
    assert last == mission.last_complete_sector(everything)
    assert default["t_end"].max() < Time.now().jd
    assert int(everything["sector"].max()) >= last


def test_orbit_table_sectors_agree_with_tess_stars2px():
    """Sky coverage and observing windows are joined on sector number alone,
    so a renumbering in either source silently pairs the wrong ones."""
    w = mission.load_orbit_windows()
    report = mission.validate_sector_alignment(w)

    assert len(report) == w["sector"].nunique()
    # Renumbering shows up as tens of days; real schedule drift is ~a week.
    assert report["offset"].abs().max() < 15.0
    assert report["offset"].abs().median() < 2.0


def test_errata_restore_the_published_sector_102_103_boundary():
    """MIT's orbit table holds sector 102 until 05/04, but its own pointing
    table and HEASARC both end it on 04/21. The errata fix that."""
    raw = mission.load_orbit_windows(max_sector=None, errata=False)
    fixed = mission.load_orbit_windows(max_sector=None, errata=True)

    raw_b, fix_b = mission.sector_bounds(raw), mission.sector_bounds(fixed)
    # published file: 38.2 d sector 102 and a stunted 12.9 d sector 103
    assert raw_b.loc[102, "t_end"] - raw_b.loc[102, "t_start"] > 35
    assert raw_b.loc[103, "t_end"] - raw_b.loc[103, "t_start"] < 15
    # corrected: both are normal-length and share the 2026-04-21 boundary
    assert 24 < fix_b.loc[102, "t_end"] - fix_b.loc[102, "t_start"] < 28
    assert 24 < fix_b.loc[103, "t_end"] - fix_b.loc[103, "t_start"] < 28
    boundary = Time("2026-04-21").jd
    assert abs(fix_b.loc[102, "t_end"] - boundary) < 1.0
    assert abs(fix_b.loc[103, "t_start"] - boundary) < 1.0

    # the correction must not invent or destroy observing time
    assert len(raw) == len(fixed)
    raw_cov = (raw["t_end"] - raw["t_start"]).sum()
    assert raw_cov == pytest.approx((fixed["t_end"] - fixed["t_start"]).sum())


def test_default_load_has_no_sector_misalignment():
    """With errata applied and in-progress sectors dropped, nothing should
    disagree with tess_stars2px."""
    w = mission.load_orbit_windows()
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        report = mission.validate_sector_alignment(w)
    assert report["offset"].abs().max() < 5.0


def test_alignment_check_catches_a_renumbering():
    w = mission.load_orbit_windows()
    shifted = w.copy()
    bumped = shifted["sector"] >= shifted["sector"].max() - 2
    shifted.loc[bumped, ["t_start", "t_end"]] += 40.0

    with pytest.warns(UserWarning, match="disagree with tess_stars2px"):
        mission.validate_sector_alignment(shifted)


# --- scoring a candidate strategy ----------------------------------------

def test_a_mono_gains_a_second_transit_and_becomes_a_duo():
    """The upgrade an extended mission most wants to deliver."""
    windows = make_windows([(0.0, 10.0)], pointings=[1])
    # a window straddling the next conjunction of the P=50 d planet
    new = make_windows([(48.0, 52.0)], pointings=[2])

    row = score_one(windows, new, [1, 2], period=50.0, phase=0.0)

    assert row["flag"] == al.MONO and row["n_observed"] == 1
    assert row["flag_new"] == al.AMBIGUOUS and row["n_observed_new"] == 2
    assert row["new_transits"] == 1
    assert row["upgraded"] and not row["resolved"]
    assert row["transition"] == "mono -> ambiguous"


def test_a_mono_can_go_all_the_way_to_solved():
    """Catching the second transit *and* covering the sub-aliases pins P."""
    windows = make_windows([(0.0, 10.0)], pointings=[1])
    # (12, 46) stays empty under P=50 but kills P/2 and P/3; (48, 52) catches
    # the second transit.
    new = make_windows([(12.0, 46.0), (48.0, 52.0)], pointings=[3, 2])

    row = score_one(windows, new, [1, 2, 3], period=50.0, phase=0.0)

    assert row["flag"] == al.MONO
    assert row["flag_new"] == al.SOLVED
    assert row["n_surviving_new"] == 1
    assert row["resolved"]
    assert abs(row["surviving_new"][0] - 50.0) < 0.1


def test_a_completely_missed_planet_can_become_a_mono():
    windows = make_windows([(20.0, 30.0)], pointings=[1])
    new = make_windows([(90.0, 100.0)], pointings=[2])

    row = score_one(windows, new, [1, 2], period=50.0, phase=0.5)

    assert row["flag"] == al.NO_TRANSIT and row["n_observed"] == 0
    assert row["flag_new"] == al.MONO and row["n_observed_new"] == 1
    assert row["transition"] == "no transit -> mono"


def test_new_windows_over_unobserved_sky_change_nothing():
    """A strategy only helps targets whose field it re-observes."""
    windows = make_windows([(0.0, 10.0)], pointings=[1])
    new = make_windows([(48.0, 52.0)], pointings=[9])   # target is not in 9

    row = score_one(windows, new, [1, 2], period=50.0, phase=0.0)

    assert row["flag_new"] == row["flag"] == al.MONO
    assert row["new_transits"] == 0
    assert not row["upgraded"]


def test_touched_sets_selects_exactly_the_revisited_positions():
    sets = make_sets([1, 2], [3, 4], [4, 5])
    new = make_windows([(0.0, 1.0), (2.0, 3.0)], pointings=[4, 7])
    assert list(experiment.touched_sets(sets, new)) == [1, 2]


def test_scoring_can_only_improve_an_outcome(real_setup):
    """Extra windows add information; no system may come out worse."""
    w, sets = real_setup
    sample = sets.sample(25, weights=sets["weight"], random_state=3)
    periods = np.geomspace(20.0, 400.0, 8)
    phases = np.linspace(0, 1, 5, endpoint=False)

    base = experiment.run_grid(periods, phases, sample, w, duration=DUR)
    last = int(w["sector"].max())
    future = mission.synthesize_windows(
        w, np.arange(last + 1, last + 14), np.arange(1, 14))
    scored = experiment.score_strategy(base, w, future, sample)

    assert (scored["flag_new"] >= scored["flag"]).all()
    assert (scored["new_transits"] >= 0).all()
    assert (scored["flag_new"] != al.ERROR).all()

    # more data can only shrink a surviving alias set, never grow it
    both = (scored["flag"] >= al.AMBIGUOUS) & (scored["flag_new"] >= al.AMBIGUOUS)
    assert (scored.loc[both, "n_surviving_new"]
            <= scored.loc[both, "n_surviving"]).all()

    # and the injected period must still be there afterwards
    survivors = scored[scored["flag_new"] >= al.AMBIGUOUS]
    assert len(survivors) > 20
    for _, r in survivors.iterrows():
        assert np.min(np.abs(np.asarray(r["surviving_new"]) - r["period"])) < 0.1


def test_transition_matrix_is_upper_triangular_and_conserves_systems(real_setup):
    w, sets = real_setup
    sample = sets.sample(15, weights=sets["weight"], random_state=5)
    base = experiment.run_grid(np.geomspace(25.0, 300.0, 6),
                               np.linspace(0, 1, 4, endpoint=False),
                               sample, w, duration=DUR)
    last = int(w["sector"].max())
    future = mission.synthesize_windows(
        w, np.arange(last + 1, last + 14), np.arange(14, 27))
    scored = experiment.score_strategy(base, w, future, sample)

    table = experiment.transition_matrix(scored, weighted=False)
    assert table.to_numpy().sum() == len(scored)
    # ordering is missed < mono < ambiguous < solved, so nothing below the
    # diagonal: an outcome can never degrade.
    order = ["no transit", "mono", "ambiguous", "solved"]
    for i, before in enumerate(order):
        if before not in table.index:
            continue
        for after in order[:i]:
            assert table.loc[before, after] == 0


def test_single_sector_completion_is_normalised():
    periods = np.arange(14.0, 300.0, 1.0)
    solved, mono, missed = targets.completion_for_single_sectors(periods)
    assert np.allclose(solved + mono + missed, 1.0)
    assert (missed[periods > 100] > 0.5).all()
