"""Transit ephemerides and which transits land inside observing windows."""

import numpy as np

HOURS = 1.0 / 24.0

# Central-transit duration of a solar-type host at P = 365 d.
SOLAR_DURATION_365 = 13.0 * HOURS


def transit_duration(period, stellar_density=1.0):
    """Central-transit duration in days for a host of given mean density.

    Follows T ~ (P/pi) * arcsin(R*/a), which for a central transit reduces to
    a P^(1/3) scaling at fixed stellar density. Long-period planets have
    proportionally longer transits, which matters when deciding whether an
    event fits inside an observing window.
    """
    period = np.asarray(period, dtype=float)
    return SOLAR_DURATION_365 * (period / 365.0) ** (1 / 3) * stellar_density ** (-1 / 3)


def resolve_duration(duration, period):
    """Allow ``duration`` to be a constant, or any callable of period."""
    return float(duration(period)) if callable(duration) else float(duration)


def transit_times_from_phase(period, phases, start_time, stop_time):
    """Transit times for planets whose conjunction falls at a fractional phase.

    Returns a (n_phases, n_transits) array. The transit count is derived from
    the requested [start_time, stop_time] window so the grid always spans the
    whole window regardless of which strategy is being tested.
    """
    phases = np.atleast_1d(phases).astype(float)
    n = int((stop_time - start_time) // period + 1)
    return (np.arange(n)[None, :] + phases[:, None]) * period + start_time


def transit_times_from_tc(period, tc, start_time, stop_time):
    """Transit times for a planet with a known conjunction time.

    ``tc`` need not be the first transit or lie inside the window.
    """
    first = start_time + ((tc - start_time) % period)
    return np.arange(first, stop_time, period)


def check_observability(transit_times, target_sectors, windows, duration=0.0):
    """Match transits against the windows in which a target is observable.

    A transit counts as observed only if the *entire* event, ``tc +/-
    duration/2``, falls inside a single window, so events straddling a downlink
    gap or a sector edge are correctly excluded.

    Args:
        transit_times: (n_transits,) or (n_phases, n_transits) array of BJD.
        target_sectors: sectors in which this target's sky position is observed.
        windows: window table (see ``mission``).
        duration: transit duration in days.

    Returns:
        transit_obs: (n_phases, n_transits) bool, was each transit caught.
        window_obs:  (n_phases, n_windows) bool, did each window catch one.
        visible: the subset of ``windows`` covering this target, row-aligned
            with the second axis of ``window_obs``.
    """
    tt = np.atleast_2d(transit_times)
    visible = windows[windows["pointing"].isin(np.asarray(target_sectors))]

    half = duration / 2.0
    lo = (visible["t_start"].values + half)[None, None, :]
    hi = (visible["t_end"].values - half)[None, None, :]
    hit = (tt[:, :, None] >= lo) & (tt[:, :, None] <= hi)

    return hit.any(axis=2), hit.any(axis=1), visible
