"""Period aliasing for transits recovered from sparsely sampled photometry.

Two transits separated by a gap ``dt`` are consistent with any period ``dt/n``.
This module enumerates those aliases and eliminates the ones the data exclude.
"""

from dataclasses import dataclass, field

import numpy as np

from .transits import transit_times_from_tc

NO_TRANSIT = 0
MONO = 1
AMBIGUOUS = 2
SOLVED = 3
ERROR = -100

FLAG_NAMES = {
    NO_TRANSIT: "no transit",
    MONO: "mono",
    AMBIGUOUS: "ambiguous",
    SOLVED: "solved",
    ERROR: "error",
}

# Two candidate periods closer than this (days) are treated as the same period.
PERIOD_MATCH_TOL = 0.1


@dataclass
class AliasResult:
    flag: int
    observed_times: np.ndarray = field(default_factory=lambda: np.array([]))
    alias_numbers: np.ndarray = field(default_factory=lambda: np.array([], int))
    alias_periods: np.ndarray = field(default_factory=lambda: np.array([]))
    ruled_out: np.ndarray = field(default_factory=lambda: np.array([], bool))

    @property
    def surviving(self):
        return self.alias_periods[~self.ruled_out]

    @property
    def n_surviving(self):
        return int((~self.ruled_out).sum())

    @property
    def tc(self):
        return self.observed_times[0] if len(self.observed_times) else np.nan


def find_all_aliases(observed_times, min_period=13.0):
    """Enumerate candidate periods consistent with the observed transit times.

    The longest possible period is the smallest gap between consecutive
    observed transits; every integer subdivision of that gap down to
    ``min_period`` is also a candidate.
    """
    observed_times = np.asarray(observed_times)
    if len(observed_times) < 2:
        return np.array([], int), np.array([])

    max_alias = np.min(np.diff(np.sort(observed_times)))
    n_max = max(int(np.floor(max_alias / min_period)), 1)
    numbers = np.arange(1, n_max + 1)
    return numbers, max_alias / numbers


def rule_out_aliases(observed_times, alias_periods, empty_windows, duration=0.0,
                     edge_tol=0.01):
    """Flag aliases the data exclude.

    An alias dies if it either predicts a detectable transit inside a window
    where nothing was seen, or fails to reproduce every transit that *was*
    seen.
    """
    alias_periods = np.atleast_1d(alias_periods)
    ruled_out = np.zeros(len(alias_periods), dtype=bool)
    if len(alias_periods) == 0 or len(observed_times) == 0:
        return ruled_out

    t0 = observed_times[0]
    half = duration / 2.0
    lo = empty_windows["t_start"].values + half
    hi = empty_windows["t_end"].values - half
    usable = hi - lo

    for i, period in enumerate(alias_periods):
        # Would this period have put a transit in a window that came up empty?
        phase_into_window = (t0 - lo) % period
        if np.any((phase_into_window < usable) & (phase_into_window > edge_tol)):
            ruled_out[i] = True
            continue

        # Does this period actually reproduce every transit we did see?
        cycles = np.round((observed_times - t0) / period)
        residual = np.abs(observed_times - (t0 + cycles * period))
        if np.any(residual > PERIOD_MATCH_TOL):
            ruled_out[i] = True

    return ruled_out


def analyze(observed_times, empty_windows, min_period=13.0, duration=0.0):
    """Classify a system and return its surviving period aliases."""
    observed_times = np.asarray(observed_times)
    n = len(observed_times)

    if n == 0:
        return AliasResult(flag=NO_TRANSIT)
    if n == 1:
        return AliasResult(flag=MONO, observed_times=observed_times)

    numbers, periods = find_all_aliases(observed_times, min_period=min_period)
    ruled_out = rule_out_aliases(observed_times, periods, empty_windows,
                                 duration=duration)

    n_surviving = int((~ruled_out).sum())
    if n_surviving == 1:
        flag = SOLVED
    elif n_surviving > 1:
        flag = AMBIGUOUS
    else:
        # The true period is always among the candidates, so this cannot
        # happen unless the transit grid and the windows disagree.
        flag = ERROR

    return AliasResult(flag=flag, observed_times=observed_times,
                       alias_numbers=numbers, alias_periods=periods,
                       ruled_out=ruled_out)


def update_with_windows(result, true_period, new_windows, duration=0.0):
    """Fold additional observing windows into an existing alias solution.

    This is the primitive for scoring a proposed future strategy: it asks, for
    each surviving alias, whether the new windows would distinguish it from the
    true period. Only ``AMBIGUOUS`` systems can be improved.
    """
    if result.flag != AMBIGUOUS:
        return result

    ruled_out = result.ruled_out.copy()
    tc = result.tc
    half = duration / 2.0
    lo = new_windows["t_start"].values + half
    hi = new_windows["t_end"].values - half

    truth = _windows_hit(true_period, tc, lo, hi)

    for i, period in enumerate(result.alias_periods):
        if ruled_out[i]:
            continue
        predicted = _windows_hit(period, tc, lo, hi)
        if not np.array_equal(truth, predicted):
            ruled_out[i] = True

    n_surviving = int((~ruled_out).sum())
    flag = SOLVED if n_surviving == 1 else (AMBIGUOUS if n_surviving > 1 else ERROR)

    return AliasResult(flag=flag, observed_times=result.observed_times,
                       alias_numbers=result.alias_numbers,
                       alias_periods=result.alias_periods, ruled_out=ruled_out)


def _windows_hit(period, tc, lo, hi):
    """Which of the windows [lo, hi] contain a transit of this ephemeris."""
    phase_into_window = (tc - lo) % period
    return phase_into_window <= (hi - lo)
