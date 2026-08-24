"""TESS window function for long-period transiting planets.

Typical use::

    from tess_window import mission, targets, experiment

    windows = mission.load_orbit_windows()          # sectors 1..last complete
    mission.validate_sector_alignment(windows)      # sector numbering sanity
    last = int(windows["sector"].max())
    sets = targets.unique_sector_sets(targets.load_targets(max_sector=last))
    results = experiment.run_grid(periods, phases, sets, windows, duration=4/24)
    summary = experiment.summarize(results)

    # score a candidate extended mission: re-observe the year-1 fields
    future = mission.synthesize_windows(windows, np.arange(last + 1, last + 14),
                                        pointings=np.arange(1, 14))
    scored = experiment.score_strategy(results, windows, future, sets)
    experiment.transition_matrix(scored)             # before -> after, in stars
"""

from . import aliases, experiment, mission, plots, targets, transits
from .aliases import AMBIGUOUS, ERROR, MONO, NO_TRANSIT, SOLVED, AliasResult

__all__ = ["aliases", "experiment", "mission", "plots", "targets", "transits",
           "AliasResult", "NO_TRANSIT", "MONO", "AMBIGUOUS", "SOLVED", "ERROR"]
