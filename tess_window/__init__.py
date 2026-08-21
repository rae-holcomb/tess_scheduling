"""TESS window function for long-period transiting planets.

Typical use::

    from tess_window import mission, targets, experiment

    windows = mission.load_orbit_windows()
    sets = targets.unique_sector_sets(targets.load_targets(max_sector=98))
    results = experiment.run_grid(periods, phases, sets, windows, duration=4/24)
    summary = experiment.summarize(results)
"""

from . import aliases, experiment, mission, plots, targets, transits
from .aliases import AMBIGUOUS, ERROR, MONO, NO_TRANSIT, SOLVED, AliasResult

__all__ = ["aliases", "experiment", "mission", "plots", "targets", "transits",
           "AliasResult", "NO_TRANSIT", "MONO", "AMBIGUOUS", "SOLVED", "ERROR"]
