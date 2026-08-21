"""Diagnostic plots for single systems and for batch results."""

import matplotlib.pyplot as plt
import numpy as np

from . import aliases as al
from .transits import transit_times_from_tc

OUTCOME_COLORS = {"solved": "tab:green", "ambiguous": "tab:orange",
                  "mono": "tab:blue", "missed": "tab:red"}


def plot_windows(ax, visible, window_obs, t_offset=0.0):
    """Shade observing windows green if they caught a transit, red if not."""
    for hit, (_, w) in zip(window_obs, visible.iterrows()):
        ax.axvspan(w["t_start"] - t_offset, w["t_end"] - t_offset,
                   alpha=0.35, lw=0, color="green" if hit else "red")


def plot_transits(transit_times, transit_obs, visible, window_obs,
                  t_offset=2450000, title=None, ax=None):
    """Timeline of injected transits against the observing windows."""
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 3))

    plot_windows(ax, visible, window_obs, t_offset)
    tt = np.asarray(transit_times) - t_offset
    ax.scatter(tt, np.ones_like(tt), marker="|", c="grey", label="missed transit")
    ax.scatter(tt[transit_obs], np.ones_like(tt[transit_obs]), marker="o", s=18,
               c="k", label="observed transit")
    ax.set_yticks([])
    ax.set_xlabel(f"BJD - {t_offset:.0f}")
    ax.legend(loc="upper right", fontsize=8)
    if title:
        ax.set_title(title)
    return ax


def plot_aliases(result, visible, window_obs, stop_time, t_offset=2450000,
                 title=None, ax=None):
    """One row per candidate period, showing which the windows rule out."""
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 8))

    plot_windows(ax, visible, window_obs, t_offset)
    start_time = visible["t_start"].min()

    for n, period, dead in zip(result.alias_numbers, result.alias_periods,
                               result.ruled_out):
        times = transit_times_from_tc(period, result.tc, start_time, stop_time)
        ax.scatter(times - t_offset, np.full_like(times, n), marker="|", s=6,
                   c="red" if dead else "blue")

    ax.set_xlabel(f"BJD - {t_offset:.0f}")
    ax.set_ylabel("alias number")
    ax.set_title(title or f"{result.n_surviving} surviving period(s)")
    return ax


def plot_completion(summary, x="period", ax=None, title="TESS completion"):
    """Stacked bar chart of outcome fractions against period."""
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))

    bottom = np.zeros(len(summary))
    for name in ["solved", "ambiguous", "mono", "missed"]:
        ax.bar(summary[x], summary[name], bottom=bottom,
               width=np.diff(summary[x]).min() * 0.9 if len(summary) > 1 else 1,
               label=name, color=OUTCOME_COLORS[name])
        bottom += summary[name].values

    ax.set_xlabel("period (d)")
    ax.set_ylabel("fraction of injected planets")
    ax.set_ylim(0, 1)
    ax.set_title(title)
    ax.legend(loc="lower right", fontsize=8)
    return ax
