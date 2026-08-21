"""Observing windows for the TESS mission.

The fundamental data structure throughout this package is a *window table*: a
DataFrame with one row per contiguous block of observing, and columns

    sector    int    slot label, unique per ~27 day pointing
    pointing  int    which sector's sky footprint applies to this slot
    t_start   float  BJD
    t_end     float  BJD
    predicted bool   True if synthesized rather than measured

For the real mission ``pointing == sector``. They differ only for hypothetical
future strategies that re-use an earlier pointing at a later time.
"""

import numpy as np
import pandas as pd
from astropy.time import Time

# TESS_orbit_times.csv is complete only through this sector; sector 99 is
# present but truncated mid-sector in the current file.
LAST_COMPLETE_SECTOR = 98

WINDOW_COLUMNS = ["sector", "pointing", "t_start", "t_end", "predicted"]


def load_orbit_windows(path="TESS_orbit_times.csv", max_sector=LAST_COMPLETE_SECTOR):
    """Read the measured TESS orbit table into a window table."""
    df = pd.read_csv(path, comment="#")
    df = df.rename(columns={"Sector": "sector", "Orbit": "orbit"})
    df["sector"] = df["sector"].astype(int)
    df["t_start"] = Time(list(df["Start of Orbit"].values)).jd
    df["t_end"] = Time(list(df["End of Orbit"].values)).jd

    if max_sector is not None:
        df = df[df["sector"] <= max_sector]

    df["pointing"] = df["sector"]
    df["predicted"] = False
    return df[WINDOW_COLUMNS].sort_values("t_start").reset_index(drop=True)


def sector_bounds(windows):
    """Collapse a window table to one row per sector with overall start/stop."""
    return windows.groupby("sector").agg(
        pointing=("pointing", "first"),
        t_start=("t_start", "min"),
        t_end=("t_end", "max"),
        n_windows=("t_start", "size"),
    )


def mission_span(windows):
    return windows["t_start"].min(), windows["t_end"].max()


def synthesize_windows(windows, sectors, pointings=None, n_template=20,
                       n_segments=None, gap=None, cadence=None, span=None):
    """Extrapolate observing windows for sectors beyond the measured table.

    Segment count, inter-segment gap, sector duration and sector-to-sector
    cadence all default to the median of the last ``n_template`` measured
    sectors, so the synthesized schedule inherits the current duty cycle.

    ``sectors`` labels the new slots; ``pointings`` says which sector's sky
    footprint each slot re-observes (defaults to ``sectors`` itself).
    """
    sectors = np.atleast_1d(sectors)
    if pointings is None:
        pointings = sectors
    pointings = np.atleast_1d(pointings)
    if len(pointings) != len(sectors):
        raise ValueError("sectors and pointings must be the same length")

    bounds = sector_bounds(windows)
    recent = bounds.iloc[-n_template:]

    if cadence is None:
        cadence = float(np.median(np.diff(recent["t_start"].values)))
    if span is None:
        span = float(np.median(recent["t_end"].values - recent["t_start"].values))
    if n_segments is None:
        n_segments = int(recent["n_windows"].median())
    if gap is None:
        per_sector_gap = [
            _mean_gap(windows[windows["sector"] == s]) for s in recent.index
        ]
        gap = float(np.nanmedian(per_sector_gap))

    # Anchor off the end of the measured schedule rather than the last sector
    # start, so an unusually long final sector cannot overlap the first
    # synthesized one.
    sector_gap = float(np.median(
        bounds["t_start"].values[1:] - bounds["t_end"].values[:-1]))
    first_start = windows["t_end"].max() + sector_gap
    span = min(span, cadence - sector_gap)

    seg_len = (span - gap * (n_segments - 1)) / n_segments

    rows = []
    for i, (sec, point) in enumerate(zip(sectors, pointings)):
        t0 = first_start + cadence * i
        for j in range(n_segments):
            start = t0 + j * (seg_len + gap)
            rows.append((int(sec), int(point), start, start + seg_len, True))

    return pd.DataFrame(rows, columns=WINDOW_COLUMNS)


def _mean_gap(sector_windows):
    if len(sector_windows) < 2:
        return np.nan
    w = sector_windows.sort_values("t_start")
    return np.mean(w["t_start"].values[1:] - w["t_end"].values[:-1])


def extend_windows(windows, sectors, pointings=None, **kwargs):
    """Append synthesized windows to a measured window table."""
    future = synthesize_windows(windows, sectors, pointings=pointings, **kwargs)
    return pd.concat([windows, future], ignore_index=True)
