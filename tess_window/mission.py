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

import warnings

import numpy as np
import pandas as pd
from astropy.time import Time

# Refreshed periodically by the TESS project; re-download when extending the
# baseline. Sector numbering here must stay consistent with tess_stars2px --
# see validate_sector_alignment.
ORBIT_TABLE_URL = "https://tess.mit.edu/public/files/TESS_orbit_times.csv"

# Last sector fully observed in the bundled orbit table, as of the 2026-08-20
# refresh. Recorded for reference only -- load_orbit_windows derives this with
# last_complete_sector() so it cannot go stale on the next refresh.
LAST_COMPLETE_SECTOR = 106

WINDOW_COLUMNS = ["sector", "pointing", "t_start", "t_end", "predicted"]


def load_orbit_windows(path="TESS_orbit_times.csv", max_sector="auto",
                       errata=True):
    """Read the measured TESS orbit table into a window table.

    ``max_sector="auto"`` keeps only sectors that have finished observing, so a
    sector in progress cannot contribute a truncated set of segments -- that
    would look like real coverage with a spurious gap at the end. Pass an int
    to cap explicitly, or None to keep everything including planned sectors.

    ``errata=True`` applies the corrections in ``ORBIT_TABLE_ERRATA``; set it
    False to read the published file verbatim.
    """
    df = pd.read_csv(path, comment="#")
    df = df.rename(columns={"Sector": "sector", "Orbit": "orbit"})
    df["sector"] = df["sector"].astype(int)
    df["t_start"] = Time(list(df["Start of Orbit"].values)).jd
    df["t_end"] = Time(list(df["End of Orbit"].values)).jd
    df = df.sort_values("t_start").reset_index(drop=True)

    if errata:
        df = _apply_orbit_table_errata(df)

    df["pointing"] = df["sector"]
    df["predicted"] = False
    df = df[WINDOW_COLUMNS].sort_values("t_start").reset_index(drop=True)

    if isinstance(max_sector, str):
        if max_sector != "auto":
            raise ValueError(f"max_sector must be an int, None or 'auto'; "
                             f"got {max_sector!r}")
        max_sector = last_complete_sector(df)

    if max_sector is not None:
        df = df[df["sector"] <= max_sector].reset_index(drop=True)

    return df


# Known defects in the published orbit table, with the fix and the evidence.
#
# As of the 2026-08-20 file, orbit 215 appears FOUR times under sector 102 --
# the only orbit in 329 rows to appear more than twice -- and two of those rows
# span 12.50 d and 12.13 d, the only segments over 10 d in the 4-segment era
# (mean 6.64 d, sigma 0.68) and each almost exactly two normal segments. Orbits
# 216 and 217 are absent. Reading the four rows chronologically as orbits 215,
# 216 (sector 102) and orbit 217 split in two (sector 103) makes orbit
# numbering contiguous and reproduces the sector boundaries published in MIT's
# own pointing table and at HEASARC -- 102: 03/27-04/21, 103: 04/21-05/17 --
# which the raw file contradicts by holding 102 until 05/04. It also drops the
# tess_stars2px midtime offsets from +7.12/+7.41 d to +0.44/+0.70 d, in line
# with sectors 101 (+0.25) and 104 (+0.90). tess_stars2px is NOT at fault here;
# its pointings match MIT's table to 0.000 deg.
#
# Each entry names the rows by (sector, orbit) and is applied only when the
# signature still matches, so it silently becomes a no-op once MIT republishes.
ORBIT_TABLE_ERRATA = [
    {
        "match": {"sector": 102, "orbit": 215},
        "n_rows": 4,
        # per row, in time order: (corrected sector, corrected orbit)
        "relabel": [(102, 215), (102, 216), (103, 217), (103, 217)],
        "note": "sector 102/103 boundary; orbits 216-217 mislabelled as 215",
    },
]


def _apply_orbit_table_errata(df, verbose=False):
    """Correct known mislabelled rows, but only if the defect is still there."""
    df = df.copy()
    for erratum in ORBIT_TABLE_ERRATA:
        sel = np.ones(len(df), dtype=bool)
        for col, val in erratum["match"].items():
            sel &= (df[col] == val).to_numpy()
        idx = df.index[sel]

        if len(idx) != erratum["n_rows"]:
            continue  # already fixed upstream, or the file changed shape

        for row, (sector, orbit) in zip(idx, erratum["relabel"]):
            df.loc[row, "sector"] = sector
            df.loc[row, "orbit"] = orbit
        if verbose:
            print(f"applied orbit-table erratum: {erratum['note']}")
    return df


def last_complete_sector(windows, reference_time=None):
    """Highest sector whose observing has finished by ``reference_time``.

    Derives the cutoff from the table instead of trusting
    ``LAST_COMPLETE_SECTOR``, which goes stale every time the orbit table is
    refreshed. Sectors still in progress or merely planned are excluded, since
    their listed segments may be incomplete.
    """
    if reference_time is None:
        reference_time = Time.now().jd

    bounds = sector_bounds(windows)
    done = bounds.index[bounds["t_end"] < reference_time]
    if len(done) == 0:
        raise ValueError("no sector in the table has finished by "
                         f"JD {reference_time}")

    if bounds["t_end"].max() < reference_time:
        warnings.warn(
            "Every sector in the orbit table is already in the past "
            f"(last ends JD {bounds['t_end'].max():.1f}). The table is stale; "
            f"re-download it from {ORBIT_TABLE_URL}.", stacklevel=2)

    return int(done.max())


def validate_sector_alignment(windows, tolerance=5.0, verbose=False):
    """Check orbit-table sector numbering against ``tess_stars2px``.

    Sky coverage comes from ``tess_stars2px`` and timing comes from the orbit
    table; the two are joined on sector number alone, so a renumbering in
    either source silently pairs a target list with the wrong time windows.
    That happened for real at sectors 97-99 when the extended mission merged
    pointings. Comparing per-sector midtimes catches it.

    Returns a DataFrame indexed by sector with the two midtimes and their
    offset in days, and warns for any sector off by more than ``tolerance``.
    """
    from tess_stars2px import TESS_Spacecraft_Pointing_Data

    bounds = sector_bounds(windows)
    table_mid = 0.5 * (bounds["t_start"] + bounds["t_end"])

    pointing = TESS_Spacecraft_Pointing_Data()
    reference = pd.Series(np.asarray(pointing.midtimes, dtype=float),
                          index=np.asarray(pointing.sectors, dtype=int))

    shared = table_mid.index.intersection(reference.index)
    report = pd.DataFrame({
        "table_midtime": table_mid.loc[shared],
        "pointing_midtime": reference.loc[shared],
    })
    report["offset"] = report["table_midtime"] - report["pointing_midtime"]

    bad = report[report["offset"].abs() > tolerance]
    if len(bad):
        worst = bad["offset"].abs().idxmax()
        warnings.warn(
            f"{len(bad)} sector(s) disagree with tess_stars2px by more than "
            f"{tolerance} d; worst is sector {worst} at "
            f"{report.loc[worst, 'offset']:+.2f} d. Sky coverage and observing "
            "windows may be joined on mismatched sector numbers -- check "
            f"whether {ORBIT_TABLE_URL} has renumbered sectors.", stacklevel=2)
    elif verbose:
        print(f"{len(report)} sectors aligned within {tolerance} d "
              f"(max |offset| = {report['offset'].abs().max():.2f} d)")

    missing = table_mid.index.difference(reference.index)
    if len(missing):
        warnings.warn(
            f"Sectors {sorted(missing)} have observing windows but no "
            "tess_stars2px pointing, so no targets will be assigned to them.",
            stacklevel=2)

    return report


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
