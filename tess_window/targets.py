"""Target sky positions and the sectors that cover them."""

from ast import literal_eval

import numpy as np
import pandas as pd

EXTENDED_CATALOG = "target_df_shortform_extended_mission.csv"
BASE_CATALOG = "target_df_shortform.csv"


def load_targets(path=EXTENDED_CATALOG, max_sector=None):
    """Load the per-target sector list produced by tess_stars2px."""
    df = pd.read_csv(path)
    df = df.drop(columns=df.columns[0])
    df["sec"] = df["sec"].apply(literal_eval).apply(np.array)

    if max_sector is not None:
        df["sec"] = df["sec"].apply(lambda s: s[s <= max_sector])
        df = df[df["sec"].apply(len) > 0].reset_index(drop=True)

    df["num_s"] = df["sec"].apply(len)
    return df


def unique_sector_sets(targets, min_sectors=1):
    """Collapse targets onto their distinct sector combinations.

    Most of the ~8000 catalog targets share a sector combination with some
    other target, so the window-function calculation only needs to run once per
    combination. ``weight`` carries the number of targets represented, which is
    what makes a sky-averaged completeness meaningful.
    """
    keep = targets[targets["num_s"] >= min_sectors]
    rows = {}
    for _, row in keep.iterrows():
        key = tuple(sorted(row["sec"]))
        if key not in rows:
            rows[key] = {"sectors": np.array(key), "weight": 0,
                         "num_s": len(key), "ec_lat": [], "ec_long": []}
        rows[key]["weight"] += 1
        rows[key]["ec_lat"].append(row["ec_lat"])
        rows[key]["ec_long"].append(row["ec_long"])

    out = pd.DataFrame([
        {"sectors": v["sectors"], "num_s": v["num_s"], "weight": v["weight"],
         "ec_lat": float(np.mean(v["ec_lat"])),
         "ec_long": float(np.mean(v["ec_long"]))}
        for v in rows.values()
    ])
    return out.sort_values("weight", ascending=False).reset_index(drop=True)


def completion_for_single_sectors(periods, sector_length=27.0):
    """Analytic window function for targets observed in exactly one sector.

    Single-sector targets need no simulation: for P below the sector length the
    outcome depends only on where the conjunction falls, so the solved, mono
    and missed fractions are closed-form.
    """
    periods = np.atleast_1d(periods).astype(float)
    solved = np.full(periods.shape, np.nan)
    mono = np.full(periods.shape, np.nan)
    missed = np.full(periods.shape, np.nan)

    short = (periods >= sector_length / 2) & (periods < sector_length)
    solved[short] = (sector_length - periods[short]) / sector_length * 2
    mono[short] = (2 * periods[short] - sector_length) / sector_length
    missed[short] = 0.0

    long = periods >= sector_length
    solved[long] = 0.0
    mono[long] = sector_length / periods[long]
    missed[long] = (periods[long] - sector_length) / periods[long]

    return solved, mono, missed
