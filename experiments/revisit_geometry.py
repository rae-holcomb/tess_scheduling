"""Revisit geometry: does *when* you return matter more than *how long* you wait?

A sky-free toy model. One system, two observed transits separated by a gap
``dt``, and one or more future observing windows placed at offset ``x * dt``
after the second transit. Nothing here touches sector footprints, star
catalogues or the real schedule -- the only variables are the arithmetic of the
gap and the placement.

The central object is the **signature** of an alias: the set of transit times it
predicts to be *fully detectable* inside the candidate windows. An alias
survives the new observation if and only if its signature exactly matches the
signature of whichever alias is true --

  * a signature entry the truth lacks  -> predicted a transit that was not seen
  * a truth entry the signature lacks  -> failed to reproduce an observed transit

which are precisely the two rule-out criteria in ``tess_window.aliases``.

So a candidate window does nothing more or less than **partition the ladder into
equivalence classes**. If the truth lies in class ``c``, exactly ``|c|`` aliases
survive. Under a flat prior over the ladder this gives, with no need to loop
over which alias is true:

    E[surviving]     = sum_c |c|^2 / N
    P(solved)        = #{c : |c| = 1} / N
    E[bits left]     = sum_c |c| log2|c| / N

Adding a second window refines the partition. The minimum number of windows
needed to drive every class to size 1 is the "how many more observations?"
question, answered exactly rather than by simulation.

Run with the codeastro env:
    /opt/anaconda3/envs/codeastro/bin/python experiments/revisit_geometry.py
"""

from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

# ---------------------------------------------------------------- parameters

SECTOR = 27.0        # width of one TESS sector, days
TOL = 0.1            # PERIOD_MATCH_TOL from tess_window.aliases
FIG = "figures"

# Fiducial system: a "biennial duotransit" -- two transits two years apart.
DT_FIDUCIAL = 730.0

GOLDEN = (np.sqrt(5.0) - 1.0) / 2.0   # 0.6180..., the most irrational number


def duration(period):
    """Central transit of a solar-density host. D6 in METHODS.md."""
    return (13.0 / 24.0) * (period / 365.0) ** (1.0 / 3.0)


# ------------------------------------------------------------------ the core

def times_in_window(period, t0, w_start, w_end, dur):
    """Transit times of `period` (anchored at t0) that fit *entirely* in a window."""
    lo = w_start + dur / 2.0
    hi = w_end - dur / 2.0
    if hi <= lo:
        return np.empty(0)
    k_lo = int(np.ceil((lo - t0) / period))
    k_hi = int(np.floor((hi - t0) / period))
    if k_hi < k_lo:
        return np.empty(0)
    return t0 + np.arange(k_lo, k_hi + 1) * period


def signature(period, windows, dur):
    """Discretised set of detectable transit times across all windows."""
    out = []
    for w_start, w_end in windows:
        t = times_in_window(period, 0.0, w_start, w_end, dur)
        out.extend(np.round(t / TOL).astype(np.int64).tolist())
    return tuple(out)


def viable_ladder(dt, w_orig=SECTOR, min_period=5.0):
    """Aliases consistent with seeing *exactly two* transits, in the two original windows.

    The two discovery windows are centred on the two transits. An alias whose
    period is short enough to have put a second transit inside one of those
    windows is already dead -- which is where the ~13 d floor of D7 comes from,
    as an emergent consequence of the sector width rather than an assumption.
    """
    win_a = (-w_orig / 2.0, w_orig / 2.0)
    win_b = (dt - w_orig / 2.0, dt + w_orig / 2.0)
    target = (0, int(np.round(dt / TOL)))
    keep = []
    for n in range(1, int(dt // min_period) + 1):
        if signature(dt / n, [win_a, win_b], duration(dt / n)) == target:
            keep.append(n)
    return np.array(keep, dtype=int)


def partition(dt, alias_n, windows):
    """Group aliases by signature. Returns a list of arrays of alias numbers."""
    groups = defaultdict(list)
    for n in alias_n:
        period = dt / n
        groups[signature(period, windows, duration(period))].append(n)
    return [np.array(v) for v in groups.values()]


def metrics(classes, n_total):
    """Expected survivors, P(solved) and expected remaining entropy, flat prior."""
    sizes = np.array([len(c) for c in classes], dtype=float)
    return (
        float((sizes ** 2).sum() / n_total),
        float(sizes[sizes == 1].sum() / n_total),
        float((sizes * np.log2(sizes)).sum() / n_total),
    )


def new_window(dt, x, width=SECTOR):
    """A candidate observing window *centred* x*dt after the second transit.

    Centre rather than start, because the degeneracies land on integer x that
    way and the plots become readable. x = 1 means "centred one full gap after
    the second transit", i.e. on the time the longest-period alias predicts.
    """
    centre = dt + x * dt
    return (centre - width / 2.0, centre + width / 2.0)


def sweep(dt, alias_n, xs, width=SECTOR):
    """Expected survivors / P(solved) / bits, as a function of placement x."""
    out = np.empty((len(xs), 3))
    for i, x in enumerate(xs):
        out[i] = metrics(partition(dt, alias_n, [new_window(dt, x, width)]),
                         len(alias_n))
    return out


# ------------------------------------------------------------------- figures

def fig_setup(dt, alias_n, path):
    """Explain the experiment: the ladder, the placement, and the signature."""
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.5, 1.4], hspace=0.45)
    x_demo = GOLDEN

    # --- panel 1: the timeline -------------------------------------------
    ax = fig.add_subplot(gs[0])
    for w0, w1 in [(-SECTOR / 2, SECTOR / 2), (dt - SECTOR / 2, dt + SECTOR / 2)]:
        ax.add_patch(Rectangle((w0, 0.15), w1 - w0, 0.7, color="0.75", zorder=1))
    c0, c1 = new_window(dt, x_demo)
    ax.add_patch(Rectangle((c0, 0.15), c1 - c0, 0.7, color="#d95f02", alpha=0.8,
                           zorder=1))
    for t in (0.0, dt):
        ax.plot([t, t], [0.15, 0.85], color="k", lw=2.5, zorder=3)
    ax.annotate("", xy=(0, 1.02), xytext=(dt, 1.02),
                arrowprops=dict(arrowstyle="<->", color="k", lw=1.2))
    ax.text(dt / 2, 1.10, f"gap  $\\Delta t$ = {dt:.0f} d", ha="center", fontsize=11)
    ax.annotate("", xy=(dt, -0.12), xytext=(c0, -0.12),
                arrowprops=dict(arrowstyle="<->", color="#d95f02", lw=1.2))
    ax.text((dt + c0) / 2, -0.34, f"offset  $x\\,\\Delta t$   ($x$ = {x_demo:.3f})",
            ha="center", color="#d95f02", fontsize=11)
    ax.text(0, 0.95, "transit 1", ha="center", va="bottom", fontsize=9)
    ax.text(dt, 0.95, "transit 2", ha="center", va="bottom", fontsize=9)
    ax.text(c1 + 20, 0.5, "candidate\nnew window", color="#d95f02", va="center",
            fontsize=10)
    ax.set_xlim(-120, dt * (1 + x_demo) + 300)
    ax.set_ylim(-0.5, 1.35)
    ax.set_yticks([])
    ax.set_xlabel("time (days)")
    ax.set_title("1 · The setup — two transits, then one new window placed at $x\\,\\Delta t$",
                 loc="left", fontsize=12, fontweight="bold")
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)

    # --- panel 2: the alias ladder ---------------------------------------
    ax = fig.add_subplot(gs[1])
    show = alias_n[alias_n <= 14]
    for row, n in enumerate(show):
        period = dt / n
        t = np.arange(0, dt * (1 + x_demo) + SECTOR, period)
        ax.plot(t, np.full_like(t, row), ".", color="0.4", ms=4)
        ax.plot([0, dt], [row, row], "|", color="k", ms=11, mew=2)
        ax.text(-60, row, f"$\\Delta t/{n}$", ha="right", va="center", fontsize=8)
    ax.add_patch(Rectangle((c0, -0.7), c1 - c0, len(show) - 0.6,
                           color="#d95f02", alpha=0.18, zorder=0))
    ax.set_xlim(-120, dt * (1 + x_demo) + 300)
    ax.set_ylim(-0.8, len(show) - 0.2)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlabel("time (days)")
    ax.set_title("2 · Every alias $P=\\Delta t/n$ reproduces both transits — they differ "
                 "only in where *else* they predict one",
                 loc="left", fontsize=12, fontweight="bold")
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)

    # --- panel 3: the signature inside the new window --------------------
    ax = fig.add_subplot(gs[2])
    classes = partition(dt, alias_n, [new_window(dt, x_demo)])
    colour_of = {}
    palette = plt.cm.tab20(np.linspace(0, 1, 20))
    for i, c in enumerate(classes):
        for n in c:
            colour_of[n] = palette[i % 20] if len(c) > 1 else (0.15, 0.15, 0.15, 1.0)
    for row, n in enumerate(alias_n):
        period = dt / n
        t = times_in_window(period, 0.0, c0, c1, duration(period))
        ax.plot([c0, c1], [row, row], "-", color="0.9", lw=1, zorder=0)
        if len(t):
            ax.plot(t, np.full_like(t, row), "o", color=colour_of[n], ms=5, zorder=2)
    ax.set_xlim(c0 - 1, c1 + 1)
    ax.set_ylim(-1, len(alias_n))
    ax.invert_yaxis()
    ax.set_ylabel("alias  $n$")
    ax.set_xlabel("time inside the new window (days)")
    ax.set_title("3 · The signature — coloured = shares its pattern with another alias "
                 "(indistinguishable);  black = unique (resolved)",
                 loc="left", fontsize=12, fontweight="bold")
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)

    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_partitions(dt, alias_n, path):
    """Contrast placements: rational x collapses the ladder, irrational x shatters it."""
    cases = [(1.0, "$x = 1$\n(centred one gap later)"),
             (0.5, "$x = 1/2$"),
             (1.0 / 3.0, "$x = 1/3$"),
             (GOLDEN, "$x = 0.618$\n(golden ratio)")]
    fig, axes = plt.subplots(len(cases), 1, figsize=(13, 7.5), sharex=True)
    n_total = len(alias_n)
    for ax, (x, label) in zip(axes, cases):
        classes = sorted(partition(dt, alias_n, [new_window(dt, x)]),
                         key=len, reverse=True)
        exp_s, p_solved, bits = metrics(classes, n_total)
        palette = plt.cm.tab20(np.linspace(0, 1, 20))
        pos = {n: i for i, n in enumerate(alias_n)}
        for i, c in enumerate(classes):
            col = palette[i % 20] if len(c) > 1 else (0.2, 0.2, 0.2, 1.0)
            for n in c:
                ax.add_patch(Rectangle((pos[n] - 0.5, 0), 1, 1, color=col))
        ax.set_xlim(-0.5, n_total - 0.5)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=10)
        ax.text(1.005, 0.5,
                f"$E$[survivors] = {exp_s:4.1f}   "
                f"$P$(solved) = {p_solved:4.0%}   "
                f"{bits:.2f} bits left",
                transform=ax.transAxes, va="center", fontsize=9, family="monospace")
        for s in ("left", "right", "top", "bottom"):
            ax.spines[s].set_visible(False)
    axes[-1].set_xlabel("alias $n$   —   coloured block = a set of mutually indistinguishable "
                        "aliases;   dark grey = uniquely resolved")
    fig.suptitle("One window, four placements — the same 27-day observation buys "
                 "wildly different information",
                 fontsize=13, fontweight="bold", x=0.09, ha="left", y=1.0)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_sweep(dt, alias_n, xs, res, path):
    """The money plot: expected survivors against placement."""
    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=True,
                             gridspec_kw=dict(height_ratios=[2, 1], hspace=0.12))
    n_total = len(alias_n)

    ax = axes[0]
    ax.plot(xs, res[:, 0], lw=0.7, color="#1b3f8b", zorder=3)
    ax.axhline(n_total, color="0.6", ls=":", lw=1)
    ax.text(xs[0], n_total, "no information gained  ", va="bottom", ha="left",
            fontsize=9, color="0.4")
    ax.set_ylim(0, n_total * 1.16)
    for p, q in [(1, 1), (2, 1), (1, 2), (3, 2), (1, 3), (2, 3), (4, 3), (5, 3),
                 (1, 4), (3, 4), (5, 4), (7, 4)]:
        x = p / q
        if xs[0] <= x <= xs[-1]:
            ax.axvline(x, color="#d95f02", lw=0.8, alpha=0.35, zorder=0)
            ax.text(x, n_total * 1.02, f"$\\frac{{{p}}}{{{q}}}$", ha="center",
                    va="bottom", fontsize=9, color="#d95f02")
    ax.set_ylabel("expected surviving aliases")
    ax.set_title("Expected surviving aliases vs. where the new window is placed  "
                 f"($\\Delta t$ = {dt:.0f} d, {n_total} aliases, one {SECTOR:.0f} d window)",
                 loc="left", fontsize=12, fontweight="bold")
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)

    ax = axes[1]
    ax.plot(xs, res[:, 1] * 100, lw=0.7, color="#0b7d5f")
    ax.set_ylabel("P(solved)  [%]")
    ax.set_xlabel("placement  $x$  =  (window centre $-$ transit 2) / $\\Delta t$")
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)

    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_width(dt, alias_n, path, w_disc=SECTOR):
    """Why a *longer* revisit window breaks the degeneracy — the S97/98 question.

    A revisit window of the same width as the discovery window, centred on a
    multiple of dt, resolves nothing at all (see the docstring theorem). The
    fix is not to wait longer, it is to stare longer: only a window wider than
    the discovery window can fit a second transit of the shortest alias.
    """
    n_total = len(alias_n)
    widths = np.linspace(w_disc * 0.6, w_disc * 3.2, 200)
    centred = [metrics(partition(dt, alias_n,
                                 [(2 * dt - w / 2, 2 * dt + w / 2)]), n_total)[0]
               for w in widths]
    generic = [metrics(partition(dt, alias_n,
                                 [new_window(dt, GOLDEN, w)]), n_total)[0]
               for w in widths]

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))

    ax = axes[0]
    ax.plot(widths, centred, lw=2, color="#b2182b",
            label="centred on $2\\Delta t$  (the naive choice)")
    ax.plot(widths, generic, lw=2, color="#1b3f8b",
            label="generic placement ($x=0.618$)")
    ax.axvline(w_disc, color="0.5", ls="--", lw=1)
    ax.text(w_disc, n_total * 0.97, " standard sector", fontsize=9, color="0.35")
    ax.axvline(2 * w_disc, color="#0b7d5f", ls="--", lw=1)
    ax.text(2 * w_disc, n_total * 0.97, " S97/98 (4-orbit)", fontsize=9,
            color="#0b7d5f")
    ax.set_xlabel("width of the revisit window (days)")
    ax.set_ylabel("expected surviving aliases")
    ax.set_title("Staring longer beats waiting longer", loc="left",
                 fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)

    ax = axes[1]
    xs = np.linspace(0.06, 2.0, 3200)
    for w, col, lab in [(w_disc, "#b2182b", f"{w_disc:.0f} d (1 sector)"),
                        (2 * w_disc, "#0b7d5f", f"{2 * w_disc:.0f} d (double)")]:
        ax.plot(xs, sweep(dt, alias_n, xs, width=w)[:, 0], lw=0.7, color=col,
                label=lab)
    ax.set_xlabel("placement  $x$")
    ax.set_ylabel("expected surviving aliases")
    ax.set_title("A double-length window blunts the spikes — but does not remove them",
                 loc="left", fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)

    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_two_windows(dt, alias_n, path, n_grid=260):
    """Two windows: does the pair need to be arithmetically dissimilar?"""
    xs = np.linspace(0.02, 1.5, n_grid)
    grid = np.full((n_grid, n_grid), np.nan)
    n_total = len(alias_n)
    for i, x1 in enumerate(xs):
        for j, x2 in enumerate(xs):
            if x2 <= x1:
                continue
            wins = [new_window(dt, x1), new_window(dt, x2)]
            grid[j, i] = metrics(partition(dt, alias_n, wins), n_total)[0]

    fig, ax = plt.subplots(figsize=(8.6, 7))
    im = ax.imshow(np.log10(grid), origin="lower", cmap="viridis_r", aspect="auto",
                   extent=[xs[0], xs[-1], xs[0], xs[-1]])
    cb = fig.colorbar(im, ax=ax)
    ticks = np.array([1, 2, 3, 5, 10, 20, 40])
    cb.set_ticks(np.log10(ticks))
    cb.set_ticklabels([str(t) for t in ticks])
    cb.set_label("expected surviving aliases   (dark = bad)")
    ax.set_xlabel("placement of window 1,  $x_1$")
    ax.set_ylabel("placement of window 2,  $x_2$")
    ax.set_title("Two new windows — dark ridges are wasted pairs\n"
                 "vertical/horizontal: one window at a degenerate $x$.   "
                 "diagonal: $x_2-x_1$ a simple rational, so the second window "
                 "re-asks the first one's question",
                 loc="left", fontsize=10.5, fontweight="bold")
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return xs, grid


# ---------------------------------------------------------------------- main

def main():
    dt = DT_FIDUCIAL
    alias_n = viable_ladder(dt)
    n_total = len(alias_n)
    shortest = dt / alias_n.max()
    print(f"dt = {dt:.0f} d   viable aliases: {n_total}"
          f"   periods {shortest:.2f} .. {dt:.0f} d")
    print(f"  (the {shortest:.1f} d floor is emergent from the {SECTOR:.0f} d "
          f"discovery window, not assumed)")

    fig_setup(dt, alias_n, f"{FIG}/revisit_setup.png")
    fig_partitions(dt, alias_n, f"{FIG}/revisit_partitions.png")

    print("\nzero-information theorem: revisit window of the SAME width as the")
    print("discovery window, centred on k*dt, should leave every alias alive.")
    for d in (365.0, 730.0, 1460.0):
        a = viable_ladder(d)
        for k in (2, 3):
            c = k * d
            e = metrics(partition(d, a, [(c - SECTOR / 2, c + SECTOR / 2)]),
                        len(a))[0]
            status = "OK" if abs(e - len(a)) < 1e-9 else "VIOLATED"
            print(f"  dt={d:6.0f}  k={k}  N={len(a):3d}  E[surv]={e:6.2f}  {status}")

    xs = np.linspace(0.06, 2.0, 9000)
    res = sweep(dt, alias_n, xs)
    fig_sweep(dt, alias_n, xs, res, f"{FIG}/revisit_sweep.png")
    fig_width(dt, alias_n, f"{FIG}/revisit_width.png")

    order = np.argsort(res[:, 0])
    print("\nbest placements (fewest expected survivors):")
    for k in order[:6]:
        print(f"  x = {xs[k]:.4f}   E[surv] = {res[k, 0]:5.2f}   "
              f"P(solved) = {res[k, 1]:5.1%}")
    print("worst placements:")
    for k in order[::-1][:6]:
        print(f"  x = {xs[k]:.4f}   E[surv] = {res[k, 0]:5.2f}   "
              f"P(solved) = {res[k, 1]:5.1%}")

    print("\nnamed placements:")
    for x, lab in [(1.0, "x=1"), (0.5, "x=1/2"), (1.0 / 3, "x=1/3"),
                   (0.25, "x=1/4"), (GOLDEN, "golden"), (2.0, "x=2")]:
        e, p, b = metrics(partition(dt, alias_n, [new_window(dt, x)]), n_total)
        print(f"  {lab:>7}  E[surv] = {e:5.2f}   P(solved) = {p:5.1%}   "
              f"bits left = {b:.2f}")

    fig_two_windows(dt, alias_n, f"{FIG}/revisit_two_windows.png")
    print(f"\nfigures written to {FIG}/")


if __name__ == "__main__":
    main()
