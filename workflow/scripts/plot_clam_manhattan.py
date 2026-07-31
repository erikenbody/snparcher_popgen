# Send stdout and stderr to the rule's log file. Snakemake does this for R
# scripts (via sink) but not for python ones, which otherwise leaves the log
# empty and the traceback buried in the master log.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

#!/usr/bin/env python3
"""Manhattan plot of windowed Fst (top) and dxy (bottom) for each population
pair in a clam run.

One wide figure per pair: two stacked panels (Fst over dxy), windows laid out
genome-wide ordered by scaffold size (largest first), with alternating
light/dark shading per scaffold. Scaffold sizes come from the reference .fai.

Fst/dxy only exist when the run has >=2 populations; with a single population
this writes the done sentinel and produces no plots (header-only pairwise
summary). With >2 populations clam emits one (population1, population2) pair per
window-row group, and we produce one figure per pair.

Run as a Snakemake `script:` — expects snakemake.input/output/params.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DARK = "#1f4e79"
LIGHT = "#7fb1d6"


def _load(path, value_col):
    """Load a clam window TSV, keeping chrom/start/end/value. Returns None if
    the file is absent (single-population run)."""
    if not path or not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t")
    keep = ["chrom", "start", "end", "population1", "population2", value_col]
    df = df[[c for c in keep if c in df.columns]].copy()
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")
    df["pos"] = (df["start"] + df["end"]) / 2.0
    return df


def _pair_label(p1, p2):
    """Strip clam's 'pop' prefix for display/filenames."""
    s1 = str(p1)[3:] if str(p1).startswith("pop") else str(p1)
    s2 = str(p2)[3:] if str(p2).startswith("pop") else str(p2)
    return s1, s2


def make_plot(prefix, title, sub_fst, sub_dxy, contig_len, out_png):
    chroms = sorted(
        set(sub_fst["chrom"]) | set(sub_dxy["chrom"]),
        key=lambda c: contig_len.get(c, 0), reverse=True,
    )
    offset, ticks, tick_labels = {}, [], []
    cum = 0
    for i, c in enumerate(chroms):
        offset[c] = cum
        L = contig_len.get(c, 0)
        ticks.append(cum + L / 2.0)
        tick_labels.append(str(i + 1))
        cum += L
    genome_len = cum
    color_of = {c: (DARK if i % 2 == 0 else LIGHT) for i, c in enumerate(chroms)}

    def xy(df, metric):
        s = df.dropna(subset=[metric])
        x = s["pos"].values + s["chrom"].map(offset).values
        return x, s[metric].values, s["chrom"].map(color_of).values

    fig, (ax_fst, ax_dxy) = plt.subplots(
        2, 1, figsize=(22, 7), sharex=True, gridspec_kw=dict(hspace=0.12))

    xf, yf, cf = xy(sub_fst, "fst")
    ax_fst.scatter(xf, yf, c=cf, s=4, edgecolors="none", rasterized=True)
    ax_fst.set_ylabel("Fst (Hudson)")
    if len(yf):
        ax_fst.set_ylim(bottom=min(0, np.nanmin(yf)))
    ax_fst.set_title(f"{prefix}  —  {title}", fontsize=13, loc="left")

    xd, yd, cd = xy(sub_dxy, "dxy")
    ax_dxy.scatter(xd, yd, c=cd, s=4, edgecolors="none", rasterized=True)
    ax_dxy.set_ylabel("dxy")
    ax_dxy.set_ylim(bottom=0)

    for ax in (ax_fst, ax_dxy):
        ax.set_xlim(0, genome_len)
        ax.margins(x=0)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    ax_dxy.set_xticks(ticks)
    ax_dxy.set_xticklabels(tick_labels, fontsize=6, rotation=90)
    ax_dxy.set_xlabel("Scaffold (ordered by size, largest = 1)")

    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_png.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"  {title}: {len(chroms)} scaffolds, {len(xf)} fst / {len(xd)} dxy "
          f"windows -> {out_png}")


def main():
    snk = snakemake  # noqa: F821 (injected by Snakemake)
    prefix = snk.params.prefix
    outdir = snk.params.outdir
    fai = snk.input.fai

    contig_len = pd.read_csv(
        fai, sep="\t", header=None, usecols=[0, 1],
        names=["chrom", "length"]).set_index("chrom")["length"].to_dict()

    fst = _load(snk.params.fst, "fst")
    dxy = _load(snk.params.dxy, "dxy")

    if fst is None and dxy is None:
        print(f"[{prefix}] single-population run (no fst/dxy) -- no Manhattan plot")
        return

    # Enumerate population pairs present in either metric.
    def pairs_of(df):
        if df is None or "population1" not in df.columns:
            return set()
        return set(map(tuple, df[["population1", "population2"]].drop_duplicates().values))

    pairs = sorted(pairs_of(fst) | pairs_of(dxy))
    print(f"[{prefix}] {len(pairs)} population pair(s)")

    def slice_pair(df, p1, p2):
        if df is None:
            return pd.DataFrame(columns=["chrom", "pos"])
        if "population1" in df.columns:
            return df[(df["population1"] == p1) & (df["population2"] == p2)]
        return df

    for p1, p2 in pairs:
        s1, s2 = _pair_label(p1, p2)
        title = f"{s1} vs {s2}"
        suffix = "" if len(pairs) == 1 else f"_{s1}_vs_{s2}"
        out_png = os.path.join(outdir, f"{prefix}_manhattan_fst_dxy{suffix}.png")
        make_plot(prefix, title,
                  slice_pair(fst, p1, p2), slice_pair(dxy, p1, p2),
                  contig_len, out_png)


main()
