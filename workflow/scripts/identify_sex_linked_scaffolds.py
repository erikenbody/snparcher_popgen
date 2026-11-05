"""Scaffold QC: which scaffolds to leave out of every per-population analysis.

Input: genotype counts in evenly spaced windows of the analysis set (rule
sex_id_extract): per window and sample, called and heterozygous genotypes and
summed depth. All SNPs in those windows are used, not the LD-pruned set: pruning
thins low-recombination scaffolds harder, which made whole macrochromosomes look
low-het per bp (41-Artemisiospiza).

Per sample and scaffold >= 1 Mb:
  het    heterozygous genotypes per window bp, over the scaffold's windows minus het
         hotspots (windows where samples average > HOT x their usual het: paralogs,
         collapsed repeats)
  depth  median over the scaffold's windows of mean DP
Each is divided by the sample's typical value: its median over scaffolds >= 5 Mb
that carry signal (het >= 0.2 x the sample's 90th percentile). A median over every
scaffold >= 5 Mb sits far below the chromosomes in assemblies with many near-empty
scaffolds (41-Reithrodontomys: 27x), and an X then looks normal.

Scaffolds of at least sex_id_rule_min_size are excluded when they are
  zero het        median normalised het < 0.1: a Y or W, haplotigs
  low het in all  >= 85% of samples below 0.25
  X/Z             a subset of samples is hemizygous, as the heterogametic sex is
                  on its X or Z (or carriers of a large deletion), by either of
                  half depth: >= 5 samples (and >= 5%) below 0.7 depth while >= 5
                    are above 0.85; the half-depth group's median depth <= 0.6 and
                    >= 40% of it below 0.25 het. Subsets at 0.6-0.7 depth with
                    normal het are batches or low-coverage libraries, not
                    hemizygosity (17-Lottia, 91-Azolla)
                  het split: >= 5 samples (and >= 10%) below 0.25 het and >= 5
                    above, the low-het ones at <= 0.9 x the others' median depth
                    (partly sex-linked scaffolds: 41-Melanerpes, low-het group at
                    0.80 vs 0.91). The larger minimum subset keeps a few
                    odd samples from removing a scaffold for everyone
  half depth, low het in most
                  median depth < 0.7 and >= 30% of samples below 0.25 het: a
                  scaffold present once per genome in most samples (58-Sceloporus
                  "sex_linked_1")
and reported but kept when
  hemizygous pattern, full depth   30-85% of samples below 0.25 het (and >= 3x the
                  project's usual share + 0.1) at normal depth: an inversion, an ROH
                  hotspot or an artefact, not a sex chromosome. Look before
                  excluding it with sex_id_exclude.
Scaffolds < 1 Mb are always excluded. Config sex_id_exclude adds scaffolds by name
and sex_id_keep rescues them. An assembly with no scaffold >= 5 Mb is too
fragmented to normalise: nothing is excluded but sex_id_exclude.

This replaces a test that needed two depth clusters (k-means, silhouette > 0.8,
>= 10 samples each) with different het (Welch t, Bonferroni p < 1e-4). It never saw
a Y, W or haplotig (one depth level, zero het in everyone), missed sex chromosomes
with < 10 samples of one sex or noisy depth (102-Colias, 41-Perognathus_ext,
58-Masticophis), and excluded scaffolds whose depth "clusters" differed by 4-15%
(sequencing batches: 77-Linanthus, 91-Azolla, 74-Chrysomela).
"""
# Send stdout and stderr to the rule's log file. Snakemake does this for R
# scripts (via sink) but not for python ones, which otherwise leaves the log
# empty and the traceback buried in the master log.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MIN_LEN = 1_000_000     # shorter scaffolds are always excluded
REF_LEN = 5_000_000     # scaffolds that set each sample's typical het and depth
HOT = 4.0               # het hotspot: window het / sample's mean, averaged over samples
ZERO_HET = 0.1
LOW_HET = 0.25
ALL_LOW = 0.85
HEMI_LOW = 0.30
HALF_DEPTH = 0.70       # a sample below this depth on a scaffold is in its half-depth group
FULL_DEPTH = 0.85
HALF_GROUP_DEPTH = 0.60  # X/Z: the half-depth group's median depth at most this
HALF_GROUP_LOWHET = 0.40  # X/Z: at least this share of the half-depth group below LOW_HET
LOWHET_DEPTH_RATIO = 0.90  # X/Z: low-het samples' median depth / normal-het samples'
LOWHET_SUBSET = 0.10    # X/Z by het split: low-het samples at least this share
MIN_GROUP = 5
XZ = "hemizygous in a subset (X/Z)"
HALF_MOST = "half depth, low het in most"
EXCLUDE_CALLS = ("zero het", "low het in all", HALF_MOST, XZ)

CALL_COLUMNS = [
    "chrom", "length", "call", "excluded", "reason", "n_samples", "het_med",
    "frac_low_het", "depth_med", "depth_lowhet", "depth_normhet", "n_half_depth",
    "n_full_depth", "depth_half_depth", "frac_half_depth_lowhet", "het_half_depth",
    "het_full_depth", "n_lowhet", "n_normhet", "n_windows", "n_hot_windows",
]
DATA_COLUMNS = ["Sample", "Chromosome", "Depth", "Heterozygosity",
                "Normalized_Depth", "Normalized_Heterozygosity"]


def window_ids(contigs, window, step):
    """(chrom, block) of every window, as sex_id_extract lays them out."""
    o = (step - window) // 2
    rows = []
    for chrom, length in contigs.itertuples(index=False):
        if length < MIN_LEN:
            continue
        b = 0
        while b + o + window <= length:
            rows.append((chrom, b))
            b += step
    return pd.DataFrame(rows, columns=["chrom", "block"])


def nanmedian(x, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(x, axis=axis)


def median_or_nan(x):
    return float(np.median(x)) if len(x) else np.nan


def write_outputs(out, exclude, calls, data, plot_fn=None):
    with open(out["exclude"], "w") as f:
        for c in sorted(exclude):
            f.write(f"{c}\n")
    calls = calls.reindex(columns=CALL_COLUMNS)
    calls.to_csv(out["calls"], sep="\t", index=False, na_rep="NA", float_format="%.4g")
    calls[calls["excluded"] == True].to_csv(out["report"], index=False, na_rep="NA",
                                            float_format="%.4g")
    data.reindex(columns=DATA_COLUMNS).to_csv(out["data"], index=False)
    if plot_fn is None:
        open(out["plot"], "w").close()
    else:
        plot_fn(out["plot"])


def main(counts_path, samples_path, contigs_path, out, window, step, rule_min, extra, keep):
    contigs = pd.read_csv(contigs_path, sep="\t", header=None, usecols=[0, 1],
                          names=["chrom", "len"], dtype={"chrom": str})
    length = dict(zip(contigs["chrom"], contigs["len"]))
    extra, keep = set(extra or []), set(keep or [])
    for name, s in (("sex_id_exclude", extra), ("sex_id_keep", keep)):
        missing = sorted(s - set(length))
        if missing:
            print(f"Warning: {name} names {len(missing)} scaffold(s) not in the assembly: "
                  f"{missing[:5]}")
    small = set(contigs.loc[contigs["len"] < MIN_LEN, "chrom"])
    empty_calls = pd.DataFrame(columns=CALL_COLUMNS)
    empty_data = pd.DataFrame(columns=DATA_COLUMNS)

    if not (contigs["len"] >= REF_LEN).any():
        print(f"No scaffold >= {REF_LEN:,} bp: too fragmented to normalise depth and "
              "heterozygosity. Nothing excluded but sex_id_exclude.")
        write_outputs(out, (extra & set(length)) - keep, empty_calls, empty_data)
        return

    samples = [line.strip() for line in open(samples_path) if line.strip()]
    nS = len(samples)
    win = window_ids(contigs, window, step)
    win["w"] = np.arange(len(win))
    cnt = pd.read_csv(counts_path, sep="\t", header=None,
                      names=["chrom", "block", "s", "called", "het", "dp_sum", "dp_n"],
                      dtype={"chrom": str})
    cnt = cnt.merge(win, on=["chrom", "block"], how="inner")
    nW = len(win)
    H = np.zeros((nW, nS))
    D = np.zeros((nW, nS))
    N = np.zeros((nW, nS))
    idx = (cnt["w"].to_numpy(), cnt["s"].to_numpy())
    np.add.at(H, idx, cnt["het"].to_numpy())
    np.add.at(D, idx, cnt["dp_sum"].to_numpy())
    np.add.at(N, idx, cnt["dp_n"].to_numpy())
    del cnt
    print(f"{nS} samples, {nW:,} windows of {window:,} bp every {step:,} bp on "
          f"{win['chrom'].nunique():,} scaffolds >= {MIN_LEN:,} bp")

    # het hotspots: windows where samples, on average, sit far above their usual rate
    wlen = win["chrom"].map(length).to_numpy()
    m = H[wlen >= REF_LEN].mean(axis=0)
    has_het = m > 0
    if has_het.any():
        hot = (H[:, has_het] / m[has_het]).mean(axis=1) > HOT
    else:
        print("No sample has heterozygous genotypes (haploid?): heterozygosity rules skipped.")
        hot = np.zeros(nW, dtype=bool)
    print(f"{hot.sum():,} het hotspot windows ({hot.mean():.1%}) left out")

    codes, scaffolds = pd.factorize(win["chrom"])
    scaffolds = np.asarray(scaffolds)
    nC = len(scaffolds)
    use = ~hot
    n_win = np.bincount(codes, minlength=nC)
    n_use = np.bincount(codes[use], minlength=nC)
    Hs = np.zeros((nC, nS))
    np.add.at(Hs, codes[use], H[use])
    with np.errstate(invalid="ignore", divide="ignore"):
        het_rate = Hs / (n_use[:, None] * window)
    het_rate[n_use == 0] = np.nan
    dep = np.where(N > 0, D / np.maximum(N, 1), np.nan)
    order = np.argsort(codes, kind="stable")
    bounds = np.searchsorted(codes[order], np.arange(nC + 1))
    depth = np.full((nC, nS), np.nan)
    for i in range(nC):
        rows = order[bounds[i]:bounds[i + 1]]
        rows = rows[use[rows]]
        if len(rows):
            depth[i] = nanmedian(dep[rows], axis=0)
    del H, D, N, dep

    # each sample's typical het and depth: scaffolds >= REF_LEN that carry signal
    sc_len = np.array([length[c] for c in scaffolds])
    forced = np.isin(scaffolds, list(extra))
    ref = (sc_len >= REF_LEN) & ~forced
    if ref.sum() < 3:
        ref = ~forced
    h_base = np.full(nS, np.nan)
    d_base = np.full(nS, np.nan)
    for j in range(nS):
        x, y = het_rate[ref, j], depth[ref, j]
        ok = np.isfinite(x)
        if not ok.any():
            continue
        q = np.quantile(x[ok], 0.9)
        good = (ok & (x >= 0.2 * q)) if q > 0 else ok
        if q > 0:
            h_base[j] = np.median(x[good])
        d_base[j] = nanmedian(y[good])
    with np.errstate(invalid="ignore", divide="ignore"):
        nh = het_rate / h_base
        nd = depth / d_base
    nh[:, ~(h_base > 0)] = np.nan

    stats = []
    for i, c in enumerate(scaffolds):
        h, d = nh[i], nd[i]
        hv = h[np.isfinite(h)]
        both = np.isfinite(h) & np.isfinite(d)
        half = both & (d < HALF_DEPTH)
        full = both & (d > FULL_DEPTH)
        stats.append(dict(
            chrom=c, length=int(sc_len[i]), n_samples=len(hv),
            het_med=median_or_nan(hv),
            frac_low_het=float(np.mean(hv < LOW_HET)) if len(hv) else np.nan,
            depth_med=float(nanmedian(d)) if np.isfinite(d).any() else np.nan,
            depth_lowhet=median_or_nan(d[both & (h < LOW_HET)]),
            depth_normhet=median_or_nan(d[both & (h >= LOW_HET)]),
            n_lowhet=int((both & (h < LOW_HET)).sum()),
            n_normhet=int((both & (h >= LOW_HET)).sum()),
            n_half_depth=int(half.sum()), n_full_depth=int(full.sum()),
            depth_half_depth=median_or_nan(d[half]),
            frac_half_depth_lowhet=float(np.mean(h[half] < LOW_HET)) if half.any() else np.nan,
            het_half_depth=median_or_nan(h[half]), het_full_depth=median_or_nan(h[full]),
            n_windows=int(n_win[i]), n_hot_windows=int(n_win[i] - n_use[i])))
    calls = pd.DataFrame(stats)
    usual = calls.loc[(calls["length"] >= REF_LEN) & (calls["het_med"] >= ZERO_HET)
                      & (calls["frac_low_het"] < ALL_LOW), "frac_low_het"]
    bg = float(usual.median()) if len(usual) else 0.0

    def call(r):
        if r.length < rule_min:
            return "not tested"
        if not has_het.any() or r.n_samples == 0:
            return "no data"
        if r.het_med < ZERO_HET:
            return "zero het"
        if r.frac_low_het >= ALL_LOW:
            return "low het in all"
        subset = max(MIN_GROUP, 0.05 * r.n_samples)
        if ((r.n_half_depth >= subset and r.n_full_depth >= MIN_GROUP
                and r.depth_half_depth <= HALF_GROUP_DEPTH
                and r.frac_half_depth_lowhet >= HALF_GROUP_LOWHET)
                or (r.n_lowhet >= max(MIN_GROUP, LOWHET_SUBSET * r.n_samples)
                    and r.n_normhet >= MIN_GROUP
                    and r.depth_lowhet <= LOWHET_DEPTH_RATIO * r.depth_normhet)):
            return XZ
        if r.depth_med < HALF_DEPTH and r.frac_low_het >= HEMI_LOW:
            return HALF_MOST
        if r.frac_low_het >= HEMI_LOW and r.frac_low_het >= 3 * bg + 0.1:
            return "hemizygous pattern, full depth"
        return "normal"

    calls["call"] = [call(r) for r in calls.itertuples(index=False)]
    by_rule = calls["call"].isin(EXCLUDE_CALLS)
    by_config = calls["chrom"].isin(extra)
    kept = calls["chrom"].isin(keep)
    calls["excluded"] = (by_rule | by_config) & ~kept
    calls["reason"] = np.select([kept & (by_rule | by_config), by_config, by_rule],
                                ["sex_id_keep", "sex_id_exclude", "rule"], "")
    exclude = (small | set(calls.loc[calls["excluded"], "chrom"]) | (extra & set(length))) - keep

    print(f"usual share of samples with low het on a scaffold: {bg:.3f}")
    print(calls["call"].value_counts().to_string())
    flagged = calls[calls["call"].isin(EXCLUDE_CALLS + ("hemizygous pattern, full depth",))
                    | by_config | kept]
    with pd.option_context("display.width", 220, "display.max_rows", 500,
                           "display.max_columns", 30):
        print(flagged.drop(columns=["n_windows"]).round(3).to_string(index=False))
    print(f"Excluding {len(exclude):,} scaffolds: {len(small & exclude):,} < {MIN_LEN:,} bp, "
          f"{int(calls['excluded'].sum())} tested")

    data = pd.DataFrame({
        "Sample": np.tile(samples, nC),
        "Chromosome": np.repeat(scaffolds, nS),
        "Depth": depth.ravel(), "Heterozygosity": het_rate.ravel(),
        "Normalized_Depth": nd.ravel(), "Normalized_Heterozygosity": nh.ravel(),
    })

    def plot(path):
        colours = {"zero het": "#984ea3", "low het in all": "#e66101",
                   HALF_MOST: "#a65628", XZ: "#377eb8",
                   "hemizygous pattern, full depth": "#4daf4a"}
        fig, ax = plt.subplots(1, 2, figsize=(11, 5))
        c = calls[calls["call"] != "not tested"]
        ax[0].scatter(c["depth_med"].clip(0, 2), c["het_med"].clip(0, 3), s=6, c="grey",
                      alpha=0.4, label="kept")
        for k, col in colours.items():
            s = c[c["call"] == k]
            if len(s):
                ax[0].scatter(s["depth_med"].clip(0, 2), s["het_med"].clip(0, 3), s=14, c=col,
                              label=k)
        ax[0].set_xlabel("median normalised depth (per scaffold)")
        ax[0].set_ylabel("median normalised het")
        ax[0].legend(fontsize=7, frameon=False)
        d = data.assign(call=np.repeat(calls["call"].to_numpy(), nS))
        rest = d[~d["call"].isin(list(colours))]
        if len(rest) > 50_000:
            rest = rest.sample(50_000, random_state=1)
        ax[1].scatter(rest["Normalized_Depth"].clip(0, 2),
                      rest["Normalized_Heterozygosity"].clip(0, 3), s=2, c="grey", alpha=0.2)
        for k, col in colours.items():
            s = d[d["call"] == k]
            if len(s):
                ax[1].scatter(s["Normalized_Depth"].clip(0, 2),
                              s["Normalized_Heterozygosity"].clip(0, 3), s=3, c=col, alpha=0.5)
        ax[1].set_xlabel("normalised depth (sample x scaffold)")
        ax[1].set_ylabel("normalised het")
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)

    write_outputs(out, exclude, calls, data, plot)


if __name__ == "__main__":
    p = snakemake.params
    main(snakemake.input["counts"], snakemake.input["samples"], snakemake.input["contigs"],
         dict(report=snakemake.output["report"], data=snakemake.output["data"],
              plot=snakemake.output["plot"], exclude=snakemake.output["exclude_scaffolds"],
              calls=snakemake.output["calls"]),
         int(p["window"]), int(p["step"]), int(p["rule_min_size"]),
         list(p["extra_exclude"] or []), list(p["keep"] or []))
