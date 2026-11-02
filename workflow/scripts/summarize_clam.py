"""Aggregate clam per-window TSVs into two genome-wide CSV summaries.

Everything is restricted to scaffolds >= min_scaffold_size, so that pi, dxy
and Fst are all measured over the same sequence.

Per-sample summary (sums across windows; clam loci runs --per-sample, so
heterozygosity.tsv has one row per sample per window):
  heterozygosity          = Σhet_total      / Σcallable_total  [het.tsv]
  heterozygosity_not_roh  = Σhet_not_in_roh / Σcallable_not_in_roh

Per-pop summary (sums across windows):
  pi                      = Σdifferences / Σcomparisons        [pi.tsv]
  heterozygosity          = pooled: Σhet over all samples / Σcallable
  heterozygosity_not_roh  = pooled, as above, outside ROH
  heterozygosity_mean_ind = mean of the samples' heterozygosity, each sample
  heterozygosity_not_roh_mean_ind   weighted equally

Pairwise summary (only when >=2 pops; header-only otherwise):
  dxy          = Σdifferences / Σcomparisons       [dxy.tsv]
  fst_hudson   = (dxy - (pi1+pi2)/2) / dxy         Hudson's Fst as a ratio of
                 averages (Bhatia et al. 2013): the genome-wide estimate
  fst_win_mean = mean of per-window fst            [fst.tsv]; descriptive
                 only, it weights every window equally and is biased
"""

# Send stdout and stderr to the rule's log file. Snakemake does this for R
# scripts (via sink) but not for python ones, which otherwise leaves the log
# empty and the traceback buried in the master log.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path

import pandas as pd


PER_SAMPLE_COLS = [
    "sample", "population",
    "heterozygosity", "heterozygosity_not_roh",
    "het_sites", "callable_sites", "het_sites_not_roh", "callable_sites_not_roh",
]
PER_POP_COLS = [
    "population", "n_samples",
    "pi", "heterozygosity", "heterozygosity_not_roh",
    "heterozygosity_mean_ind", "heterozygosity_not_roh_mean_ind",
    "callable_sites", "callable_sites_not_roh",
]
PAIRWISE_COLS = [
    "population1", "population2",
    "dxy", "fst_hudson", "fst_win_mean",
    "n_windows",
]


pi_path = Path(snakemake.input.pi)
het_path = Path(snakemake.input.het)
dxy_path = Path(snakemake.params.dxy)
fst_path = Path(snakemake.params.fst)
n_pops = int(Path(snakemake.input.n_pops).read_text().strip())
min_scaffold_size = int(snakemake.params.min_scaffold_size)

pop_n = (
    pd.read_csv(snakemake.input.samples_tsv, sep="\t", usecols=["population"])
    .groupby("population", as_index=False).size().rename(columns={"size": "n_samples"})
)

# Scaffolds >= min_scaffold_size; every statistic below is summed over these
fai = pd.read_csv(
    snakemake.input.fai, sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"]
)
big_scaffolds = set(fai.loc[fai["length"] >= min_scaffold_size, "chrom"])

# --- per-population pi ----------------------------------------------------
pi_df = pd.read_csv(pi_path, sep="\t")
pi_df = pi_df[pi_df["chrom"].isin(big_scaffolds)]
pi_agg = pi_df.groupby("population", as_index=False).agg(
    differences=("differences", "sum"),
    comparisons=("comparisons", "sum"),
)
pi_agg["pi"] = pi_agg["differences"] / pi_agg["comparisons"].replace(0, pd.NA)

# --- per-population het + het_not_roh -------------------------------------
het_df = pd.read_csv(het_path, sep="\t")
het_df = het_df[het_df["chrom"].isin(big_scaffolds)]
het_agg = het_df.groupby("population", as_index=False).agg(
    het_total=("het_total", "sum"),
    callable_sites=("callable_total", "sum"),
    het_total_not_roh=("het_not_in_roh", "sum"),
    callable_sites_not_roh=("callable_not_in_roh", "sum"),
)
het_agg["heterozygosity"] = het_agg["het_total"] / het_agg["callable_sites"].replace(0, pd.NA)
het_agg["heterozygosity_not_roh"] = (
    het_agg["het_total_not_roh"] / het_agg["callable_sites_not_roh"].replace(0, pd.NA)
)

# --- per-sample het + het_not_roh -----------------------------------------
per_sample = het_df.groupby(["sample", "population"], as_index=False).agg(
    het_sites=("het_total", "sum"),
    callable_sites=("callable_total", "sum"),
    het_sites_not_roh=("het_not_in_roh", "sum"),
    callable_sites_not_roh=("callable_not_in_roh", "sum"),
)
per_sample["heterozygosity"] = (
    per_sample["het_sites"] / per_sample["callable_sites"].replace(0, pd.NA)
)
per_sample["heterozygosity_not_roh"] = (
    per_sample["het_sites_not_roh"] / per_sample["callable_sites_not_roh"].replace(0, pd.NA)
)
per_sample = per_sample[PER_SAMPLE_COLS].sort_values(["population", "sample"])
per_sample.to_csv(snakemake.output.per_sample, index=False, float_format="%.6g")

mean_ind = per_sample.groupby("population", as_index=False).agg(
    heterozygosity_mean_ind=("heterozygosity", "mean"),
    heterozygosity_not_roh_mean_ind=("heterozygosity_not_roh", "mean"),
)

per_pop = (
    pop_n
    .merge(pi_agg[["population", "pi"]], on="population", how="left")
    .merge(
        het_agg[["population", "heterozygosity", "heterozygosity_not_roh",
                 "callable_sites", "callable_sites_not_roh"]],
        on="population", how="left",
    )
    .merge(mean_ind, on="population", how="left")
)
per_pop = per_pop[PER_POP_COLS].sort_values("population").reset_index(drop=True)
per_pop.to_csv(snakemake.output.per_pop, index=False, float_format="%.6g")

# --- pairwise dxy + fst ---------------------------------------------------
if n_pops >= 2 and dxy_path.exists():
    dxy_df = pd.read_csv(dxy_path, sep="\t")
    dxy_df = dxy_df[dxy_df["chrom"].isin(big_scaffolds)]
    pw = dxy_df.groupby(["population1", "population2"], as_index=False).agg(
        differences=("differences", "sum"),
        comparisons=("comparisons", "sum"),
    )
    pw["dxy"] = pw["differences"] / pw["comparisons"].replace(0, pd.NA)

    pi_lookup = pi_agg.set_index("population")["pi"]

    def hudson(row):
        pi1 = pi_lookup.get(row["population1"], float("nan"))
        pi2 = pi_lookup.get(row["population2"], float("nan"))
        d = row["dxy"]
        if pd.isna(d) or d == 0:
            return float("nan")
        return (d - (pi1 + pi2) / 2.0) / d

    pw["fst_hudson"] = pw.apply(hudson, axis=1)

    if fst_path.exists():
        fst_df = pd.read_csv(fst_path, sep="\t")
        fst_df = fst_df[fst_df["chrom"].isin(big_scaffolds)]
        fst_agg = fst_df.groupby(["population1", "population2"], as_index=False).agg(
            n_windows=("fst", "count"),
            fst_win_mean=("fst", "mean"),
        )
        pw = pw.merge(fst_agg, on=["population1", "population2"], how="left")
    else:
        pw["n_windows"] = 0
        pw["fst_win_mean"] = float("nan")

    pairwise = pw[PAIRWISE_COLS].sort_values(["population1", "population2"]).reset_index(drop=True)
else:
    pairwise = pd.DataFrame(columns=PAIRWISE_COLS)

pairwise.to_csv(snakemake.output.pairwise, index=False, float_format="%.6g")
