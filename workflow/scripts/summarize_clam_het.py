"""Genome-wide heterozygosity outside ROH from one clam heterozygosity table.

For the extra ROH cutoffs (clam_extra_roh_min_sizes). Same sums and the same
scaffold restriction (>= min_scaffold_size) as summarize_clam.py, so the
numbers line up with its per-sample and per-pop columns. One CSV:

  level       sample, or "population" for the pooled row of each population
  heterozygosity          = Σhet_total      / Σcallable_total
  heterozygosity_not_roh  = Σhet_not_in_roh / Σcallable_not_in_roh
  roh_fraction_callable   = 1 - Σcallable_not_in_roh / Σcallable_total
"""

import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

import pandas as pd


COLS = [
    "level", "sample", "population",
    "heterozygosity", "heterozygosity_not_roh", "roh_fraction_callable",
    "het_sites", "callable_sites", "het_sites_not_roh", "callable_sites_not_roh",
]

fai = pd.read_csv(
    snakemake.input.fai, sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"]
)
big = set(fai.loc[fai["length"] >= int(snakemake.params.min_scaffold_size), "chrom"])

het = pd.read_csv(snakemake.input.het, sep="\t")
het = het[het["chrom"].isin(big)]
sums = dict(
    het_sites=("het_total", "sum"),
    callable_sites=("callable_total", "sum"),
    het_sites_not_roh=("het_not_in_roh", "sum"),
    callable_sites_not_roh=("callable_not_in_roh", "sum"),
)
per_sample = het.groupby(["sample", "population"], as_index=False).agg(**sums)
per_sample["level"] = "sample"
per_pop = het.groupby("population", as_index=False).agg(**sums)
per_pop["level"] = "population"
per_pop["sample"] = ""

df = pd.concat([per_sample, per_pop], ignore_index=True)
df["heterozygosity"] = df["het_sites"] / df["callable_sites"].replace(0, pd.NA)
df["heterozygosity_not_roh"] = df["het_sites_not_roh"] / df["callable_sites_not_roh"].replace(0, pd.NA)
df["roh_fraction_callable"] = 1 - df["callable_sites_not_roh"] / df["callable_sites"].replace(0, pd.NA)
# clam writes the not-in-ROH counts as floats (NaN-able); keep counts exact.
for c in ["het_sites", "callable_sites", "het_sites_not_roh", "callable_sites_not_roh"]:
    df[c] = df[c].round().astype("Int64")
df = df[COLS].sort_values(["level", "population", "sample"], ascending=[False, True, True])
df.to_csv(snakemake.output.summary, index=False, float_format="%.6g")
