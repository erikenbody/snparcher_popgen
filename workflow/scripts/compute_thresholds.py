"""Compute per-chromosome Poisson-SD depth thresholds from a `clam collect`
depths zarr.

For each contig, the depth array has shape (chrom_length, n_samples). We take
its overall mean (averaged across samples and positions) and build Poisson
bounds:
    min_depth = max(0, mean - sd_scale * sqrt(mean))
    max_depth =        mean + sd_scale * sqrt(mean)

Output is the tab-separated file expected by `clam loci --thresholds-file`:
    contig\tmin_depth\tmax_depth

This matches snpArcher's per-chromosome cov_filter formula
(create_coverage_bed.py: stdev = sqrt(mean), bounds = mean +/- N*stdev),
differing only in that clam applies the per-chrom bounds per-individual
rather than to the per-position cohort mean.

Pattern lifted from clam/notebooks/working_with_collect_zarr.ipynb.
"""

# Send stdout and stderr to the rule's log file. Snakemake does this for R
# scripts (via sink) but not for python ones, which otherwise leaves the log
# empty and the traceback buried in the master log.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

import math
from pathlib import Path
from typing import Any, cast

import dask.array as da
import numpy as np
import pandas as pd
import zarr


def _read_excluded(path):
    if not path:
        return set()
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return set()
    return {line.strip() for line in p.read_text().splitlines() if line.strip()}


def main(snakemake):
    zarr_path = Path(snakemake.input.depths_zarr)
    out_path = Path(snakemake.output.thresholds)
    sd_scale = float(snakemake.params.sd_scale)
    excluded = _read_excluded(getattr(snakemake.input, "exclude_file", None))

    store = zarr.open(str(zarr_path), mode="r")
    metadata = cast(dict[str, Any], store.attrs["clam_metadata"])
    contigs = [c["name"] for c in metadata["contigs"] if c["name"] not in excluded]

    # Per-contig mean computed in parallel (threaded scheduler).
    means_lazy = {c: da.from_zarr(str(zarr_path / c)).mean() for c in contigs}
    means = da.compute(means_lazy)[0]

    rows = []
    for contig in contigs:
        m = float(means[contig])
        if not np.isfinite(m) or m <= 0:
            # Empty / no-coverage contig — leave the bounds wide so clam keeps it
            # callable wherever individual samples have any depth at all.
            rows.append({"contig": contig, "min_depth": 0.0, "max_depth": float("inf")})
            continue
        sd = math.sqrt(m)
        rows.append({
            "contig": contig,
            "min_depth": max(0.0, m - sd_scale * sd),
            "max_depth": m + sd_scale * sd,
        })

    df = pd.DataFrame(rows, columns=["contig", "min_depth", "max_depth"])
    df.to_csv(out_path, sep="\t", index=False, float_format="%.6f")


main(snakemake)  # noqa: F821
