"""Per-sample sequencing depth from a `clam collect` depths zarr.

The depths zarr stores, per contig, an array of shape (contig_length, n_samples)
of per-base depth (GVCF DP). Sample names are in
metadata["clam_metadata"]["column_names"] (axis-1 order). The contigs present
are already post-exclusion (clam collect was run with --exclude-file), so we
simply iterate every contig in the metadata.

Depth is averaged over ALL reference positions (including zero-depth /
non-callable bases) -> standard mean genome coverage per individual.

Outputs:
  - sample_depth.csv          : one row per sample
        sample, mean_depth, std_depth, n_positions
  - sample_depth_per_contig.csv : wide matrix, one row per contig
        contig, length, <sample_1>, <sample_2>, ...   (mean depth)

Single streaming pass over the zarr via dask (same read pattern as
compute_thresholds.py), reducing over the position axis while keeping the
sample axis. Per-contig mean and mean-of-squares are accumulated lazily and
computed in one da.compute call; the genome-wide mean/std are then a
length-weighted combination, so no second pass is needed.
"""

# Send stdout and stderr to the rule's log file. Snakemake does this for R
# scripts (via sink) but not for python ones, which otherwise leaves the log
# empty and the traceback buried in the master log.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path
from typing import Any, cast

import dask.array as da
import numpy as np
import pandas as pd
import zarr


def main(snakemake):
    zarr_path = Path(snakemake.input.depths_zarr)
    out_summary = Path(snakemake.output.summary)
    out_per_contig = Path(snakemake.output.per_contig)

    store = zarr.open(str(zarr_path), mode="r")
    metadata = cast(dict[str, Any], store.attrs["clam_metadata"])
    samples = list(metadata["column_names"])
    contigs = [(c["name"], int(c["length"])) for c in metadata["contigs"]]

    # Lazily build per-contig per-sample mean and mean-of-squares (reduce over
    # the position axis, keep the sample axis). cast to float64 so the squared
    # sum can't overflow the stored integer depth dtype.
    lazy = {}
    for name, _length in contigs:
        arr = da.from_zarr(str(zarr_path / name)).astype("float64")
        lazy[name] = {"mean": arr.mean(axis=0), "meansq": (arr ** 2).mean(axis=0)}
    res = da.compute(lazy)[0]

    # Per-contig matrix (mean depth, samples as columns).
    per_contig_rows = []
    for name, length in contigs:
        row = {"contig": name, "length": length}
        row.update({s: float(m) for s, m in zip(samples, res[name]["mean"])})
        per_contig_rows.append(row)
    per_contig_df = pd.DataFrame(per_contig_rows, columns=["contig", "length"] + samples)
    per_contig_df.to_csv(out_per_contig, index=False, float_format="%.4f")

    # Genome-wide per-sample: length-weighted combination of per-contig stats.
    #   mean   = sum_c length_c * mean_c   / sum_c length_c
    #   meansq = sum_c length_c * meansq_c / sum_c length_c
    #   std    = sqrt(max(0, meansq - mean^2))
    nsamp = len(samples)
    total_len = sum(length for _n, length in contigs)
    wsum_mean = np.zeros(nsamp, dtype="float64")
    wsum_meansq = np.zeros(nsamp, dtype="float64")
    for name, length in contigs:
        wsum_mean += length * np.asarray(res[name]["mean"], dtype="float64")
        wsum_meansq += length * np.asarray(res[name]["meansq"], dtype="float64")
    mean = wsum_mean / total_len
    meansq = wsum_meansq / total_len
    std = np.sqrt(np.clip(meansq - mean ** 2, 0, None))

    summary_df = pd.DataFrame({
        "sample": samples,
        "mean_depth": mean,
        "std_depth": std,
        "n_positions": total_len,
    })
    summary_df.to_csv(out_summary, index=False, float_format="%.4f")

    print(f"[sample_depth] {nsamp} samples over {len(contigs)} contigs "
          f"({total_len:,} bp) -> {out_summary}")


main(snakemake)  # noqa: F821
