# Send stdout and stderr to the rule's log file. Snakemake does this for R
# scripts (via sink) but not for python ones, which otherwise leaves the log
# empty and the traceback buried in the master log.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import statistics
import gzip
import tempfile
import psutil
import os

def read_excluded(path):
    """Scaffold names from the sex_id exclude list (one per line; may be empty)."""
    if not path:
        return set()
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def calc_roh(roh, fai, output, population, min_roh_size, min_scaffold_size, exclude=None,
             output_frac=None, min_roh_frac=None):
    dffai = pd.read_table(fai, sep='\t', header=None, dtype={0: str})
    # Restrict to long scaffolds: F_ROH is a fraction of the genome, and short
    # scaffolds contribute length to the denominator while being too short to
    # detect a run of homozygosity in.
    dffai = dffai[dffai[1] > int(min_scaffold_size)]
    # Scaffolds excluded by sex_id are not in the population VCF, so no ROH can be
    # called on them: their length must leave the denominator too. Keeping it
    # there biased F_ROH low by the excluded share of the long scaffolds.
    excluded = read_excluded(exclude)
    n_long = len(dffai)
    dffai = dffai[~dffai[0].isin(excluded)]
    print(f"{n_long} scaffolds > {min_scaffold_size} bp; {n_long - len(dffai)} of them "
          f"excluded by sex_id; denominator {dffai[1].sum()} bp")
    chroms = dffai[0].values

    # A fragmented assembly has no scaffold long enough to measure F_ROH on.
    # Write empty outputs, as when no sample has a long ROH, so the rest of the
    # run (pi, heterozygosity, the ROH plot) goes ahead without it.
    if dffai.empty:
        print(f"No scaffolds longer than min_scaffold_size ({min_scaffold_size} bp) "
              f"left in {fai}; F_ROH not computed.")
        for path in (output, output.replace(".froh", "_top.froh"), output_frac):
            if path:
                with open(path, "w"):
                    pass
        return False

    glength = dffai[1].sum()
    
    dfroh = pd.read_table(roh, sep='\t', header=None, names=['chrom', 'start', 'end', 'sample'],
                          dtype={'chrom': str, 'sample': str})
    dfroh['length'] = dfroh['end'] - dfroh['start']
    if output_frac:
        froh_relative(dfroh, dffai, glength, output_frac, float(min_roh_frac))
    dfroh = dfroh[dfroh['length'] > int(min_roh_size)]  # only look at 500kbp ROH and longer
    dfroh = dfroh[dfroh['chrom'].isin(chroms)]

    dffroh = dfroh.groupby(
        ['sample']
    ).agg(
        {
            'length': "sum"
        }
    ).div(glength)
    dffroh = dffroh.reset_index()

    if dffroh.empty:
        print("No individuals found with ROH greater than 500kbp.")
        with open(output, 'w'):
            pass
        with open(output.replace(".froh", "_top.froh"), "w"):
            pass
        return False

    dffroh.to_csv(output, sep='\t', index=False, header=None)
    print(dffroh)

    if len(dffroh) < 8:
        combined_df = dffroh
    else:
        # get the largest 2 ROH
        df_lg = dffroh.nlargest(2, 'length')
        # get 8 random rows so we have 10 individuals total
        random_rows = dffroh.sample(n=8)
        combined_df = pd.concat([df_lg, random_rows], axis=0)
    combined_df.to_csv(output.replace(".froh", "_top.froh"), sep='\t', index=False, header=False, columns=[combined_df.columns[0]])
    #combined_df.to_csv(output, sep='\t', index=False, header=False, columns=[combined_df.columns[0]])
    return True


def froh_relative(dfroh, dffai, glength, output, min_frac):
    """F_ROH counting a merged ROH as long when it covers at least min_frac of
    its scaffold (1% by default; in humans 1 Mb is about 1% of a chromosome),
    rather than an absolute length. Same scaffolds and denominator as the
    absolute rule, and the same two columns (samples with no qualifying ROH
    are absent, as in .froh)."""
    size = dict(zip(dffai[0], dffai[1]))
    d = dfroh[dfroh['chrom'].isin(size)].copy()
    d = d[d['length'] / d['chrom'].map(size) >= min_frac]
    res = d.groupby('sample')['length'].sum().div(glength).reset_index()
    res.to_csv(output, sep='\t', index=False, header=None)
    print(f"F_ROH (ROH >= {min_frac:g} of scaffold): {len(res)} samples")


def main():
    roh_input = snakemake.input["rg"]
    fai_input = snakemake.input["fai"]
    froh_output = snakemake.output["froh"]
    population = snakemake.params["population"]
    min_roh_size = snakemake.params["min_roh_size"]
    min_scaffold_size = snakemake.params["min_scaffold_size"]
    exclude = snakemake.input.get("exclude_scaffolds")

    calc_roh(roh_input, fai_input, froh_output, population, min_roh_size,
             min_scaffold_size, exclude,
             output_frac=snakemake.output.get("froh_frac"),
             min_roh_frac=snakemake.params.get("min_roh_frac"))

if __name__ == "__main__":
    main()
