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

def calc_roh(roh, fai, output, population, min_roh_size, min_scaffold_size):
    dffai = pd.read_table(fai, sep='\t', header=None)
    # Restrict to long scaffolds: F_ROH is a fraction of the genome, and short
    # scaffolds contribute length to the denominator while being too short to
    # detect a run of homozygosity in.
    dffai = dffai[dffai[1] > int(min_scaffold_size)]
    chroms = dffai[0].values

    if dffai.empty:
        raise SystemExit(
            f"no scaffolds longer than min_scaffold_size ({min_scaffold_size} bp) "
            f"in {fai}; lower min_scaffold_size for this assembly"
        )

    glength = dffai[1].sum()
    
    dfroh = pd.read_table(roh, sep='\t', header=None, names=['chrom', 'start', 'end', 'sample'])
    dfroh['length'] = dfroh['end'] - dfroh['start']
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


def main():
    roh_input = snakemake.input["rg"]
    fai_input = snakemake.input["fai"]
    froh_output = snakemake.output["froh"]
    population = snakemake.params["population"]
    min_roh_size = snakemake.params["min_roh_size"]
    min_scaffold_size = snakemake.params["min_scaffold_size"]

    calc_roh(roh_input, fai_input, froh_output, population, min_roh_size,
             min_scaffold_size)

if __name__ == "__main__":
    main()
