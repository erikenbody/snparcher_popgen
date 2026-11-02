"""Pick, for each K, the ADMIXTURE replicate with the highest log-likelihood.

Replicates differ only in their random starting point, so the one that reached
the highest likelihood is the best fit at that K (as CLUMPAK and pong do). Its
Q matrix and log are copied to fixed per-K names, so plots and population
assignment never depend on which seed happened to be listed first.
"""

# Send stdout and stderr to the rule's log file. Snakemake does this for R
# scripts (via sink) but not for python ones, which otherwise leaves the log
# empty and the traceback buried in the master log.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

import re
import shutil
from pathlib import Path


def final_loglik(log):
    """The converged log-likelihood ADMIXTURE prints after 'Summary:'."""
    for line in Path(log).read_text().splitlines():
        m = re.match(r"Loglikelihood:\s*(\S+)", line)
        if m:
            return float(m.group(1))
    raise ValueError(f"no final Loglikelihood line in {log}")


def main(logs, qfiles, ks, seeds, out_q, out_log, table):
    log_of = {}
    q_of = {}
    for log in logs:
        m = re.search(r"_K(\d+)_seed(\d+)\.log$", log)
        log_of[int(m.group(1)), int(m.group(2))] = log
    for q in qfiles:
        m = re.search(r"\.(\d+)\.(\d+)\.Q$", q)
        q_of[int(m.group(1)), int(m.group(2))] = q

    rows = []
    for k, q_dest, log_dest in zip(ks, out_q, out_log):
        lik = {s: final_loglik(log_of[k, s]) for s in seeds}
        best = max(seeds, key=lambda s: lik[s])
        shutil.copyfile(q_of[k, best], q_dest)
        shutil.copyfile(log_of[k, best], log_dest)
        spread = max(lik.values()) - min(lik.values())
        print(f"K={k}: seed {best} (loglik {lik[best]:.3f}; range across seeds {spread:.3f})")
        for s in seeds:
            rows.append(f"{k}\t{s}\t{lik[s]:.6f}\t{str(s == best).lower()}")

    with open(table, "w") as fh:
        fh.write("K\tseed\tloglikelihood\tbest\n")
        fh.write("\n".join(rows) + "\n")


main(
    logs=list(snakemake.input.logs),
    qfiles=list(snakemake.input.qfiles),
    ks=[int(k) for k in snakemake.params.ks],
    seeds=[int(s) for s in snakemake.params.seeds],
    out_q=list(snakemake.output.q),
    out_log=list(snakemake.output.logs),
    table=snakemake.output.table,
)
