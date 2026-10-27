"""Collect GONE2 and currentNe2 runs into long tables and plot GONE2 Ne.

Each run directory is ne/{tool}/population_{k}_{mode}/ and holds an
exit_status file. Failed runs are kept (status table) rather than dropped, so
a project with no usable Ne estimate still finishes; filter downstream.

Outputs:
  gone2       project, population, mode, generation, ne
  currentne2  project, population, mode, section, parameter, value
              (every "# label:" / value pair of the output file; section is
              the "# Ne estimation ..." heading a value falls under, if any)
  status      project, tool, population, mode, exit_status, has_output
  plot        GONE2 Ne trajectories, one line per population and mode
"""

# Send stdout and stderr to the rule's log file.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

RUN_DIR = re.compile(r"population_(\d+)_(standard|metapop)$")


def run_info(rundir):
    rundir = Path(rundir)
    population, mode = RUN_DIR.search(rundir.name).groups()
    status_file = rundir / "exit_status"
    status = status_file.read_text().strip() if status_file.exists() else "NA"
    return rundir, int(population), mode, status


def read_gone2(rundir, population):
    """The _GONE2_Ne table (_GONE2_mix_Ne under -x), or None if absent."""
    for name in (f"population_{population}_GONE2_Ne", f"population_{population}_GONE2_mix_Ne"):
        path = rundir / name
        if path.exists():
            rows = []
            for line in path.read_text().splitlines():
                fields = line.split()
                try:
                    rows.append((int(float(fields[0])), float(fields[1])))
                except (ValueError, IndexError):
                    continue  # header
            return pd.DataFrame(rows, columns=["generation", "ne"])
    return None


def read_currentne2(rundir):
    """Every "# label:" line followed by a value line, or None if absent."""
    paths = list(rundir.glob("*_currentNe2_OUTPUT.txt"))
    if not paths:
        return None
    rows, section, label = [], "", None
    for line in paths[0].read_text().splitlines():
        if line.startswith("# Ne estimation"):
            section = line.lstrip("# ").rstrip(".")
            label = None
        elif line.startswith("#"):
            text = line.lstrip("# ").strip()
            label = text[:-1] if text.endswith(":") else None
        elif label is not None and line.strip():
            rows.append((section, label, line.strip()))
            label = None
    return pd.DataFrame(rows, columns=["section", "parameter", "value"])


def main():
    project = snakemake.params["project_id"]
    readers = {"gone2": read_gone2, "currentne2": lambda d, p: read_currentne2(d)}

    status_rows = []
    tables = {"gone2": [], "currentne2": []}
    for tool in tables:
        for rundir in snakemake.input.get(tool, []):
            rundir, population, mode, status = run_info(rundir)
            df = readers[tool](rundir, population)
            has_output = df is not None and not df.empty
            status_rows.append((project, tool, population, mode, status, has_output))
            if has_output:
                df.insert(0, "mode", mode)
                df.insert(0, "population", population)
                df.insert(0, "project", project)
                tables[tool].append(df)

    columns = {
        "gone2": ["project", "population", "mode", "generation", "ne"],
        "currentne2": ["project", "population", "mode", "section", "parameter", "value"],
    }
    for tool, frames in tables.items():
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns[tool])
        df.to_csv(snakemake.output[tool], sep="\t", index=False)
        tables[tool] = df

    status = pd.DataFrame(
        status_rows,
        columns=["project", "tool", "population", "mode", "exit_status", "has_output"],
    )
    status.to_csv(snakemake.output["status"], sep="\t", index=False)
    print(status.to_string(index=False))

    gone2 = tables["gone2"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    if gone2.empty:
        ax.text(0.5, 0.5, "no successful GONE2 runs", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_axis_off()
    else:
        for (population, mode), df in gone2.groupby(["population", "mode"]):
            ax.plot(df["generation"], df["ne"], label=f"population {population} ({mode})")
        ax.set_yscale("log")
        ax.set_xlabel("generations ago")
        ax.set_ylabel("Ne")
        ax.legend(frameon=False)
    ax.set_title(f"{project}: GONE2 Ne ({snakemake.params['recomb_rate']} cM/Mb)")
    fig.tight_layout()
    fig.savefig(snakemake.output["plot"], dpi=150)


main()
