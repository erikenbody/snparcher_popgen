"""Collect GONE2 and currentNe2 runs into tables and decide, per population,
whether the panmictic or the metapopulation (-x) model applies.

Run directories are ne/gone2/population_{k}_{mode}/seed_{seed}/ (GONE2 is run
once per seed, as the GONE2 README recommends for gauging run-to-run
variation) and ne/currentne2/population_{k}_{mode}/. Each holds an
exit_status file. Failed runs are kept (status table) rather than dropped, so
a project with no usable Ne estimate still finishes; filter downstream.

Model choice. The -x runs fit the sample's excess of homozygotes (Fit) as
differentiation among equal-sized subpopulations (Fst) plus deviation from
Hardy-Weinberg within them (Fis). A population is called structured when
that fit puts Fst at or above ne_structure_min_fst; otherwise the -x model
has found no subdivision to correct for and the standard run is used. The
Fst comes from currentNe2 -x, or from GONE2 -x (median over seeds) when
currentNe2 -x has no output.

Outputs:
  gone2          project, population, mode, seed, generation, ne
  gone2_summary  project, population, mode, generation, n_seeds, mean, sd,
                 median, min, max (across seeds)
  gone2_change   project, population, mode, generation, ne (mean across
                 seeds), ne_change (to the next listed generation back)
  currentne2     project, population, mode, section, parameter, value
                 (every "# label:" / value pair of the output file; section is
                 the "# Ne estimation ..." heading a value falls under, if any)
  estimates      project, population, mode, method, estimate, lo50, hi50,
                 lo90, hi90, interval (one row per Ne point estimate)
  model          project, population, n_samples, chosen_mode, structured,
                 fst, fis, fit, migration, n_subpops, fst_source, gone2_fst,
                 min_fst, reason
  status         project, tool, population, mode, seed, exit_status, has_output
"""

# Send stdout and stderr to the rule's log file.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

import re
from pathlib import Path

import numpy as np
import pandas as pd

GONE2_RUN = re.compile(r"population_(\d+)_(standard|metapop)/seed_(\d+)$")
CURRENTNE2_RUN = re.compile(r"population_(\d+)_(standard|metapop)$")
# Generations reported in the change table, as in the original CCGP GONE
# summary (GONE2 reports 150 generations, so 200 is gone).
CHANGE_GENERATIONS = [1, 10, 50, 100, 150]


def exit_status(rundir):
    path = rundir / "exit_status"
    return path.read_text().strip() if path.exists() else "NA"


def label_values(path):
    """Every "# label:" line followed by a value line, as (section, label,
    value). section is the latest "# Ne estimation ..." heading, if any."""
    rows, section, label = [], "", None
    for line in path.read_text().splitlines():
        if line.startswith("# Ne estimation"):
            section = line.lstrip("# ").rstrip(".:")
            label = None
        elif line.startswith("#"):
            text = line.lstrip("# ").strip()
            label = text[:-1] if text.endswith(":") else None
        elif label is not None and line.strip():
            rows.append((section, label, line.strip()))
            label = None
    return rows


def read_gone2(rundir, population):
    """The _GONE2_Ne table (_GONE2_mix_Ne under -x) and the _STATS values,
    or (None, {}) if absent."""
    for stem in ("GONE2", "GONE2_mix"):
        ne_path = rundir / f"population_{population}_{stem}_Ne"
        if ne_path.exists():
            rows = []
            for line in ne_path.read_text().splitlines():
                fields = line.split()
                try:
                    rows.append((int(float(fields[0])), float(fields[1])))
                except (ValueError, IndexError):
                    continue  # header
            stats_path = rundir / f"population_{population}_{stem}_STATS"
            stats = (
                {label: value for _, label, value in label_values(stats_path)}
                if stats_path.exists() else {}
            )
            return pd.DataFrame(rows, columns=["generation", "ne"]), stats
    return None, {}


def read_currentne2(rundir):
    paths = list(rundir.glob("*_OUTPUT.txt"))
    if not paths:
        return None
    return pd.DataFrame(label_values(paths[0]), columns=["section", "parameter", "value"])


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def stat(stats, prefix):
    """The value of the first label starting with prefix."""
    for label, value in stats.items():
        if label.startswith(prefix):
            return number(value)
    return np.nan


def currentne2_estimates(df):
    """Point estimates with 50%/90% CIs from one currentNe2 output. Labels
    vary ("Lower limit of 50% CI", "Upper bound of 50% CI", "Lower 50% limit
    of the N_T estimate"), so the bounds are matched loosely."""
    out = []
    for section, block in df.groupby("section", sort=False):
        values = dict(zip(block["parameter"], block["value"]))
        point = next((v for k, v in values.items()
                      if k == "Ne point estimate" or k.startswith("N_T of the metapopulation")), None)
        if point is None:
            continue
        bounds = {}
        for label, value in values.items():
            m = re.match(r"(Lower|Upper)\b.*?\b(50|90)%", label)
            if m:
                bounds[("lo" if m.group(1) == "Lower" else "hi") + m.group(2)] = number(value)
        if "between chromosomes" in section:
            method = "currentNe2, between chromosomes"
        elif section:
            method = "currentNe2"
        else:
            method = "currentNe2 -x (N_T)"
        out.append(dict(method=method, estimate=number(point), interval="CI", **bounds))
    return out


def choose_model(population, cn2, gone2_stats, modes, min_fst):
    """Structured or panmictic, from the -x fits. Returns one model row."""
    cn2_x = cn2.get((population, "metapop"))
    values = dict(zip(cn2_x["parameter"], cn2_x["value"])) if cn2_x is not None else {}
    fst = stat(values, "Fst")
    fis = stat(values, "Fis")
    fit = stat(values, "Fit")
    migration = stat(values, "Migration rate")
    n_subpops = stat(values, "Number of subpopulations")
    gone2_x = [s for (p, m), runs in gone2_stats.items() if p == population and m == "metapop"
               for s in runs]
    gone2_fst = np.nanmedian([stat(s, "Fst") for s in gone2_x]) if gone2_x else np.nan
    source = "currentNe2 -x"
    if np.isnan(fst) and not np.isnan(gone2_fst):
        fst, source = gone2_fst, "GONE2 -x (median over seeds)"
        migration = np.nanmedian([stat(s, "Migration rate") for s in gone2_x])
        n_subpops = np.nanmedian([stat(s, "Number of subpopulations") for s in gone2_x])
    if np.isnan(fst):
        source = "none"

    if "metapop" not in modes:
        structured, chosen, reason = np.nan, "standard", "only standard runs configured"
    elif "standard" not in modes:
        structured, chosen, reason = np.nan, "metapop", "only metapop runs configured"
    elif np.isnan(fst):
        structured, chosen, reason = np.nan, "standard", "no -x fit available; defaulting to standard"
    elif fst >= min_fst:
        structured, chosen = True, "metapop"
        reason = f"-x Fst {fst:.3f} >= {min_fst}: sample spans differentiated subpopulations"
    else:
        structured, chosen = False, "standard"
        reason = f"-x Fst {fst:.3f} < {min_fst}: no subdivision for -x to correct"
    return dict(chosen_mode=chosen, structured=structured, fst=fst, fis=fis, fit=fit,
                migration=migration, n_subpops=n_subpops, fst_source=source,
                gone2_fst=gone2_fst, min_fst=min_fst, reason=reason)


def main():
    project = snakemake.params["project_id"]
    min_fst = float(snakemake.params["min_fst"])
    modes = {"gone2": list(snakemake.params["gone2_modes"]),
             "currentne2": list(snakemake.params["currentne2_modes"])}
    all_modes = set(modes["gone2"]) | set(modes["currentne2"])

    status_rows, gone2_frames, gone2_stats, n_samples = [], [], {}, {}
    for rundir in map(Path, snakemake.input.get("gone2", [])):
        population, mode, seed = GONE2_RUN.search(str(rundir)).groups()
        population, seed = int(population), int(seed)
        df, stats = read_gone2(rundir, population)
        has_output = df is not None and not df.empty
        status_rows.append((project, "gone2", population, mode, seed, exit_status(rundir), has_output))
        if has_output:
            df.insert(0, "seed", seed)
            df.insert(0, "mode", mode)
            df.insert(0, "population", population)
            df.insert(0, "project", project)
            gone2_frames.append(df)
            gone2_stats.setdefault((population, mode), []).append(stats)
            n = stat(stats, "Total number of individuals")
            if not np.isnan(n):
                n_samples[population] = int(n)

    cn2, cn2_frames, estimate_rows = {}, [], []
    for rundir in map(Path, snakemake.input.get("currentne2", [])):
        population, mode = CURRENTNE2_RUN.search(rundir.name).groups()
        population = int(population)
        df = read_currentne2(rundir)
        has_output = df is not None and not df.empty
        status_rows.append((project, "currentne2", population, mode, "", exit_status(rundir), has_output))
        if has_output:
            cn2[(population, mode)] = df
            for est in currentne2_estimates(df):
                estimate_rows.append(dict(project=project, population=population, mode=mode, **est))
            n = stat(dict(zip(df["parameter"], df["value"])), "Total number of individuals")
            if not np.isnan(n):
                n_samples.setdefault(population, int(n))
            df = df.copy()
            df.insert(0, "mode", mode)
            df.insert(0, "population", population)
            df.insert(0, "project", project)
            cn2_frames.append(df)

    gone2 = (pd.concat(gone2_frames, ignore_index=True) if gone2_frames else
             pd.DataFrame(columns=["project", "population", "mode", "seed", "generation", "ne"]))
    gone2.to_csv(snakemake.output["gone2"], sep="\t", index=False)

    summary = (
        gone2.groupby(["project", "population", "mode", "generation"])["ne"]
        .agg(n_seeds="count", mean="mean", sd="std", median="median", min="min", max="max")
        .reset_index()
    )
    summary.to_csv(snakemake.output["gone2_summary"], sep="\t", index=False, float_format="%.4g")

    change = summary[summary["generation"].isin(CHANGE_GENERATIONS)][
        ["project", "population", "mode", "generation", "mean"]].rename(columns={"mean": "ne"})
    change = change.sort_values(["population", "mode", "generation"])
    change["ne_change"] = change.groupby(["population", "mode"])["ne"].diff(-1)
    change.to_csv(snakemake.output["gone2_change"], sep="\t", index=False, float_format="%.1f")

    cn2_cols = ["project", "population", "mode", "section", "parameter", "value"]
    (pd.concat(cn2_frames, ignore_index=True) if cn2_frames else pd.DataFrame(columns=cn2_cols)
     ).to_csv(snakemake.output["currentne2"], sep="\t", index=False)

    # GONE2's most recent generation, as a point estimate comparable to
    # currentNe2's contemporary Ne; the interval is the range across seeds.
    for (population, mode), df in summary[summary["generation"] == 1].groupby(["population", "mode"]):
        row = df.iloc[0]
        estimate_rows.append(dict(
            project=project, population=population, mode=mode,
            method="GONE2, generation 1" + (" (N_T)" if mode == "metapop" else ""),
            estimate=row["median"], lo90=row["min"], hi90=row["max"], interval="seed range",
        ))
    estimates = pd.DataFrame(estimate_rows, columns=[
        "project", "population", "mode", "method", "estimate", "lo50", "hi50", "lo90", "hi90", "interval"])
    estimates.sort_values(["population", "mode", "method"]).to_csv(
        snakemake.output["estimates"], sep="\t", index=False, float_format="%.4g")

    populations = sorted({r[2] for r in status_rows})
    model = pd.DataFrame([
        dict(project=project, population=p, n_samples=n_samples.get(p, np.nan),
             **choose_model(p, cn2, gone2_stats, all_modes, min_fst))
        for p in populations
    ], columns=["project", "population", "n_samples", "chosen_mode", "structured", "fst", "fis",
                "fit", "migration", "n_subpops", "fst_source", "gone2_fst", "min_fst", "reason"])
    model.to_csv(snakemake.output["model"], sep="\t", index=False, float_format="%.4g")

    status = pd.DataFrame(status_rows, columns=[
        "project", "tool", "population", "mode", "seed", "exit_status", "has_output"])
    status.to_csv(snakemake.output["status"], sep="\t", index=False)
    print(status.to_string(index=False))
    print(model.to_string(index=False))


main()
