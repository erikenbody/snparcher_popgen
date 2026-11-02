"""Plot the Ne tables written by summarize_ne.py.

Outputs:
  gone2_all    GONE2 trajectories for every population and mode (one panel
               per population; thin lines are seeds, thick lines their mean)
  gone2_final  publication figure: each population's trajectory under its
               chosen model over generations ne_plot_generations, mean across
               seeds with a +/- 1 SD band (as the original CCGP GONE plot),
               coloured as in the ancestry figures
  estimates    contemporary Ne from every method, one panel per population;
               the chosen model's estimates are coloured, the other greyed
"""

# Send stdout and stderr to the rule's log file.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter

# The ancestry palette (ancestry_assignment.R: chaparral1 reordered), so a
# population has the same colour here as on the ancestry map.
CHAPARRAL = ["#DCC27A", "#B0B9BE", "#63605F", "#985E5C", "#AEBFA8", "#F19B34"]
POP_COLORS = [CHAPARRAL[i] for i in (0, 3, 4, 1, 2, 5)]
MODE_COLORS = {"standard": "#2F5F98", "metapop": "#C8642A"}
MODE_LABELS = {"standard": "standard (panmictic)", "metapop": "-x (metapopulation, N_T)"}
GREY = "#B5B5B5"
INK = "#222222"
# Darkens the edges of the palette's pale colours so lines read on white.
OUTLINE = [pe.Stroke(linewidth=3.2, foreground="#4A4A4A"), pe.Normal()]
METHODS = [
    ("standard", "currentNe2"),
    ("standard", "currentNe2, between chromosomes"),
    ("standard", "GONE2, generation 1"),
    ("metapop", "currentNe2 -x (N_T)"),
    ("metapop", "GONE2, generation 1 (N_T)"),
]

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.edgecolor": "#555555",
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})


def pop_color(population):
    return POP_COLORS[(int(population) - 1) % len(POP_COLORS)]


def pop_label(population, model):
    row = model[model["population"] == population]
    n = row["n_samples"].iloc[0] if len(row) else np.nan
    return f"Population {population}" + ("" if pd.isna(n) else f" (n = {int(n)})")


def chosen_mode(population, model):
    row = model[model["population"] == population]
    return row["chosen_mode"].iloc[0] if len(row) else "standard"


def log_axis(ax, axis="y"):
    getattr(ax, f"set_{axis}scale")("log")
    target = ax.yaxis if axis == "y" else ax.xaxis
    target.set_minor_locator(LogLocator(subs=np.arange(2, 10)))
    target.set_minor_formatter(NullFormatter())


def empty(ax, text):
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, color="#666666")
    ax.set_axis_off()


def save(fig, paths):
    for path in paths:
        fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_gone2_all(gone2, summary, model, project, rate, paths):
    populations = sorted(set(model["population"]))
    ncol = min(3, max(1, len(populations)))
    nrow = max(1, math.ceil(len(populations) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 2.8 * nrow), squeeze=False,
                             sharex=True)
    for ax in axes.flat[len(populations):]:
        ax.set_axis_off()
    for ax, population in zip(axes.flat, populations):
        pop_runs = gone2[gone2["population"] == population]
        if pop_runs.empty:
            empty(ax, f"{pop_label(population, model)}\nno successful GONE2 runs")
            continue
        for mode, color in MODE_COLORS.items():
            for _, seed_df in pop_runs[pop_runs["mode"] == mode].groupby("seed"):
                ax.plot(seed_df["generation"], seed_df["ne"], color=color, lw=0.6, alpha=0.35)
            mean = summary[(summary["population"] == population) & (summary["mode"] == mode)]
            if not mean.empty:
                ax.plot(mean["generation"], mean["mean"], color=color, lw=1.8)
        log_axis(ax)
        ax.set_title(f"{pop_label(population, model)}; chosen: {chosen_mode(population, model)}",
                     loc="left")
        ax.grid(axis="y", which="major", color="#E6E6E6", lw=0.6)
        ax.set_xlim(0, gone2["generation"].max())
    for ax in axes[-1]:
        ax.set_xlabel("Generations ago")
    for ax in axes[:, 0]:
        ax.set_ylabel("Ne")
    handles = [Line2D([], [], color=c, lw=1.8, label=f"{MODE_LABELS[m]}: mean of seeds")
               for m, c in MODE_COLORS.items()]
    handles.append(Line2D([], [], color="#777777", lw=0.6, alpha=0.6, label="individual seeds"))
    fig.suptitle(f"{project}: GONE2, all runs ({rate} cM/Mb)", x=0.01, ha="left")
    fig.tight_layout()
    # Below the panels, never over the data.
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3,
               frameon=False)
    save(fig, paths)


def plot_gone2_final(summary, model, generations, paths):
    first, last = generations
    summary = summary[summary["generation"].between(first, last)]
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    drawn = False
    for population in sorted(set(model["population"])):
        mode = chosen_mode(population, model)
        df = summary[(summary["population"] == population) & (summary["mode"] == mode)]
        if df.empty:
            continue
        drawn = True
        color = pop_color(population)
        sd = df["sd"].fillna(0)
        lo = (df["mean"] - sd).clip(lower=df["min"])  # keep the band positive on a log axis
        ax.fill_between(df["generation"], lo, df["mean"] + sd, color=color, alpha=0.3, lw=0)
        label = pop_label(population, model) + (", metapopulation N_T" if mode == "metapop" else "")
        ax.plot(df["generation"], df["mean"], color=color, lw=2, label=label, path_effects=OUTLINE)
    if not drawn:
        empty(ax, "no successful GONE2 runs")
    else:
        log_axis(ax)
        ax.set_xlim(first, last)
        gap = (last - first) * 0.08  # keep automatic ticks clear of the end labels
        ax.set_xticks(sorted({first, last, *(t for t in ax.get_xticks()
                                              if first + gap < t < last - gap)}))
        ax.set_xlabel("Generations ago")
        ax.set_ylabel("Effective population size (Ne)")
        ax.grid(axis="y", which="major", color="#E6E6E6", lw=0.6)
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False,
                  title="Mean of GONE2 seeds ± 1 SD", title_fontsize=8, alignment="left")
    save(fig, paths)


def plot_estimates(estimates, model, project, paths):
    populations = sorted(set(model["population"]))
    methods = [m for m in METHODS if ((estimates["mode"] == m[0]) & (estimates["method"] == m[1])).any()]
    if not populations or not methods:
        fig, ax = plt.subplots(figsize=(6.5, 2))
        empty(ax, "no Ne estimates")
        save(fig, paths)
        return
    ypos = {m: i for i, m in enumerate(reversed(methods))}
    fig, axes = plt.subplots(1, len(populations), figsize=(2.3 + 2.4 * len(populations), 0.45 * len(methods) + 1.3),
                             sharey=True, squeeze=False)
    finite = estimates[["estimate", "lo90", "hi90", "lo50", "hi50"]].to_numpy(dtype=float)
    finite = finite[np.isfinite(finite) & (finite > 0)]
    xlim = (finite.min() / 1.6, finite.max() * 1.6) if finite.size else (1, 10)
    for ax, population in zip(axes[0], populations):
        chosen = chosen_mode(population, model)
        color = pop_color(population)
        for (mode, method), y in ypos.items():
            row = estimates[(estimates["population"] == population) & (estimates["mode"] == mode)
                            & (estimates["method"] == method)]
            if row.empty or not np.isfinite(row["estimate"].iloc[0]):
                ax.text(xlim[0] * 1.15, y, "no estimate", va="center", fontsize=7, color="#888888")
                continue
            r = row.iloc[0]
            on = mode == chosen
            line = INK if on else GREY
            ax.hlines(y, r["lo90"], r["hi90"], color=line, lw=1)
            if np.isfinite(r["lo50"]) and np.isfinite(r["hi50"]):
                ax.hlines(y, r["lo50"], r["hi50"], color=line, lw=3.5)
            ax.plot(r["estimate"], y, "o", ms=7, mfc=color if on else "white",
                    mec=INK if on else GREY, mew=1)
            ax.annotate(f"{r['estimate']:.0f}", (r["estimate"], y), xytext=(0, 6),
                        textcoords="offset points", ha="center", fontsize=7,
                        color=INK if on else "#999999")
        if any(m == "metapop" for m, _ in methods) and any(m == "standard" for m, _ in methods):
            split = min(y for (m, _), y in ypos.items() if m == "standard") - 0.5
            ax.axhline(split, color="#DDDDDD", lw=0.8, ls="--")
        log_axis(ax, "x")
        ax.set_xlim(*xlim)
        ax.set_ylim(-0.6, len(methods) - 0.4)
        ax.grid(axis="x", which="major", color="#EEEEEE", lw=0.6)
        ax.set_title(f"{pop_label(population, model)}\nchosen: {MODE_LABELS[chosen]}", loc="left")
        ax.set_xlabel("Ne")
        ax.tick_params(axis="y", length=0)
    axes[0][0].set_yticks(list(ypos.values()))
    axes[0][0].set_yticklabels([m for _, m in ypos])
    handles = [
        Line2D([], [], marker="o", ls="", ms=7, mfc="#999999", mec=INK, label="chosen model"),
        Line2D([], [], marker="o", ls="", ms=7, mfc="white", mec=GREY, label="other model"),
        Line2D([], [], color=INK, lw=3.5, label="50% CI"),
        Line2D([], [], color=INK, lw=1, label="90% CI (GONE2: range across seeds)"),
    ]
    fig.suptitle(f"{project}: contemporary Ne by method", x=0.01, ha="left")
    fig.tight_layout()
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=4, frameon=False)
    save(fig, paths)


def main():
    out = snakemake.output
    project = snakemake.params["project_id"]
    gone2 = pd.read_csv(snakemake.input["gone2"], sep="\t")
    summary = pd.read_csv(snakemake.input["gone2_summary"], sep="\t")
    estimates = pd.read_csv(snakemake.input["estimates"], sep="\t")
    model = pd.read_csv(snakemake.input["model"], sep="\t")
    plot_gone2_all(gone2, summary, model, project, snakemake.params["recomb_rate"], [out["gone2_all"]])
    plot_gone2_final(summary, model, snakemake.params["generations"],
                     [out["gone2_final_pdf"], out["gone2_final_png"]])
    plot_estimates(estimates, model, project, [out["estimates_pdf"], out["estimates_png"]])


main()
