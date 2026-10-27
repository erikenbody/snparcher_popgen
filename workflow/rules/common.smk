"""Shared setup: config validation, mode gating, path helpers, target list.

Inputs are named in the config by path, not derived from wildcards, so the same
workflow runs standalone (point `vcf` at any filtered SNP VCF) or as a snpArcher
module (snpArcher passes its own output paths in). Every output lives under
`outdir` (default results/popgen) and every log under `logdir`.

Knobs that change the shape of the DAG are resolved once, here:

  use_coords    Whether sample coordinates are available at all. When false,
                TESS and every map panel drop out of the workflow and populations
                are defined from ADMIXTURE alone.
  run_ne        GONE2 + currentNe2 (rules/ne.smk).
  run_clam      GVCF-based clam analyses (rules/clam.smk; needs gvcf_dir).
  run_trackhub  UCSC track hub (rules/trackhub.smk).

Everything else is a parameter, not a branch.
"""

from pathlib import Path

from snakemake.utils import validate
from snakemake.exceptions import WorkflowError


# current_basedir, not basedir: under snpArcher's `module` directive basedir is
# snpArcher's workflow directory, not this one.
WORKFLOW_DIR = Path(workflow.current_basedir).parent
validate(config, str(WORKFLOW_DIR / "schemas" / "config.schema.yaml"))

PREFIX = config["final_prefix"]
OUTDIR = config["outdir"].rstrip("/")
LOGDIR = config["logdir"].rstrip("/")


def out(path):
    """Path under outdir. `{prefix}` is filled in; other wildcards are left alone."""
    return f"{OUTDIR}/{path}".replace("{prefix}", PREFIX)


def logp(path):
    """Path under logdir, filled in the same way as out()."""
    return f"{LOGDIR}/{path}".replace("{prefix}", PREFIX)


# --------------------------------------------------------------------- input
# The filtered SNP set is this workflow's input: no rule here produces it. It is
# not checked for existence at parse time, because inside snpArcher it is built
# by the postprocess module in the same run.
INPUT_VCF = config["vcf"]

# Contig names and lengths come from the reference .fai when one is given and
# from the VCF header's ##contig lines otherwise; see rule contig_lengths.
INPUT_FAI = config["fai"]
CONTIGS = out("prepare/{prefix}.contigs.tsv")

# --------------------------------------------------------------------- modes
CALLABLE_BED = config["callable_sites_bed"]
USE_BED = bool(CALLABLE_BED)
USE_COORDS = bool(config["use_coords"])
# TESS is inherently spatial, so it can only run with coordinates. Rather than
# fail late inside R, force it off here and say so.
RUN_TESS = bool(config["run_tess"]) and USE_COORDS
RUN_NE = bool(config["run_ne"])
RUN_CLAM = bool(config["run_clam"])
RUN_TRACKHUB = bool(config["run_trackhub"])
ASSIGN_METHOD = config["method"]

if USE_COORDS and not config["coords_sheet"]:
    raise WorkflowError(
        "use_coords is true but coords_sheet is not set. Point coords_sheet at a "
        "tab- or comma-delimited file with a header and the columns sample_id, "
        "longitude, latitude (see config/coords.tsv.example), or set "
        "use_coords: false."
    )

if ASSIGN_METHOD == "tess" and not RUN_TESS:
    reason = (
        "use_coords is false" if not USE_COORDS else "run_tess is false"
    )
    raise WorkflowError(
        f"method is 'tess' but TESS is not running because {reason}. "
        "Either enable coordinates and run_tess, or set method: admixture."
    )

if config["run_tess"] and not USE_COORDS:
    logger.warning(
        "run_tess is true but use_coords is false — TESS and the pie/map panels "
        "will be skipped. Populations will be defined from ADMIXTURE."
    )

if RUN_CLAM and not config["gvcf_dir"]:
    raise WorkflowError(
        "run_clam is true but gvcf_dir is not set. Point gvcf_dir at a directory "
        "of per-sample GVCFs named {sample}.g.vcf.gz (with .tbi indexes), or set "
        "run_clam: false."
    )

if RUN_TRACKHUB and not (USE_BED and config["trackhub_email"]):
    raise WorkflowError(
        "run_trackhub needs callable_sites_bed (for the non-callable track) and "
        "trackhub_email (UCSC requires a contact address in hub.txt)."
    )

# ADMIXTURE replicates. The first seed is the "plot seed": the replicate whose Q
# matrices are shown in the composite figure and used for population assignment.
ADMIX_K = list(range(1, int(config["admixture_kmax"]) + 1))
ADMIX_SEEDS = list(config["admixture_seeds"])
PLOT_SEED = ADMIX_SEEDS[0]

# manual_k_assignment: null means "let the data choose". The R scripts want a
# number or the string "NA".
MANUAL_K = config["manual_k_assignment"]
MANUAL_K_PARAM = "NA" if MANUAL_K is None else int(MANUAL_K)


# --------------------------------------------------------------------- paths
def gvcf_dir():
    """gvcf_dir with a {prefix} placeholder filled in."""
    return str(config["gvcf_dir"]).replace("{prefix}", PREFIX)


def depths_zarr_dir():
    """Where `clam collect` writes its depth store."""
    tpl = config["clam_depths_zarr_dir"] or out("clam/depths.zarr")
    return tpl.replace("{prefix}", PREFIX)


# ------------------------------------------------------------------ helpers
def as_bool(value):
    """Coerce a config value to bool. --config passes everything as a string."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("true", "1", "yes", "on")


def population_files():
    """{k: [samples...]} for every population_<k>.txt that has enough samples.

    Reads the checkpoint output directory. Populations below
    min_samples_per_pop are dropped: bcftools roh, GONE2 and clam all misbehave
    on a single sample, and a population of one is not a population.
    """
    pops_dir = Path(out("populations/{prefix}_populations"))
    if not pops_dir.is_dir():
        return {}
    min_n = int(config["min_samples_per_pop"])
    found = {}
    for f in sorted(pops_dir.iterdir()):
        if not (f.is_file() and f.name.startswith("population_") and f.suffix == ".txt"):
            continue
        k = f.stem.split("_", 1)[1]
        # final ancestry assignment writes population_NA.txt for samples below
        # qvalue_threshold; those are deliberately not a population.
        if not k.isdigit():
            continue
        samples = [s.strip() for s in f.read_text().splitlines() if s.strip()]
        if len(samples) < min_n:
            logger.warning(
                f"population_{k} has {len(samples)} sample(s), below "
                f"min_samples_per_pop={min_n} — skipping it."
            )
            continue
        found[k] = samples
    return found


def active_k(wc=None):
    """Population labels surviving min_samples_per_pop, gating on the checkpoint."""
    checkpoints.assign_populations.get()
    return sorted(population_files(), key=int)


# ------------------------------------------------------------------ targets
def workflow_targets():
    targets = [
        # structure: ancestry inference and the figure you use to pick K
        out("structure/{prefix}_structure_composite.pdf"),
        # scaffold QC: sex-linked and anomalous scaffolds
        out("populations/{prefix}_sex_id/{prefix}_sex_chromosome_candidates.csv"),
        out("populations/{prefix}_sex_id/{prefix}_scaffold_anomalies.csv"),
        # population assignment figure + per-population statistics
        out("populations/{prefix}_ancestry_final.pdf"),
        out("populations/{prefix}_stats.done"),
        # pairwise genetic distances
        out("prepare/{prefix}_dist.dist"),
    ]
    if RUN_TESS:
        targets.append(out("structure/{prefix}_TESS_pieplot.pdf"))
    if RUN_NE:
        targets.append(out("ne/{prefix}_ne.done"))
    if RUN_CLAM:
        targets.append(out("clam/{prefix}_manhattan.done"))
    if RUN_TRACKHUB:
        targets.append(out("trackhub/hub.txt"))
    return targets
