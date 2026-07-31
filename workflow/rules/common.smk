"""Shared setup: config validation, mode gating, path helpers, target list.

Two knobs change the shape of the DAG and are resolved once, here:

  use_coords  Whether sample coordinates are available at all. When false,
              TESS and every map panel drop out of the workflow and populations
              are defined from ADMIXTURE alone.
  run_clam    Whether the GVCF-based clam analyses run (needs gvcf_dir).

Everything else is a parameter, not a branch.
"""

from pathlib import Path

from snakemake.utils import validate
from snakemake.exceptions import WorkflowError


validate(config, os.path.join(workflow.basedir, "schemas/config.schema.yaml"))

PREFIX = config["final_prefix"]
REF = config["refgenome"]

# --------------------------------------------------------------------- input
# The filtered SNP set is this workflow's input: no rule produces it. Checking it
# here means a missing or mistyped path is reported as itself, rather than as a
# MissingInputException several rules deep.
CLEAN_SNPS = f"results/{REF}/{PREFIX}_clean_snps.vcf.gz"
if not Path(CLEAN_SNPS).exists():
    raise WorkflowError(
        f"input VCF not found: {CLEAN_SNPS}\n"
        "This workflow starts from a filtered, biallelic SNP VCF — for example "
        "snpArcher's {prefix}_clean_snps.vcf.gz. Site filtering (depth, "
        "missingness, MAF, callable regions) is the variant caller's job, not "
        "this workflow's.\n"
        "The path is built from refgenome and final_prefix, and is interpreted "
        "relative to the working directory (the --directory argument, if you "
        "passed one)."
    )

# --------------------------------------------------------------------- modes
CALLABLE_BED = config["callable_sites_bed"]
USE_BED = bool(CALLABLE_BED)
USE_COORDS = bool(config["use_coords"])
# TESS is inherently spatial, so it can only run with coordinates. Rather than
# fail late inside R, force it off here and say so.
RUN_TESS = bool(config["run_tess"]) and USE_COORDS
RUN_GONE = bool(config["run_gone"])
RUN_CLAM = bool(config["run_clam"])
ASSIGN_METHOD = config["method"]

if USE_COORDS and not config["coords_sheet"]:
    raise WorkflowError(
        "use_coords is true but coords_sheet is not set. Point coords_sheet at a "
        "tab-delimited file with a header and the columns sample_id, longitude, "
        "latitude (see config/coords.tsv.example), or set use_coords: false."
    )

if USE_COORDS and not Path(config["coords_sheet"]).exists():
    raise WorkflowError(
        f"coords_sheet not found: {config['coords_sheet']}\n"
        "The path is interpreted relative to the working directory "
        "(the --directory argument, if you passed one)."
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

# ADMIXTURE replicates. The first seed is the "plot seed": the replicate whose Q
# matrices are shown in the composite figure and used for population assignment.
ADMIX_K = list(range(1, int(config["admixture_kmax"]) + 1))
ADMIX_SEEDS = list(config["admixture_seeds"])
PLOT_SEED = ADMIX_SEEDS[0]

GONE_SEEDS = list(range(1, int(config["gone_seeds"]) + 1))

# manual_k_assignment: null means "let the data choose". The R scripts want a
# number or the string "NA".
MANUAL_K = config["manual_k_assignment"]
MANUAL_K_PARAM = "NA" if MANUAL_K is None else int(MANUAL_K)


# --------------------------------------------------------------------- paths
def gvcf_dir():
    """gvcf_dir with {refGenome}/{prefix} placeholders filled in."""
    return str(config["gvcf_dir"]).format(refGenome=REF, prefix=PREFIX)


def depths_zarr_dir():
    """Where `clam collect` writes its depth store."""
    tpl = config["clam_depths_zarr_dir"] or (
        "results/{refGenome}/clam/{prefix}/depths.zarr"
    )
    for token in ("{refGenome}", "{prefix}"):
        if token not in tpl:
            raise WorkflowError(
                f"clam_depths_zarr_dir must contain {token} (snakemake requires "
                f"both {{refGenome}} and {{prefix}} in the path). Got: {tpl!r}"
            )
    return tpl


def gone_dir():
    """Directory containing script_GONE.sh and PROGRAMMES/."""
    return config["gone_path"] or "resources/GONE/Linux"


# ------------------------------------------------------------------ helpers
def as_bool(value):
    """Coerce a config value to bool. --config passes everything as a string."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("true", "1", "yes", "on")


def population_files(ref, prefix):
    """{k: [samples...]} for every population_<k>.txt that has enough samples.

    Reads the checkpoint output directory. Populations below
    min_samples_per_pop are dropped: bcftools roh, GONE and clam all misbehave
    on a single sample, and a population of one is not a population.
    """
    pops_dir = Path(f"results/{ref}/populations/{prefix}_populations")
    if not pops_dir.is_dir():
        return {}
    min_n = int(config["min_samples_per_pop"])
    out = {}
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
        out[k] = samples
    return out


def active_k(wc):
    """Population labels surviving min_samples_per_pop, gating on the checkpoint."""
    checkpoints.assign_populations.get(**wc)
    return sorted(population_files(wc.refGenome, wc.prefix), key=int)


# ------------------------------------------------------------------ targets
def workflow_targets():
    base = f"results/{REF}"
    targets = [
        # structure: ancestry inference and the figure you use to pick K
        f"{base}/structure/{PREFIX}_structure_composite.pdf",
        # scaffold QC: sex-linked and anomalous scaffolds
        f"{base}/populations/{PREFIX}_sex_id/{PREFIX}_sex_chromosome_candidates.csv",
        f"{base}/populations/{PREFIX}_sex_id/{PREFIX}_scaffold_anomalies.csv",
        # population assignment figure + per-population statistics
        f"{base}/populations/{PREFIX}_ancestry_final.pdf",
        f"{base}/populations/{PREFIX}_stats.done",
        # pairwise genetic distances
        f"{base}/prepare/{PREFIX}_dist.dist",
    ]
    if RUN_TESS:
        targets.append(f"{base}/structure/{PREFIX}_TESS_pieplot.pdf")
    if RUN_CLAM:
        targets.append(f"{base}/clam/{PREFIX}/{PREFIX}_manhattan.done")
    return targets
