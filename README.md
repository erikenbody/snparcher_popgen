# snparcher_popgen

A Snakemake workflow that takes a joint-genotyped VCF and produces population
structure and genetic diversity results: ancestry inference, population
assignment, runs of homozygosity, nucleotide diversity, and callable-site-aware
diversity and divergence statistics from GVCFs.

It is the generalised form of the downstream analyses used by the
[California Conservation Genomics Project](https://www.ccgproject.org/), rebuilt
so that it runs on any organism with no dependency on CCGP infrastructure.
Sample coordinates come from a file you supply, and the whole workflow can run
without coordinates at all.

The natural upstream is [snpArcher](https://github.com/harvardinformatics/snpArcher),
whose output layout the input paths follow, but any joint-called VCF works.

## What it does

```
        clean biallelic SNP VCF ──┐
                                  │
  coordinate sheet ───────────────┴───► prepare ──► analysis sample set
                                              │        LD-pruned + subsampled SNPs
                                              │        distance matrices
                                              ▼
                                          structure ──► PCA
                                              │         ADMIXTURE (K, seeds)
                                              │         TESS            [coords]
                                              │         composite figure
                                              ▼
                                        populations ──► sex-linked scaffold QC
                                              │         population assignment
                                              │         ROH, F_ROH
                                              │         pi, heterozygosity
                                              │         GONE recent Ne  [optional]
                                              ▼
     per-sample GVCFs ───────────────────►  clam ─────► callable-site pi
                                                        dxy, Fst
                                                        ROH-corrected heterozygosity
```

`structure` and `populations` are separated by a deliberate manual step: run the
workflow, look at the composite figure, then pin `manual_k_assignment` in the
config and run again. Population assignment is a judgement call, and the
workflow is built to make you make it rather than to hide it.

## Inputs

All paths are relative to the working directory, which you set with
`--directory`. This lets one clone of the workflow serve many projects.

| Path | Required | Notes |
|---|---|---|
| `results/{refgenome}/{prefix}_clean_snps.vcf.gz` | yes | filtered, biallelic SNPs |
| `results/{refgenome}/data/genome/{refgenome}.fna.fai` | yes | reference index |
| a coordinate sheet | when `use_coords: true` | see below |
| a directory of GVCFs | when `run_clam: true` | `{sample}.g.vcf.gz` + `.tbi` |
| a callable-sites BED | optional | the pi/heterozygosity denominator; see below |

`{prefix}` is `final_prefix` and `{refgenome}` is `refgenome` from the config.

### Filtering is not this workflow's job

The input is an already-filtered SNP set, and there is no filtering stage. Depth,
missingness, MAF and callable-region filtering, and dropping samples that failed
QC, all belong upstream with the variant caller, where the evidence for them
lives — snpArcher already does them. Repeating that work here would mean two
places to keep consistent and two ways to get it wrong.

The workflow removes exactly two things of its own accord:

- **samples with no coordinates**, when `use_coords: true` — so that every
  analysis sees the same sample set. Allele counts and AF are recomputed
  afterwards, so sites that became invariant do not linger.
- **sex-linked and anomalous scaffolds**, detected from relative depth and
  heterozygosity in `populations`, and excluded from every per-population VCF and
  from clam.

Note that the second of those applies to the per-population statistics (ROH,
F_ROH, pi, clam) but *not* to `prepare`, so PCA, ADMIXTURE and TESS are computed
including sex-linked scaffolds. Inspect
`populations/{prefix}_sex_id/{prefix}_sex_chromosome_candidates.csv` — if a large
scaffold is flagged, it is worth knowing that the structure results saw it.

### The coordinate sheet

Tab-delimited, **with a header row**:

```
sample_id	longitude	latitude
B0624	-120.039024	37.508364
B1347	-120.421143	39.702961
```

Sample IDs must match the VCF sample names exactly. Column names are matched
loosely (`sample`/`sample_id`/`indv`/…, `lon`/`long`/`x`, `lat`/`y`) but values
are checked strictly, including a specific check for transposed longitude and
latitude — the failure mode that produces a map which looks plausible and is
wrong.

Samples missing from the sheet are dropped from the analysis by default
(`drop_samples_without_coords`), so every analysis sees the same sample set. Set
it to `false` to make missing coordinates an error instead.

### Running without coordinates

Set `use_coords: false`. TESS and every map panel leave the DAG, and populations
are defined from ADMIXTURE alone. Everything else — PCA, ADMIXTURE, population
assignment, ROH, F_ROH, pi, GONE, clam — is unaffected.

TESS is spatially explicit and genuinely cannot run without coordinates; the
workflow refuses `method: tess` in this mode rather than failing later inside R.

### The callable-sites BED

Optional, and it now does exactly one thing: it provides the denominator for pi
and per-sample heterozygosity. Without it, the denominator is the summed length
of every scaffold in the `.fai` — so positions that were never callable stay in
the denominator, and pi and heterozygosity come out biased low and are not
comparable to values computed with a callable-sites mask.

The clam module does not use it, because it derives callability per sample from
the GVCFs. So if you run clam, its statistics are properly denominated whether or
not you supply a BED. The non-clam `pi`/`het` outputs are the ones that need it,
and since clam is off by default, supplying a BED is worth it if you care about
absolute values rather than relative comparisons within one project.

## Quick start

```bash
git clone <this repo> && cd snparcher_popgen
conda env create -n snakemake -f workflow/envs/snakemake.yml   # or use your own
conda activate snakemake

# check that everything resolves, on the bundled test data
snakemake --directory .test --configfile .test/config.yaml -n

# run it
snakemake --directory .test --configfile .test/config.yaml \
          --software-deployment-method conda --cores 8
```

For your own project:

```bash
cp config/config.yaml /path/to/project/config.yaml
# edit final_prefix, refgenome, coords_sheet, gvcf_dir
snakemake --directory /path/to/project --configfile /path/to/project/config.yaml \
          --software-deployment-method conda --cores 20
```

Run these from the repository root. `profiles/default/` is picked up
automatically as the workflow profile and sets `--software-deployment-method
conda` for you, but only when Snakemake can find it — passing the flag
explicitly, as above, works from anywhere.

Conda environments are created per rule on first run. Point `--conda-prefix` at
a shared location to reuse them across projects; expect roughly 15 GB if you
build all of them, most of it the R environment used by TESS.

On a cluster, `profiles/slurm/` is a starting point (needs
`snakemake-executor-plugin-slurm`).

## Configuration

Every option, its default and what it means is documented in
`workflow/schemas/config.schema.yaml`, which also validates your config and
fills in defaults. `config/config.yaml` is a template with the defaults written
out. A minimal config is:

```yaml
final_prefix: my-project
refgenome: GCA_000000000.1
use_coords: false
```

The two optional modules are off by default: set `run_clam: true` (with
`gvcf_dir`) for callable-site statistics, and `run_gone: true` for recent Ne.
Coordinate-based analyses are on by default, so `use_coords: false` is what you
add when you have no coordinates.

The options worth thinking about before a real run:

- **`manual_k_assignment`** — the number of populations. Leave `null` for the
  first pass, read `{prefix}_structure_composite.pdf`, then set it.
- **`method`** — `admixture` or `tess`: which ancestry estimate defines
  populations.
- **`qvalue_threshold`** — samples whose top ancestry coefficient falls below
  this are treated as admixed and excluded from per-population statistics. They
  are listed in `{prefix}_populations-done.txt`, not silently dropped.
- **`min_samples_per_pop`** — populations smaller than this are skipped, with a
  warning. `bcftools roh`, GONE and clam all misbehave on one sample.
- **`admixture_kmax`** / **`admixture_seeds`** — the K range and replicates. The
  first seed is the one plotted and used for assignment, so keep it stable.
- **`clam_depths_zarr_dir`** — clam's depth store is roughly
  `n_samples × genome_length × 4` bytes before compression. Point it at bulk
  storage for large projects. It is deleted after use unless
  `clam_keep_depths_zarr` is set.

## Outputs

Under `results/{refgenome}/`:

| Path | Contents |
|---|---|
| `prepare/{prefix}_analysis_set.vcf.gz` | the sample set every analysis uses |
| `prepare/{prefix}_pruned_ld.vcf.gz` | LD-pruned SNPs |
| `prepare/{prefix}_dist.dist` | pairwise distance matrix |
| `structure/{prefix}.eigenvec` | PCA |
| `structure/admixture/Q_files/` | ancestry proportions, per K and seed |
| `structure/{prefix}_structure_composite.pdf` | **read this to choose K** |
| `structure/{prefix}_TESS_*` | TESS results and pie maps |
| `populations/{prefix}_ancestry_final.pdf` | the assignment actually used |
| `populations/{prefix}_populations/` | one sample list per population |
| `populations/{prefix}_populations-done.txt` | K, method, threshold, excluded samples |
| `populations/{prefix}_sex_id/` | candidate sex-linked and anomalous scaffolds |
| `populations/{prefix}_population_pi/` | pi and per-sample heterozygosity |
| `populations/{prefix}_merged_roh/` | ROH, F_ROH and ROH figures |
| `populations/{prefix}_GONE/plots/` | recent Ne trajectories |
| `clam/{prefix}/{prefix}_summary_per_pop.csv` | callable-site pi, heterozygosity in/out of ROH |
| `clam/{prefix}/{prefix}_summary_pairwise.csv` | dxy and Hudson's Fst per population pair |

Sex-linked and anomalous scaffolds are excluded from every per-population VCF
and from clam. Inspect `{prefix}_sex_chromosome_candidates.csv` and
`{prefix}_scaffold_anomalies.csv` — the detection is a heuristic based on
relative depth and heterozygosity, and on small or sparse datasets it flags more
than it should.

## Test data

`.test/` holds a 17-sample, 78k-SNP subset of a black bear
(*Ursus americanus*) dataset, with **synthetic** coordinates in two loose
clusters. The coordinates are invented, not real localities. Results are
meaningless; the purpose is to check that the workflow executes.

```bash
# with coordinates (runs TESS, so the first run builds the large R environment)
snakemake --directory .test --configfile .test/config.yaml \
          --software-deployment-method conda --cores 8

# without coordinates — faster, and the quickest end-to-end check
snakemake --directory .test --configfile .test/config_nocoords.yaml \
          --software-deployment-method conda --cores 8
```

Both paths have been run to completion on this data. Neither test config enables
clam or GONE: per-sample GVCFs are far too large to ship in a repository, and the
test slice has nowhere near enough linked variants for LD-based Ne inference to
mean anything. To exercise clam, point `gvcf_dir` at real GVCFs for the samples
in your own dataset.

## Software

Every rule declares its own conda environment in `workflow/envs/`, so the only
thing you install yourself is Snakemake (>= 7.32; developed against 9.13).

Three R packages are not on conda and are installed from GitHub on first use by
the `algatr_env` rule: [algatr](https://github.com/TheWangLab/algatr),
[wingen](https://github.com/AnushaPB/wingen) and
[tess3r](https://github.com/bcm-uga/TESS3_encho_sen). They are only needed when
`run_tess` is on — and be warned that this step compiles dozens of R
dependencies from source and takes a long time on a first run. It needs network
access. If you do not need TESS, `run_tess: false` avoids it entirely.

GONE is cloned and built by the workflow when `run_gone` is on, or point
`gone_path` at an existing installation.

`clam` comes from its conda package.

## Not included

Deliberately left out of this general version, since they depend on
California-specific data or CCGP infrastructure:

- **Environmental association** (RDA, LFMM, wingen) — needs regional climate
  rasters and a coastline shapefile
- **Mutation load** (SIFT, phyloP) — needs per-reference annotation tables that
  cannot be derived within the workflow
- **UCSC track hubs** — needs hosting and a hub-registered email
- **Contig renaming** — needs an assembly report mapping accessions to
  chromosome names

## Citing the tools

This workflow orchestrates other people's software. If you publish results,
cite: bcftools/samtools, PLINK, ADMIXTURE, VCFtools, scikit-allel, TESS3 and
algatr, GONE, and clam.
