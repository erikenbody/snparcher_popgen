# snparcher_popgen

A Snakemake workflow that takes a joint-genotyped VCF and produces population
structure and genetic diversity results: ancestry inference, population
assignment, runs of homozygosity, nucleotide diversity, recent and contemporary
effective population size, callable-site-aware diversity and divergence
statistics from GVCFs, and a UCSC track hub.

It is the generalised form of the downstream analyses used by the
[California Conservation Genomics Project](https://www.ccgproject.org/), rebuilt
so that it runs on any organism with no dependency on CCGP infrastructure.
Sample coordinates come from a file you supply, and the whole workflow can run
without coordinates at all.

The natural upstream is [snpArcher](https://github.com/harvardinformatics/snpArcher):
it runs standalone on snpArcher's outputs (or any joint-called VCF), and is laid
out to be imported as a snpArcher module — every input is a path in the config,
and every output lands under one directory.

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
                                              ▼
                                              ne ─────► GONE2 recent Ne        [run_ne]
                                                        currentNe2 current Ne
     per-sample GVCFs ───────────────────►  clam ─────► callable-site pi       [run_clam]
                                                        dxy, Fst
                                                        ROH-corrected heterozygosity
     clean SNPs + callable sites ───────►  trackhub ─► UCSC track hub          [run_trackhub]
```

`structure` and `populations` are separated by a deliberate manual step: run the
workflow, look at the composite figure, then pin `manual_k_assignment` in the
config and run again. Population assignment is a judgement call, and the
workflow is built to make you make it rather than to hide it.

## Inputs

Every input is named by path in the config. Relative paths are relative to the
working directory, which you set with `--directory`; this lets one clone of the
workflow serve many projects.

| Config key | Required | Notes |
|---|---|---|
| `vcf` | yes | filtered, biallelic SNPs, bgzipped and indexed |
| `fai` | no | contig lengths; read from the VCF's `##contig` header lines when absent |
| `coords_sheet` | when `use_coords: true` | see below |
| `callable_sites_bed` | no | the pi/heterozygosity denominator and the track hub's non-callable track; see below |
| `gvcf_dir` | when `run_clam: true` | `{sample}.g.vcf.gz` + `.tbi` |

You do not need the reference genome or its index. Contig names and lengths —
all the workflow ever reads from a `.fai` — come from the VCF header, which
GATK, Sentieon and bcftools fill with every reference contig, including ones
that end up with no variants. The `contig_lengths` rule fails loudly if the
header has no lengths, or if a contig carrying variants is missing from it.

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
  heterozygosity in `populations`, and excluded from every per-population VCF,
  from the Ne analyses and from clam.

Note that the second of those applies to the per-population statistics (ROH,
F_ROH, pi, Ne, clam) but *not* to `prepare`, so PCA, ADMIXTURE and TESS are
computed including sex-linked scaffolds. Inspect
`populations/{prefix}_sex_id/{prefix}_sex_chromosome_candidates.csv` — if a large
scaffold is flagged, it is worth knowing that the structure results saw it.

### The coordinate sheet

A header row, then one row per sample. Comma-delimited if the file name ends in
`.csv`, tab-delimited otherwise:

```
sample_id	longitude	latitude
B0624	-120.039024	37.508364
B1347	-120.421143	39.702961
```

Sample IDs must match the VCF sample names exactly. Column names are matched
loosely (`sample`/`sample_id`/`indv`/…, `lon`/`long`/`x`, `lat`/`y`) — so
snpArcher's `sample_metadata.csv` works as is — but values are checked strictly,
including a specific check for transposed longitude and latitude, the failure
mode that produces a map which looks plausible and is wrong.

Samples missing from the sheet are dropped from the analysis by default
(`drop_samples_without_coords`), so every analysis sees the same sample set. Set
it to `false` to make missing coordinates an error instead.

### Running without coordinates

Set `use_coords: false`. TESS and every map panel leave the DAG, and populations
are defined from ADMIXTURE alone. Everything else — PCA, ADMIXTURE, population
assignment, ROH, F_ROH, pi, Ne, clam — is unaffected.

TESS is spatially explicit and genuinely cannot run without coordinates; the
workflow refuses `method: tess` in this mode rather than failing later inside R.

### The callable-sites BED

Optional. It provides the denominator for pi and per-sample heterozygosity.
Without it, the denominator is the summed length of every contig — so positions
that were never callable stay in the denominator, and pi and heterozygosity come
out biased low and are not comparable to values computed with a callable-sites
mask. The track hub also uses it for its non-callable-regions track.

The clam module does not use it, because it derives callability per sample from
the GVCFs. So if you run clam, its statistics are properly denominated whether or
not you supply a BED.

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
# edit vcf, final_prefix, coords_sheet, callable_sites_bed, ...
snakemake --directory /path/to/project --configfile /path/to/project/config.yaml \
          --software-deployment-method conda --cores 20
```

Run these from the repository root. `profiles/default/` is picked up
automatically as the workflow profile and sets `--software-deployment-method
conda` for you, but only when Snakemake can find it — passing the flag
explicitly, as above, works from anywhere.

Conda environments are created per rule on first run. Point `--conda-prefix` at
a shared location to reuse them across projects, and point `ne_tools_dir` at a
shared location to compile GONE2 and currentNe2 once.

On a cluster, `profiles/slurm/` is a starting point (needs
`snakemake-executor-plugin-slurm`).

## Configuration

Every option, its default and what it means is documented in
`workflow/schemas/config.schema.yaml`, which also validates your config and
fills in defaults. `config/config.yaml` is a template with the defaults written
out. A minimal config is:

```yaml
vcf: path/to/clean_snps.vcf.gz
use_coords: false
```

The optional analyses are off by default: `run_ne: true` for GONE2/currentNe2,
`run_clam: true` (with `gvcf_dir`) for callable-site statistics, and
`run_trackhub: true` (with `refgenome`, `callable_sites_bed`, `trackhub_email`)
for a track hub. Coordinate-based analyses are on by default, so
`use_coords: false` is what you add when you have no coordinates.

The options worth thinking about before a real run:

- **`manual_k_assignment`** — the number of populations. Leave `null` for the
  first pass, read `{prefix}_structure_composite.pdf`, then set it.
- **`method`** — `admixture` or `tess`: which ancestry estimate defines
  populations.
- **`qvalue_threshold`** — samples whose top ancestry coefficient falls below
  this are treated as admixed and excluded from per-population statistics. They
  are listed in `{prefix}_populations-done.txt`, not silently dropped.
- **`min_samples_per_pop`** — populations smaller than this are skipped, with a
  warning. `bcftools roh`, the Ne tools and clam all misbehave on one sample.
- **`admixture_kmax`** / **`admixture_seeds`** — the K range and replicates. The
  first seed is the one plotted and used for assignment, so keep it stable.
- **`ne_recomb_rate`** — the constant cM/Mb both Ne tools use in place of a
  genetic map. Ne estimates scale with it, so think about it per taxon.
- **`clam_depths_zarr_dir`** — clam's depth store is roughly
  `n_samples × genome_length × 4` bytes before compression. Point it at bulk
  storage for large projects. It is deleted after use unless
  `clam_keep_depths_zarr` is set.

### Effective population size

`run_ne: true` adds, per population, [GONE2](https://github.com/esrud/GONE2)
(recent Ne trajectory) and [currentNe2](https://github.com/esrud/currentNe2)
(contemporary Ne). Both are compiled from pinned commits into `ne_tools_dir`,
with their locus and sample-size limits (`ne_max_snps`, `ne_max_ind`) baked in;
each population's SNPs are thinned at random to fit. Only contigs of at least
`ne_min_contig_size` are used.

Each tool runs in the modes listed in `ne_gone2_modes` / `ne_currentne2_modes`:
`standard` (one panmictic population) and `metapop` (`-x`, a structured
sample — useful when structure persists within an assigned population).

Runs never fail the workflow. Small or odd populations often cannot support an
estimate — GONE2 then writes only its `_STATS` file — so every run records its
exit status, and `{prefix}_ne_status.tsv` lists how each ended. Filter on it.
`snakemake ne_all` runs just this part once populations exist.

## Outputs

Under `outdir` (default `results/popgen/`):

| Path | Contents |
|---|---|
| `prepare/{prefix}.contigs.tsv` | contig names and lengths used throughout |
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
| `ne/{prefix}_gone2_Ne.tsv`, `ne/{prefix}_gone2_Ne.png` | recent Ne trajectories |
| `ne/{prefix}_currentne2.tsv` | contemporary Ne estimates |
| `ne/{prefix}_ne_status.tsv` | how every Ne run ended |
| `clam/{prefix}_summary_per_pop.csv` | callable-site pi, heterozygosity in/out of ROH |
| `clam/{prefix}_summary_pairwise.csv` | dxy and Hudson's Fst per population pair |
| `trackhub/` | a self-contained UCSC track hub; copy it to a web server |

Logs are under `logdir` (default `logs/popgen/`).

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
# with coordinates, TESS, Ne and the track hub
snakemake --directory .test --configfile .test/config.yaml \
          --software-deployment-method conda --cores 8

# without coordinates, using the .fai rather than the VCF header — the quickest
# end-to-end check
snakemake --directory .test --configfile .test/config_nocoords.yaml \
          --software-deployment-method conda --cores 8
```

The test configs do not enable clam: per-sample GVCFs are far too large to ship
in a repository. To exercise clam, point `gvcf_dir` at real GVCFs for the samples
in your own dataset. The test populations are far too small for meaningful Ne;
the Ne runs there only check that the tools build, run and are summarised.

## Software

Every rule declares its own conda environment in `workflow/envs/`, so the only
thing you install yourself is Snakemake (>= 7.32; developed against 9.13).

Three tools are not on conda and are built on first use, which needs network
access:

- [tess3r](https://github.com/bcm-uga/TESS3_encho_sen), installed into the TESS
  environment by the `tess_env` rule — a small Rcpp package, so a quick
  compile. Only needed when `run_tess` is on.
- [GONE2](https://github.com/esrud/GONE2) and
  [currentNe2](https://github.com/esrud/currentNe2), compiled into
  `ne_tools_dir`. Only needed when `run_ne` is on.

Each is pinned to a commit (`tess3r_ref`, `ne_gone2_ref`, `ne_currentne2_ref`).

`clam` comes from its conda package.

## Not included

Deliberately left out of this general version, since they depend on
California-specific data or other infrastructure:

- **Environmental association** (RDA, LFMM, wingen) — needs regional climate
  rasters and a coastline shapefile
- **Mutation load** (SIFT, phyloP) — needs per-reference annotation tables that
  cannot be derived within the workflow
- **Contig renaming** — needs an assembly report mapping accessions to
  chromosome names

## Citing the tools

This workflow orchestrates other people's software. If you publish results,
cite: bcftools/samtools, PLINK, ADMIXTURE, VCFtools, scikit-allel, TESS3, GONE2,
currentNe2 and clam.

## License

MIT — see [LICENSE](LICENSE).
