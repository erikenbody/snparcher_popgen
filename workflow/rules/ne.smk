"""Effective population size: recent trajectories (GONE2) and contemporary Ne
(currentNe2), per population. Included only when run_ne is true, and runnable
on its own afterwards with `snakemake ne_all`.

Both tools read an uncompressed VCF with unphased genotypes and turn physical
positions into genetic ones with a constant recombination rate (ne_recomb_rate,
cM/Mb). They are compiled once into ne_tools_dir — point it at a shared location
when running many projects — with their loci/individual limits baked in, which
is why each population's VCF is thinned to at most ne_max_snps. currentNe2
compares every pair of SNPs, so its input is thinned further, to
ne_currentne2_snps.

Runs never fail the workflow. GONE2 deliberately writes only its _STATS file
when the data cannot support an estimate, and small or structured populations
are expected to fail sometimes; each run records its exit status, and the
summary table reports it, so results can be filtered afterwards.

GONE2 is run once per seed in ne_gone2_seeds (its README's way of gauging
run-to-run variation, and what the original CCGP GONE runs did with 10 seeds);
the summary reports the mean, SD and range across seeds.

Modes (ne_gone2_modes / ne_currentne2_modes):
  standard  one panmictic population
  metapop   -x: the sample is from a structured population (GONE2: subsampled
            to ne_gone2_metapop_snps, as its README recommends, because -x is
            slow; currentNe2: two equal subpopulations, so -k is not used)

Which model applies is decided per population from the -x fits (see
summarize_ne.py): metapop when the fitted Fst is at least
ne_structure_min_fst, standard otherwise.

Outputs, under {outdir}/ne/:
  vcf/population_{k}.vcf                         thinned GT-only input (temp)
  vcf/population_{k}_currentne2.vcf              thinned further for currentNe2 (temp)
  gone2/population_{k}_{mode}/seed_{seed}/       GONE2 outputs + exit status
  currentne2/population_{k}_{mode}/              currentNe2 output + exit status
  {prefix}_gone2_Ne.tsv                          every Ne trajectory (per seed), long format
  {prefix}_gone2_Ne_summary.tsv                  mean/SD/median/range across seeds
  {prefix}_gone2_Ne_change.tsv                   Ne at generations 1/10/50/100/150
  {prefix}_currentne2.tsv                        every currentNe2 output value
  {prefix}_ne_estimates.tsv                      contemporary Ne per method, with CIs
  {prefix}_ne_model.tsv                          standard or metapop per population, and why
  {prefix}_ne_status.tsv                         every run and how it ended
  {prefix}_gone2_Ne.png                          every GONE2 run, per population
  {prefix}_gone2_Ne_final.{pdf,png}              chosen-model trajectories (publication)
  {prefix}_ne_estimates.{pdf,png}                contemporary Ne by method
"""

NE_TOOLS = config["ne_tools_dir"].rstrip("/")
NE_MAXLOCI = int(config["ne_max_snps"])
NE_MAXIND = int(config["ne_max_ind"])
# The build directory names the commit and limits, so builds with different
# settings never overwrite each other in a shared ne_tools_dir.
GONE2_BIN = (
    f"{NE_TOOLS}/gone2_{config['ne_gone2_ref'][:10]}_{NE_MAXLOCI}_{NE_MAXIND}/gone2"
)
CURRENTNE2_BIN = (
    f"{NE_TOOLS}/currentne2_{config['ne_currentne2_ref'][:10]}_{NE_MAXLOCI}_{NE_MAXIND}/currentne2"
)
NE_GONE2_SEEDS = [int(s) for s in config["ne_gone2_seeds"]]
NE_MODES = {
    "gone2": list(config["ne_gone2_modes"]),
    "currentne2": list(config["ne_currentne2_modes"]),
}


wildcard_constraints:
    ne_mode="standard|metapop",
    ne_gone2_seed=r"\d+",


rule ne_build_gone2:
    """Clone GONE2 at ne_gone2_ref and compile it with the configured limits.

    Built in a temp dir and moved into place with a no-clobber rename, so two
    projects building into the same shared ne_tools_dir cannot corrupt it.
    """
    output:
        binary=GONE2_BIN,
    params:
        ref=config["ne_gone2_ref"],
        maxloci=NE_MAXLOCI,
        maxind=NE_MAXIND,
    conda:
        "../envs/ne.yml"
    log:
        logp("ne/build_gone2.log"),
    shell:
        """
        exec > {log} 2>&1
        TMP=$(mktemp -d)
        trap 'rm -rf $TMP' EXIT
        git clone --quiet https://github.com/esrud/GONE2.git $TMP/src
        git -C $TMP/src checkout --quiet {params.ref}
        make -C $TMP/src gone CC="${{CXX:-g++}}" MAXLOCI={params.maxloci} MAXIND={params.maxind}
        mkdir -p $(dirname {output.binary})
        mv -n $TMP/src/gone2 {output.binary}
        """


rule ne_build_currentne2:
    """Clone currentNe2 at ne_currentne2_ref and compile it.

    Its limits are #defines in the source rather than make variables, so they
    are set to match GONE2's before compiling.
    """
    output:
        binary=CURRENTNE2_BIN,
    params:
        ref=config["ne_currentne2_ref"],
        maxloci=NE_MAXLOCI,
        maxind=NE_MAXIND,
    conda:
        "../envs/ne.yml"
    log:
        logp("ne/build_currentne2.log"),
    shell:
        """
        exec > {log} 2>&1
        TMP=$(mktemp -d)
        trap 'rm -rf $TMP' EXIT
        git clone --quiet https://github.com/esrud/currentNe2.git $TMP/src
        git -C $TMP/src checkout --quiet {params.ref}
        sed -i -e 's/^#define MAXLOCI .*/#define MAXLOCI {params.maxloci}/' \
               -e 's/^#define MAXIND .*/#define MAXIND {params.maxind}/' \
               $TMP/src/currentne2.cpp
        make -C $TMP/src currentne CXX="${{CXX:-g++}}"
        mkdir -p $(dirname {output.binary})
        mv -n $TMP/src/currentne2 {output.binary}
        """


rule ne_contigs:
    """Contigs used for Ne: at least ne_min_contig_size bp, largest first, at
    most ne_max_contigs (both tools cap the number of chromosomes)."""
    input:
        contigs=CONTIGS,
        exclude_scaffolds=out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt"),
    output:
        contigs=out("ne/{prefix}_ne_contigs.txt"),
    params:
        min_size=config["ne_min_contig_size"],
        max_n=config["ne_max_contigs"],
    log:
        logp("ne/{prefix}_ne_contigs.log"),
    shell:
        """
        exec 2> {log}
        # Sex-linked/anomalous scaffolds are already absent from the
        # population VCFs; leave them out of the list too so it says what was used.
        # FILENAME, not NR==FNR: with no scaffold excluded the first file is
        # empty, and NR==FNR would then read the contigs as exclusions
        # (17-Metacarcinus got no Ne contigs that way).
        awk -v m={params.min_size} -v ex_file={input.exclude_scaffolds} \
            'FILENAME == ex_file {{ ex[$1] = 1; next }} !($1 in ex) && $2 >= m' \
            {input.exclude_scaffolds} {input.contigs} \
          | sort -k2,2nr | head -n {params.max_n} | cut -f1 > {output.contigs}
        echo "$(wc -l < {output.contigs}) contigs >= {params.min_size} bp" >&2
        """


rule ne_vcf:
    """One population's Ne input: biallelic single-base SNPs that are
    polymorphic within the population, on the Ne contigs, GT only, thinned at
    random to at most ne_max_snps, and uncompressed (neither tool reads gzip)."""
    input:
        vcf=out("populations/{prefix}_population_vcf/population_{k}.vcf.gz"),
        csi=out("populations/{prefix}_population_vcf/population_{k}.vcf.gz.csi"),
        contigs=out("ne/{prefix}_ne_contigs.txt"),
    output:
        vcf=temp(out("ne/vcf/population_{k}.vcf")),
    params:
        max_snps=NE_MAXLOCI,
        seed=config["ne_seed"],
    # bcftools 1.16's `annotate -x INFO,^FORMAT/GT` aborts with heap corruption
    # on some populations (17-Megalorchestia, 17-Metacarcinus, 17-Panulirus);
    # 1.20+ does not. trackhub.yml carries bcftools 1.21.
    conda:
        "../envs/trackhub.yml"
    log:
        logp("ne/{prefix}_vcf_{k}.log"),
    shell:
        """
        exec 2> {log}
        TMP=$(mktemp -d)
        trap 'rm -rf $TMP' EXIT

        bcftools view -r $(paste -sd, {input.contigs}) -m2 -M2 -v snps \
            -e 'strlen(REF)!=1 || strlen(ALT)!=1 || AC==0 || AC==AN' {input.vcf} -Ou \
          | bcftools annotate -x INFO,^FORMAT/GT -Oz -o $TMP/all.vcf.gz
        N=$(bcftools index -n $TMP/all.vcf.gz 2>/dev/null || bcftools view -H $TMP/all.vcf.gz | wc -l)
        WINDOW=$(( (N + {params.max_snps} - 1) / {params.max_snps} ))
        echo "population {wildcards.k}: $N SNPs; keeping 1 in $WINDOW" >&2

        if [ "$WINDOW" -gt 1 ]; then
            bcftools +prune --random-seed {params.seed} -w $WINDOW -n 1 -N rand \
              -Ov -o {output.vcf} $TMP/all.vcf.gz
        else
            bcftools view -Ov -o {output.vcf} $TMP/all.vcf.gz
        fi
        """


rule gone2:
    """Recent Ne trajectory for one seed. Always succeeds; the run's exit
    status is recorded."""
    input:
        vcf=out("ne/vcf/population_{k}.vcf"),
        binary=ancient(GONE2_BIN),
    output:
        rundir=directory(out("ne/gone2/population_{k}_{ne_mode}/seed_{ne_gone2_seed}")),
    params:
        rate=config["ne_recomb_rate"],
        min_cm=20.5,
        mode_args=lambda wc: (
            f"-x -s {config['ne_gone2_metapop_snps']}" if wc.ne_mode == "metapop" else ""
        ),
    threads: config["ne_threads"]
    resources:
        mem_mb=lambda wc, attempt: 8000 * attempt,
    conda:
        "../envs/ne.yml"
    log:
        logp("ne/{prefix}_gone2_{k}_{ne_mode}_seed{ne_gone2_seed}.log"),
    benchmark:
        logp("ne/{prefix}_gone2_{k}_{ne_mode}_seed{ne_gone2_seed}.benchmark.tsv")
    shell:
        """
        mkdir -p {output.rundir}
        # GONE2 aborts if any chromosome spans 20 cM or less, so keep only
        # contigs whose SNPs span more than that at the configured rate.
        VCF={output.rundir}/input.vcf
        awk -v rate={params.rate} -v min_cm={params.min_cm} '
            NR == FNR {{
                if ($0 !~ /^#/) {{
                    if (!($1 in lo) || $2 < lo[$1]) lo[$1] = $2
                    if ($2 > hi[$1]) hi[$1] = $2
                }}
                next
            }}
            /^#/ {{ print; next }}
            (hi[$1] - lo[$1]) / 1e6 * rate > min_cm
        ' {input.vcf} {input.vcf} > $VCF
        echo "contigs spanning > {params.min_cm} cM: $(grep -v '^#' $VCF | cut -f1 | uniq | wc -l)" > {log}
        set +e
        {input.binary} -g 0 -r {params.rate} -S {wildcards.ne_gone2_seed} -t {threads} \
            {params.mode_args} -o {output.rundir}/population_{wildcards.k} \
            $VCF >> {log} 2>&1
        echo $? > {output.rundir}/exit_status
        rm -f $VCF
        exit 0
        """


rule currentne2_vcf:
    """currentNe2's input: the population's Ne VCF thinned at random (seeded,
    unlike currentNe2's own -s) to at most ne_currentne2_snps, because its run
    time grows with the square of the SNP count."""
    input:
        vcf=out("ne/vcf/population_{k}.vcf"),
    output:
        vcf=temp(out("ne/vcf/population_{k}_currentne2.vcf")),
    params:
        max_snps=config["ne_currentne2_snps"],
        seed=config["ne_seed"],
    conda:
        "../envs/bcftools.yml"
    log:
        logp("ne/{prefix}_currentne2_vcf_{k}.log"),
    shell:
        """
        exec 2> {log}
        N=$(grep -vc '^#' {input.vcf})
        WINDOW=$(( (N + {params.max_snps} - 1) / {params.max_snps} ))
        echo "population {wildcards.k}: $N SNPs; keeping 1 in $WINDOW for currentNe2" >&2
        if [ "$WINDOW" -gt 1 ]; then
            bcftools +prune --random-seed {params.seed} -w $WINDOW -n 1 -N rand \
              -Ov -o {output.vcf} {input.vcf}
        else
            cp {input.vcf} {output.vcf}
        fi
        """


rule currentne2:
    """Contemporary Ne. Always succeeds; the run's exit status is recorded."""
    input:
        vcf=out("ne/vcf/population_{k}_currentne2.vcf"),
        binary=ancient(CURRENTNE2_BIN),
    output:
        rundir=directory(out("ne/currentne2/population_{k}_{ne_mode}")),
    params:
        rate=config["ne_recomb_rate"],
        mode_args=lambda wc: (
            "-x" if wc.ne_mode == "metapop"
            else (f"-k {config['ne_currentne2_k']}" if config["ne_currentne2_k"] is not None else "")
        ),
    threads: config["ne_currentne2_threads"]
    resources:
        mem_mb=lambda wc, attempt: (NE_MAXLOCI * NE_MAXIND // 1_000_000 + 4000) * attempt,
    conda:
        "../envs/ne.yml"
    log:
        logp("ne/{prefix}_currentne2_{k}_{ne_mode}.log"),
    benchmark:
        logp("ne/{prefix}_currentne2_{k}_{ne_mode}.benchmark.tsv")
    shell:
        # currentNe2 writes its progress file next to its input, so each run
        # reads the VCF through a link in its own directory; otherwise the
        # standard and metapop runs overwrite each other's progress file.
        """
        mkdir -p {output.rundir}
        ln -sf $(realpath {input.vcf}) {output.rundir}/input.vcf
        set +e
        {input.binary} -r {params.rate} -t {threads} {params.mode_args} \
            -o {output.rundir}/population_{wildcards.k}_currentNe2_OUTPUT.txt \
            {output.rundir}/input.vcf > {log} 2>&1
        echo $? > {output.rundir}/exit_status
        rm -f {output.rundir}/input.vcf
        exit 0
        """


def ne_run_dirs(wc):
    ks = active_k()
    return {
        "gone2": expand(
            out("ne/gone2/population_{k}_{ne_mode}/seed_{ne_gone2_seed}"),
            k=ks, ne_mode=NE_MODES["gone2"], ne_gone2_seed=NE_GONE2_SEEDS,
        ),
        "currentne2": expand(
            out("ne/currentne2/population_{k}_{ne_mode}"), k=ks, ne_mode=NE_MODES["currentne2"]
        ),
    }


rule ne_summarize:
    """Collect every run into tables and choose each population's model."""
    input:
        unpack(ne_run_dirs),
        contigs=out("ne/{prefix}_ne_contigs.txt"),
    output:
        gone2=out("ne/{prefix}_gone2_Ne.tsv"),
        gone2_summary=out("ne/{prefix}_gone2_Ne_summary.tsv"),
        gone2_change=out("ne/{prefix}_gone2_Ne_change.tsv"),
        currentne2=out("ne/{prefix}_currentne2.tsv"),
        estimates=out("ne/{prefix}_ne_estimates.tsv"),
        model=out("ne/{prefix}_ne_model.tsv"),
        status=out("ne/{prefix}_ne_status.tsv"),
    params:
        project_id=PREFIX,
        min_fst=config["ne_structure_min_fst"],
        gone2_modes=NE_MODES["gone2"],
        currentne2_modes=NE_MODES["currentne2"],
    conda:
        "../envs/plot.yml"
    log:
        logp("ne/{prefix}_summarize.log"),
    script:
        "../scripts/summarize_ne.py"


rule ne_plot:
    """Plot every GONE2 run, the chosen-model trajectories and the
    contemporary Ne estimates."""
    input:
        gone2=out("ne/{prefix}_gone2_Ne.tsv"),
        gone2_summary=out("ne/{prefix}_gone2_Ne_summary.tsv"),
        estimates=out("ne/{prefix}_ne_estimates.tsv"),
        model=out("ne/{prefix}_ne_model.tsv"),
    output:
        gone2_all=out("ne/{prefix}_gone2_Ne.png"),
        gone2_final_pdf=out("ne/{prefix}_gone2_Ne_final.pdf"),
        gone2_final_png=out("ne/{prefix}_gone2_Ne_final.png"),
        estimates_pdf=out("ne/{prefix}_ne_estimates.pdf"),
        estimates_png=out("ne/{prefix}_ne_estimates.png"),
    params:
        project_id=PREFIX,
        recomb_rate=config["ne_recomb_rate"],
        generations=config["ne_plot_generations"],
    conda:
        "../envs/plot.yml"
    log:
        logp("ne/{prefix}_plot.log"),
    script:
        "../scripts/plot_ne.py"


rule ne_all:
    """Target for running just the Ne analyses: `snakemake ne_all`."""
    input:
        out("ne/{prefix}_ne_model.tsv"),
        out("ne/{prefix}_ne_status.tsv"),
        out("ne/{prefix}_gone2_Ne_final.pdf"),
        out("ne/{prefix}_ne_estimates.pdf"),
    output:
        touch(out("ne/{prefix}_ne.done")),
