"""Per-population statistics, downstream of the population assignment checkpoint.

Order of operations:
  1. sex-linked and anomalous scaffolds are identified from depth and
     heterozygosity, and excluded from every per-population VCF and from the
     F_ROH and pi denominators
  2. samples are assigned to populations (checkpoint: the number of populations
     is not known until the ancestry coefficients have been inspected)
  3. per population: ROH -> F_ROH, pi and heterozygosity, per-individual pi and
     ROH figures

Recent and contemporary Ne (GONE2, currentNe2) build on the per-population VCFs
made here but live in rules/ne.smk, behind run_ne.

Outputs, under {outdir}/populations/:
  {prefix}_sex_id/                          scaffold QC: candidates, anomalies, plots
  {prefix}_ancestry_final.pdf               the assignment actually used
  {prefix}_populations/population_{k}.txt   sample lists, one per population
  {prefix}_population_vcf/                  per-population VCFs
  {prefix}_merged_roh/                      merged ROH, F_ROH, ROH figures
  {prefix}_population_pi/                   pi and heterozygosity per population
  {prefix}_stats.done                       sentinel: all per-population work done
"""


# ---------------------------------------------------------------- scaffold QC
rule sex_id_extract:
    """Genotype counts in evenly spaced windows of the analysis set, for the
    scaffold QC below.

    Windows of sex_id_window bp, one centred in every sex_id_step bp, on scaffolds
    >= 1 Mb (shorter ones are excluded outright): ~10% of the genome by default.
    All SNPs in those windows, not the LD-pruned set, because pruning thins
    low-recombination scaffolds harder and so distorts het per bp between
    scaffolds. The window list is split across threads, each a bcftools query |
    awk stream; mawk is ~1.5x faster than gawk here and is used when present.
    """
    input:
        vcf=out("prepare/{prefix}_analysis_set.vcf.gz"),
        idx=out("prepare/{prefix}_analysis_set.vcf.gz.csi"),
        contigs=CONTIGS,
    output:
        # per window x sample: called, het, summed DP (see sex_id_window_counts.awk)
        counts=out("populations/{prefix}_sex_id/{prefix}_window_counts.tsv.gz"),
        samples=out("populations/{prefix}_sex_id/{prefix}_window_samples.txt"),
    params:
        window=config["sex_id_window"],
        step=config["sex_id_step"],
        awk_script=str(WORKFLOW_DIR / "scripts" / "sex_id_window_counts.awk"),
    threads: 4
    resources:
        mem_mb=4000,
    conda:
        "../envs/bcftools.yml"
    log:
        logp("populations/{prefix}_sex_id_extract.log"),
    shell:
        """
        exec 2> {log}
        AWK=$(command -v mawk || command -v awk)
        TMP=$(mktemp -d {output.counts}.tmp.XXXXXX)
        trap 'rm -rf "$TMP"' EXIT
        bcftools query -l {input.vcf} > {output.samples}
        # 1-based inclusive regions, one window centred in each step-sized block
        awk -v w={params.window} -v s={params.step} 'BEGIN {{OFS = "\\t"; o = int((s - w) / 2)}}
            $2 >= 1000000 {{for (b = 0; b + o + w <= $2; b += s) print $1, b + o + 1, b + o + w}}' \
            {input.contigs} > "$TMP/regions.tsv"
        echo "$(wc -l < "$TMP/regions.tsv") windows, $AWK, {threads} streams" >&2
        if [ ! -s "$TMP/regions.tsv" ]; then
            gzip < /dev/null > {output.counts}
            exit 0
        fi
        split -n l/{threads} -d "$TMP/regions.tsv" "$TMP/part."
        pids=()
        for r in "$TMP"/part.*; do
            [ -s "$r" ] || continue  # fewer windows than threads
            ( bcftools query -R "$r" -f '%CHROM\\t%POS[\\t%GT:%DP]\\n' {input.vcf} \
                | "$AWK" -v step={params.step} -f {params.awk_script} > "$r.counts" ) &
            pids+=($!)
        done
        for p in "${{pids[@]}}"; do wait "$p"; done
        cat "$TMP"/part.*.counts | gzip > {output.counts}
        """


rule sex_id_analyze:
    """Scaffolds to exclude: zero-het, low-het-in-all and half-depth (X/Z)
    scaffolds, plus every scaffold < 1 Mb (scripts/identify_sex_linked_scaffolds.py
    documents the rules).

    The exclude list this produces is applied to every per-population VCF, to the
    F_ROH and pi denominators, to Ne and to clam, so that sex-linked scaffolds,
    haplotigs and collapsed sequence do not distort diversity estimates.
    """
    input:
        counts=out("populations/{prefix}_sex_id/{prefix}_window_counts.tsv.gz"),
        samples=out("populations/{prefix}_sex_id/{prefix}_window_samples.txt"),
        contigs=CONTIGS,
    output:
        # scaffolds >= 1 Mb that were excluded, with the call and its evidence
        report=out("populations/{prefix}_sex_id/{prefix}_sex_chromosome_candidates.csv"),
        # every scaffold >= 1 Mb: call, excluded, evidence
        calls=out("populations/{prefix}_sex_id/{prefix}_scaffold_calls.tsv"),
        data=out("populations/{prefix}_sex_id/{prefix}_depth_heterozygosity.csv"),
        plot=out("populations/{prefix}_sex_id/{prefix}_chromosome_depth_heterozygosity.png"),
        exclude_scaffolds=out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt"),
    params:
        window=config["sex_id_window"],
        step=config["sex_id_step"],
        rule_min_size=config["sex_id_rule_min_size"],
        extra_exclude=config["sex_id_exclude"],
        keep=config["sex_id_keep"],
    resources:
        mem_mb=16000,
    conda:
        "../envs/sex_id.yml"
    log:
        logp("populations/{prefix}_sex_id_analyze.log"),
    script:
        "../scripts/identify_sex_linked_scaffolds.py"


rule scaffold_anomalies:
    """Scaffolds with odd depth/heterozygosity: haplotigs, SVs, contamination."""
    input:
        depth_het=out("populations/{prefix}_sex_id/{prefix}_depth_heterozygosity.csv"),
        vcf=out("prepare/{prefix}_pruned_ld.vcf.gz"),
    output:
        csv=out("populations/{prefix}_sex_id/{prefix}_scaffold_anomalies.csv"),
        txt=out("populations/{prefix}_sex_id/{prefix}_scaffold_anomalies_summary.txt"),
        png=out("populations/{prefix}_sex_id/{prefix}_anomaly_visualization.png"),
    conda:
        "../envs/sex_id.yml"
    log:
        logp("populations/{prefix}_scaffold_anomalies.log"),
    script:
        "../scripts/analyze_scaffold_anomalies.py"


# ------------------------------------------------------- population assignment
checkpoint assign_populations:
    """Split samples into populations. A checkpoint because K — and therefore how
    many downstream jobs exist — is only known once this has run."""
    input:
        fam=out("structure/admixture/{prefix}.fam"),
        cv_errors=out("structure/admixture/logs/{prefix}_best_K.txt"),
        coords=out("prepare/{prefix}.coords.txt"),
        composite=out("structure/{prefix}_structure_composite.pdf"),
        admixture_q=expand(
            out("structure/admixture/best_seed/{prefix}.{K}.Q"),
            K=ADMIX_K, allow_missing=True,
        ),
        tess=[out("structure/{prefix}_TESS_qmatrix.csv")] if RUN_TESS else [],
    output:
        plot=out("populations/{prefix}_ancestry_final.pdf"),
        outdir=directory(out("populations/{prefix}_populations")),
        all_qvalues=out("populations/{prefix}_all_qvalues.txt"),
        pops_done=out("populations/{prefix}_populations-done.txt"),
    params:
        admixture_path=out("structure/admixture/best_seed/{prefix}"),
        method=ASSIGN_METHOD,
        manual_k_assignment=MANUAL_K_PARAM,
        threshold_qval=config["qvalue_threshold"],
        exclude_samples=config.get("exclude_samples") or [],
        use_coords=USE_COORDS,
        run_tess=RUN_TESS,
        project_id=PREFIX,
    conda:
        "../envs/pie.yml"
    log:
        logp("populations/{prefix}_assign_populations.log"),
    script:
        "../scripts/ancestry_assignment.R"


rule population_vcf:
    """Per-population VCF, with sex-linked/anomalous scaffolds removed.

    Allele counts are recomputed for the population (+fill-tags) so that AF,
    F_MISSING and everything derived from them refer to this population and not
    to the whole dataset. Sites are kept even when invariant within the
    population, so per-population pi denominators stay comparable.
    """
    input:
        pops=out("populations/{prefix}_populations/population_{k}.txt"),
        pops_done=out("populations/{prefix}_populations-done.txt"),
        vcf=out("prepare/{prefix}_analysis_set.vcf.gz"),
        idx=out("prepare/{prefix}_analysis_set.vcf.gz.csi"),
        exclude_scaffolds=out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt"),
    output:
        vcf=out("populations/{prefix}_population_vcf/population_{k}.vcf.gz"),
        csi=out("populations/{prefix}_population_vcf/population_{k}.vcf.gz.csi"),
    threads: 4
    conda:
        "../envs/bcftools.yml"
    log:
        logp("populations/{prefix}_population_vcf_{k}.log"),
    shell:
        """
        exec > {log} 2>&1
        # bcftools rejects an empty -T file ("Failed to read the targets"), so
        # the exclusion is only applied when there is something to exclude.
        EXC_ARG=""
        EXC_BED={output.vcf}.exclude.bed
        if [ -s {input.exclude_scaffolds} ]; then
            awk 'BEGIN{{OFS="\\t"}} {{print $1, 1, 2147483647}}' {input.exclude_scaffolds} > $EXC_BED
            EXC_ARG="-T ^$EXC_BED"
        else
            echo "no scaffolds to exclude"
        fi

        bcftools view --threads {threads} -S {input.pops} $EXC_ARG -Ou -a -U {input.vcf} \
          | bcftools +fill-tags --threads {threads} -O z -o {output.vcf}
        rm -f $EXC_BED
        # CSI only: a .tbi cannot index contigs over 512 Mb (58-Pseudacris).
        bcftools index {output.vcf}
        """


# ------------------------------------------------------------------ ROH / F_ROH
rule population_roh:
    input:
        vcf=out("populations/{prefix}_population_vcf/population_{k}.vcf.gz"),
    output:
        # Only the RG (region) lines are used; -Orz writes just those, without
        # bcftools roh's per-site, per-sample lines (2-Doryteuthis, 98 samples
        # x 280M sites: >100 GB and >13 h with them; on a 2.7M-site contig
        # -Orz took 81 s vs 396 s with identical RG lines).
        roh=temp(out("populations/{prefix}_roh/population_{k}.roh.gz")),
        rg=out("populations/{prefix}_roh/population_{k}.rg.roh"),
    threads: 10
    resources:
        mem_mb=32000,
    conda:
        "../envs/bcftools.yml"
    log:
        logp("populations/{prefix}_roh_{k}.log"),
    shell:
        """
        exec 2> {log}
        bcftools roh --threads {threads} -G30 -Orz -o {output.roh} {input.vcf}
        zgrep -E '^RG[[:space:]]|# RG' {output.roh} > {output.rg}
        """


rule merge_roh:
    """Merge ROH separated by small gaps into single runs."""
    input:
        rg=out("populations/{prefix}_roh/population_{k}.rg.roh"),
    output:
        merged_rg_ind=out("populations/{prefix}_merged_roh/population_{k}.rg.roh"),
    params:
        gap_size=config["roh_merge_gap"],
    conda:
        "../envs/bcftools.yml"
    log:
        logp("populations/{prefix}_merge_roh_{k}.log"),
    script:
        "../scripts/merge_roh.py"


rule froh:
    """F_ROH per sample: fraction of the genome in ROH >= roh_min_size, and
    (_1pct.froh) in ROH covering >= roh_min_frac of their scaffold.

    The genome is the scaffolds > min_scaffold_size that are not excluded: ROH
    cannot be called on excluded scaffolds (they are not in the population VCF),
    so their length is left out of the denominator too.
    """
    input:
        rg=out("populations/{prefix}_merged_roh/population_{k}.rg.roh"),
        fai=CONTIGS,
        exclude_scaffolds=out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt"),
    output:
        froh=out("populations/{prefix}_merged_roh/{prefix}_population_{k}.froh"),
        top=out("populations/{prefix}_merged_roh/{prefix}_population_{k}_top.froh"),
        # Same, counting ROH covering >= roh_min_frac of their scaffold.
        froh_frac=out("populations/{prefix}_merged_roh/{prefix}_population_{k}_1pct.froh"),
    params:
        population="{k}",
        min_roh_size=config["roh_min_size"],
        min_roh_frac=config["roh_min_frac"],
        min_scaffold_size=config["min_scaffold_size"],
    conda:
        "../envs/froh.yml"
    log:
        logp("populations/{prefix}_froh_{k}.log"),
    script:
        "../scripts/froh.py"


# ------------------------------------------------------------------- pi / het
rule population_pi:
    """Genome-wide pi and per-sample heterozygosity for one population.

    The denominator is the callable-sites BED when one is configured, otherwise
    the whole assembly length, in both cases without the excluded scaffolds, whose
    variants are not in the population VCF — see scripts/calculate_pi.py.
    """
    input:
        pop_vcf=out("populations/{prefix}_population_vcf/population_{k}.vcf.gz"),
        roh_file=out("populations/{prefix}_merged_roh/population_{k}.rg.roh"),
        fai=CONTIGS,
        exclude_scaffolds=out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt"),
        bed=[CALLABLE_BED] if USE_BED else [],
    output:
        out_pi=out("populations/{prefix}_population_pi/pi_population_{k}.csv"),
        out_het=out("populations/{prefix}_population_pi/het_population_{k}.csv"),
    params:
        project_id=PREFIX,
        k_value="{k}",
    conda:
        "../envs/pi.yml"
    log:
        logp("populations/{prefix}_pi_{k}.log"),
    script:
        "../scripts/calculate_pi.py"


rule individual_pi:
    input:
        vcf=out("populations/{prefix}_population_vcf/population_{k}.vcf.gz"),
    output:
        # vcftools cannot write compressed output; the table is gzipped after.
        pi=out("populations/{prefix}_individual_pi/{k}/{sample}.windowed.pi.gz"),
    params:
        out=lambda wc, output: output.pi[: -len(".windowed.pi.gz")],
        window=config["pi_window"],
    conda:
        "../envs/vcftools.yml"
    log:
        logp("populations/{prefix}_individual_pi_{k}_{sample}.log"),
    shell:
        """
        vcftools --gzvcf {input.vcf} --indv {wildcards.sample} \
          --window-pi {params.window} --out {params.out} > {log} 2>&1
        gzip -f {params.out}.windowed.pi
        """


def individual_pi_files(wc):
    """Per-individual pi files for one population, from the checkpoint output."""
    pops = population_files()
    return expand(
        out("populations/{prefix}_individual_pi/{k}/{sample}.windowed.pi.gz"),
        k=wc.k, sample=pops.get(wc.k, []),
    )


rule roh_figures:
    input:
        pi_files=individual_pi_files,
        roh=out("populations/{prefix}_merged_roh/population_{k}.rg.roh"),
        top=out("populations/{prefix}_merged_roh/{prefix}_population_{k}_top.froh"),
    output:
        roh1=out("populations/{prefix}_merged_roh/pop_{k}_all_indiv_roh.pdf"),
        roh2=out("populations/{prefix}_merged_roh/pop_{k}_pi_roh_top.pdf"),
        roh3=out("populations/{prefix}_merged_roh/pop_{k}_roh_pi_ridges.pdf"),
    params:
        population="{k}",
    conda:
        "../envs/rplots.yml"
    log:
        logp("populations/{prefix}_roh_figures_{k}.log"),
    script:
        "../scripts/roh_viz_updated.R"


# -------------------------------------------------------------------- gather
def per_population_targets(wc):
    """Everything that has to exist, once per surviving population."""
    ks = active_k()
    return [
        *expand(out("populations/{prefix}_population_pi/pi_population_{k}.csv"), k=ks),
        *expand(out("populations/{prefix}_population_pi/het_population_{k}.csv"), k=ks),
        *expand(out("populations/{prefix}_merged_roh/{prefix}_population_{k}_top.froh"), k=ks),
        *expand(out("populations/{prefix}_merged_roh/pop_{k}_all_indiv_roh.pdf"), k=ks),
    ]


rule population_stats:
    """Sentinel: every per-population statistic has been computed."""
    input:
        per_population_targets,
    output:
        touch(out("populations/{prefix}_stats.done")),
