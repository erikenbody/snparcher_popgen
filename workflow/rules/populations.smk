"""Per-population statistics, downstream of the population assignment checkpoint.

Order of operations:
  1. sex-linked and anomalous scaffolds are identified from depth and
     heterozygosity, and excluded from every per-population VCF
  2. samples are assigned to populations (checkpoint: the number of populations
     is not known until the ancestry coefficients have been inspected)
  3. per population: ROH -> F_ROH, pi and heterozygosity, per-individual pi and
     ROH figures, and optionally GONE

Outputs, under results/{refGenome}/populations/:
  {prefix}_sex_id/                          scaffold QC: candidates, anomalies, plots
  {prefix}_ancestry_final.pdf               the assignment actually used
  {prefix}_populations/population_{k}.txt   sample lists, one per population
  {prefix}_population_vcf/                  per-population VCFs
  {prefix}_merged_roh/                      merged ROH, F_ROH, ROH figures
  {prefix}_population_pi/                   pi and heterozygosity per population
  {prefix}_GONE/                            recent Ne trajectories (run_gone)
  {prefix}_stats.done                       sentinel: all per-population work done
"""


# ---------------------------------------------------------------- scaffold QC
rule sex_id_extract:
    """Pull CHROM/SAMPLE/DP/GT out of the pruned VCF for the scaffold QC below."""
    input:
        vcf="results/{refGenome}/prepare/{prefix}_pruned_ld.vcf.gz",
        idx="results/{refGenome}/prepare/{prefix}_pruned_ld.vcf.gz.csi",
    output:
        raw_data="results/{refGenome}/populations/{prefix}_sex_id/{prefix}_raw_data.tsv.gz",
    params:
        query_format=r"[%CHROM\t%SAMPLE\t%DP\t%GT\n]",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/{refGenome}/populations/{prefix}_sex_id_extract.log",
    shell:
        "bcftools query -f '{params.query_format}' {input.vcf} | gzip > {output.raw_data} 2> {log}"


rule sex_id_analyze:
    """Flag candidate sex-linked scaffolds from normalised depth + heterozygosity.

    The exclude list this produces is applied to every per-population VCF and to
    clam, so that sex-linked scaffolds do not distort diversity estimates.
    """
    input:
        raw_data="results/{refGenome}/populations/{prefix}_sex_id/{prefix}_raw_data.tsv.gz",
        vcf_path_for_lengths="results/{refGenome}/prepare/{prefix}_pruned_ld.vcf.gz",
    output:
        report="results/{refGenome}/populations/{prefix}_sex_id/{prefix}_sex_chromosome_candidates.csv",
        data="results/{refGenome}/populations/{prefix}_sex_id/{prefix}_depth_heterozygosity.csv",
        plot="results/{refGenome}/populations/{prefix}_sex_id/{prefix}_chromosome_depth_heterozygosity.png",
        exclude_scaffolds="results/{refGenome}/populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt",
    conda:
        "../envs/sex_id.yml"
    log:
        "logs/{refGenome}/populations/{prefix}_sex_id_analyze.log",
    script:
        "../scripts/identify_sex_linked_scaffolds.py"


rule scaffold_anomalies:
    """Scaffolds with odd depth/heterozygosity: haplotigs, SVs, contamination."""
    input:
        depth_het="results/{refGenome}/populations/{prefix}_sex_id/{prefix}_depth_heterozygosity.csv",
        vcf="results/{refGenome}/prepare/{prefix}_pruned_ld.vcf.gz",
    output:
        csv="results/{refGenome}/populations/{prefix}_sex_id/{prefix}_scaffold_anomalies.csv",
        txt="results/{refGenome}/populations/{prefix}_sex_id/{prefix}_scaffold_anomalies_summary.txt",
        png="results/{refGenome}/populations/{prefix}_sex_id/{prefix}_anomaly_visualization.png",
    conda:
        "../envs/sex_id.yml"
    log:
        "logs/{refGenome}/populations/{prefix}_scaffold_anomalies.log",
    script:
        "../scripts/analyze_scaffold_anomalies.py"


# ------------------------------------------------------- population assignment
checkpoint assign_populations:
    """Split samples into populations. A checkpoint because K — and therefore how
    many downstream jobs exist — is only known once this has run."""
    input:
        unpack(
            lambda wc: {
                "fam": f"results/{wc.refGenome}/structure/admixture/{wc.prefix}.fam",
                "cv_errors": f"results/{wc.refGenome}/structure/admixture/logs/{wc.prefix}_best_K.txt",
                "coords": f"results/{wc.refGenome}/prepare/{wc.prefix}.coords.txt",
                "composite": f"results/{wc.refGenome}/structure/{wc.prefix}_structure_composite.pdf",
                "admixture_q": expand(
                    "results/{refGenome}/structure/admixture/Q_files/{prefix}.{K}.{seed}.Q",
                    refGenome=wc.refGenome, prefix=wc.prefix, K=ADMIX_K, seed=[PLOT_SEED],
                ),
                **(
                    {"tess": f"results/{wc.refGenome}/structure/{wc.prefix}_TESS_qmatrix.csv"}
                    if RUN_TESS
                    else {}
                ),
            }
        ),
    output:
        plot="results/{refGenome}/populations/{prefix}_ancestry_final.pdf",
        outdir=directory("results/{refGenome}/populations/{prefix}_populations"),
        all_qvalues="results/{refGenome}/populations/{prefix}_all_qvalues.txt",
        pops_done="results/{refGenome}/populations/{prefix}_populations-done.txt",
    params:
        admixture_path=lambda wc: f"results/{wc.refGenome}/structure/admixture/Q_files/{wc.prefix}",
        plot_seed=PLOT_SEED,
        method=ASSIGN_METHOD,
        manual_k_assignment=MANUAL_K_PARAM,
        threshold_qval=config["qvalue_threshold"],
        use_coords=USE_COORDS,
        run_tess=RUN_TESS,
        project_id=lambda wc: wc.prefix,
    conda:
        "../envs/pie.yml"
    log:
        "logs/{refGenome}/populations/{prefix}_assign_populations.log",
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
        pops="results/{refGenome}/populations/{prefix}_populations/population_{k}.txt",
        pops_done="results/{refGenome}/populations/{prefix}_populations-done.txt",
        vcf="results/{refGenome}/prepare/{prefix}_analysis_set.vcf.gz",
        idx="results/{refGenome}/prepare/{prefix}_analysis_set.vcf.gz.csi",
        exclude_scaffolds="results/{refGenome}/populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt",
    output:
        vcf="results/{refGenome}/populations/{prefix}_population_vcf/population_{k}.vcf.gz",
        csi="results/{refGenome}/populations/{prefix}_population_vcf/population_{k}.vcf.gz.csi",
        tbi="results/{refGenome}/populations/{prefix}_population_vcf/population_{k}.vcf.gz.tbi",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/{refGenome}/populations/{prefix}_population_vcf_{k}.log",
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

        bcftools view -S {input.pops} $EXC_ARG -Ou -a -U {input.vcf} \
          | bcftools +fill-tags -O z -o {output.vcf}
        rm -f $EXC_BED
        bcftools index {output.vcf}
        bcftools index -t {output.vcf}
        """


# ------------------------------------------------------------------ ROH / F_ROH
rule population_roh:
    input:
        vcf="results/{refGenome}/populations/{prefix}_population_vcf/population_{k}.vcf.gz",
    output:
        roh="results/{refGenome}/populations/{prefix}_roh/population_{k}.roh.gz",
        rg="results/{refGenome}/populations/{prefix}_roh/population_{k}.rg.roh",
    threads: 10
    resources:
        mem_mb=32000,
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/{refGenome}/populations/{prefix}_roh_{k}.log",
    shell:
        """
        exec 2> {log}
        bcftools roh --threads {threads} -G30 -Oz -o {output.roh} {input.vcf}
        zgrep -E '^RG[[:space:]]|# RG' {output.roh} > {output.rg}
        """


rule merge_roh:
    """Merge ROH separated by small gaps into single runs."""
    input:
        rg="results/{refGenome}/populations/{prefix}_roh/population_{k}.rg.roh",
    output:
        merged_rg_ind="results/{refGenome}/populations/{prefix}_merged_roh/population_{k}.rg.roh",
    params:
        gap_size=config["roh_merge_gap"],
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/{refGenome}/populations/{prefix}_merge_roh_{k}.log",
    script:
        "../scripts/merge_roh.py"


rule froh:
    """F_ROH per sample: fraction of the genome in ROH >= roh_min_size."""
    input:
        rg="results/{refGenome}/populations/{prefix}_merged_roh/population_{k}.rg.roh",
        fai="results/{refGenome}/data/genome/{refGenome}.fna.fai",
    output:
        froh="results/{refGenome}/populations/{prefix}_merged_roh/{prefix}_population_{k}.froh",
        top="results/{refGenome}/populations/{prefix}_merged_roh/{prefix}_population_{k}_top.froh",
    params:
        population="{k}",
        min_roh_size=config["roh_min_size"],
        min_scaffold_size=config["min_scaffold_size"],
    conda:
        "../envs/froh.yml"
    log:
        "logs/{refGenome}/populations/{prefix}_froh_{k}.log",
    script:
        "../scripts/froh.py"


# ------------------------------------------------------------------- pi / het
rule population_pi:
    """Genome-wide pi and per-sample heterozygosity for one population.

    The denominator is the callable-sites BED when one is configured, otherwise
    the whole assembly length — see scripts/calculate_pi.py.
    """
    input:
        unpack(
            lambda wc: {
                "pop_vcf": f"results/{wc.refGenome}/populations/{wc.prefix}_population_vcf/population_{wc.k}.vcf.gz",
                "roh_file": f"results/{wc.refGenome}/populations/{wc.prefix}_merged_roh/population_{wc.k}.rg.roh",
                "fai": f"results/{wc.refGenome}/data/genome/{wc.refGenome}.fna.fai",
                **({"bed": CALLABLE_BED} if USE_BED else {}),
            }
        ),
    output:
        out_pi="results/{refGenome}/populations/{prefix}_population_pi/pi_population_{k}.csv",
        out_het="results/{refGenome}/populations/{prefix}_population_pi/het_population_{k}.csv",
    params:
        project_id="{prefix}",
        k_value="{k}",
    conda:
        "../envs/pi.yml"
    log:
        "logs/{refGenome}/populations/{prefix}_pi_{k}.log",
    script:
        "../scripts/calculate_pi.py"


rule individual_pi:
    input:
        vcf="results/{refGenome}/populations/{prefix}_population_vcf/population_{k}.vcf.gz",
    output:
        pi="results/{refGenome}/populations/{prefix}_individual_pi/{k}/{sample}.windowed.pi",
    params:
        out=lambda wc, output: output.pi[: -len(".windowed.pi")],
        window=config["pi_window"],
    conda:
        "../envs/vcftools.yml"
    log:
        "logs/{refGenome}/populations/{prefix}_individual_pi_{k}_{sample}.log",
    shell:
        """
        vcftools --gzvcf {input.vcf} --indv {wildcards.sample} \
          --window-pi {params.window} --out {params.out} > {log} 2>&1
        """


def individual_pi_files(wc):
    """Per-individual pi files for one population, from the checkpoint output."""
    pops = population_files(wc.refGenome, wc.prefix)
    return expand(
        "results/{refGenome}/populations/{prefix}_individual_pi/{k}/{sample}.windowed.pi",
        refGenome=wc.refGenome, prefix=wc.prefix, k=wc.k, sample=pops.get(wc.k, []),
    )


rule roh_figures:
    input:
        pi_files=individual_pi_files,
        roh="results/{refGenome}/populations/{prefix}_merged_roh/population_{k}.rg.roh",
        top="results/{refGenome}/populations/{prefix}_merged_roh/{prefix}_population_{k}_top.froh",
    output:
        roh1="results/{refGenome}/populations/{prefix}_merged_roh/pop_{k}_all_indiv_roh.pdf",
        roh2="results/{refGenome}/populations/{prefix}_merged_roh/pop_{k}_pi_roh_top.pdf",
        roh3="results/{refGenome}/populations/{prefix}_merged_roh/pop_{k}_roh_pi_ridges.pdf",
    params:
        population="{k}",
    conda:
        "../envs/rplots.yml"
    log:
        "logs/{refGenome}/populations/{prefix}_roh_figures_{k}.log",
    script:
        "../scripts/roh_viz_updated.R"


# ----------------------------------------------------------------------- GONE
if RUN_GONE:

    if not config["gone_path"]:

        rule gone_setup:
            """Fetch GONE. It ships prebuilt Linux binaries under Linux/PROGRAMMES."""
            output:
                script=f"{gone_dir()}/script_GONE.sh",
                params_file=f"{gone_dir()}/INPUT_PARAMETERS_FILE",
            params:
                dest=lambda wc, output: str(Path(output.script).parent.parent),
            log:
                "logs/gone_setup.log",
            shell:
                """
                exec > {log} 2>&1
                rm -rf {params.dest}
                git clone --depth 1 https://github.com/esrud/GONE.git {params.dest}
                chmod -R +x {params.dest}/Linux
                """

    rule gone_contigs:
        """The gone_n_chroms largest scaffolds; GONE expects a handful of chromosomes."""
        input:
            fai="results/{refGenome}/data/genome/{refGenome}.fna.fai",
        output:
            chrom_names="results/{refGenome}/populations/{prefix}_GONE/{prefix}_GONE_contigs.txt",
        params:
            n=config["gone_n_chroms"],
        log:
            "logs/{refGenome}/populations/{prefix}_gone_contigs.log",
        shell:
            """
            awk '{{print $2, $1}}' {input.fai} 2> {log} \
              | sort -k1,1nr | head -n {params.n} | awk '{{print $2}}' > {output.chrom_names}
            """

    rule gone_subset_vcf:
        input:
            vcf="results/{refGenome}/populations/{prefix}_population_vcf/population_{k}.vcf.gz",
            chrom_names="results/{refGenome}/populations/{prefix}_GONE/{prefix}_GONE_contigs.txt",
        output:
            vcf=temp("results/{refGenome}/populations/{prefix}_GONE/subsets/population_{k}.vcf.gz"),
            csi=temp("results/{refGenome}/populations/{prefix}_GONE/subsets/population_{k}.vcf.gz.csi"),
        conda:
            "../envs/bcftools.yml"
        log:
            "logs/{refGenome}/populations/{prefix}_gone_subset_{k}.log",
        shell:
            """
            exec 2> {log}
            bcftools view -t $(paste -sd, {input.chrom_names}) -Oz -o {output.vcf} {input.vcf}
            bcftools index {output.vcf}
            """

    rule gone_clean_snps:
        """Drop multi-base and fixed sites.

        GATK writes SNPs that overlap an indel with multi-base REF/ALT (GTG>GCG).
        plink turns those into malformed ped rows and GONE then fails, so they
        have to go. Sites fixed within the population carry no LD information.
        """
        input:
            vcf="results/{refGenome}/populations/{prefix}_GONE/subsets/population_{k}.vcf.gz",
            csi="results/{refGenome}/populations/{prefix}_GONE/subsets/population_{k}.vcf.gz.csi",
        output:
            vcf="results/{refGenome}/populations/{prefix}_GONE/vcfs/population_{k}.vcf.gz",
            csi="results/{refGenome}/populations/{prefix}_GONE/vcfs/population_{k}.vcf.gz.csi",
            positions=temp("results/{refGenome}/populations/{prefix}_GONE/vcfs/population_{k}_positions.txt.gz"),
        conda:
            "../envs/bcftools.yml"
        log:
            "logs/{refGenome}/populations/{prefix}_gone_clean_{k}.log",
        shell:
            """
            exec 2> {log}
            bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {input.vcf} \
              | awk 'length($3) == 1 && length($4) == 1 {{print $1"\\t"$2}}' \
              | bgzip -c > {output.positions}
            tabix --csi -s1 -b2 -e2 {output.positions}

            bcftools view -T {output.positions} {input.vcf} -O u -a -U \
              | bcftools +fill-tags -Ou \
              | bcftools view -e 'AF==0 | AF==1 | ref="N" | ALT="."' -O z -o {output.vcf}
            bcftools index {output.vcf}
            """

    rule gone_thin_snps:
        """Thin to ~30k SNPs per scaffold per replicate; GONE is LD-based and
        does not need every site, and each seed samples independently."""
        input:
            vcf="results/{refGenome}/populations/{prefix}_GONE/vcfs/population_{k}.vcf.gz",
            csi="results/{refGenome}/populations/{prefix}_GONE/vcfs/population_{k}.vcf.gz.csi",
            chrom_names="results/{refGenome}/populations/{prefix}_GONE/{prefix}_GONE_contigs.txt",
        output:
            vcf=temp("results/{refGenome}/populations/{prefix}_GONE/vcfs/population_{k}_seed_{seed}.vcf.gz"),
            csi=temp("results/{refGenome}/populations/{prefix}_GONE/vcfs/population_{k}_seed_{seed}.vcf.gz.csi"),
        conda:
            "../envs/bcftools.yml"
        log:
            "logs/{refGenome}/populations/{prefix}_gone_thin_{k}_{seed}.log",
        shell:
            """
            exec 2> {log}
            TMPDIR=$(mktemp -d)
            trap "rm -rf $TMPDIR" EXIT

            for SCAFFOLD in $(cat {input.chrom_names}); do
                SCAF_VCF=$TMPDIR/${{SCAFFOLD}}.vcf.gz
                SCAF_PRUNED=$TMPDIR/${{SCAFFOLD}}_pruned.vcf.gz
                bcftools view -r $SCAFFOLD {input.vcf} -O z -o $SCAF_VCF

                ALLSITES=$(bcftools query -f '%CHROM\\t%POS\\n' $SCAF_VCF | wc -l)
                WINDOW=$(( ALLSITES / 30000 ))
                if [ "$WINDOW" -gt 1 ]; then
                    bcftools +prune --random-seed {wildcards.seed} -w $WINDOW -n 1 -N rand \
                      -O z -o $SCAF_PRUNED $SCAF_VCF
                else
                    cp $SCAF_VCF $SCAF_PRUNED
                fi
                bcftools index $SCAF_PRUNED
            done

            bcftools concat -O z -o {output.vcf} $TMPDIR/*_pruned.vcf.gz
            bcftools index {output.vcf}
            """

    rule gone_plink:
        input:
            vcf="results/{refGenome}/populations/{prefix}_GONE/vcfs/population_{k}_seed_{seed}.vcf.gz",
        output:
            bed=temp("results/{refGenome}/populations/{prefix}_GONE/plink/{prefix}_pop_{k}_seed_{seed}.bed"),
            map_orig=temp("results/{refGenome}/populations/{prefix}_GONE/plink/{prefix}_pop_{k}_seed_{seed}.orig.map"),
            ped="results/{refGenome}/populations/{prefix}_GONE/plink/{prefix}_pop_{k}_seed_{seed}.ped",
        params:
            out=lambda wc, output: output.bed[: -len(".bed")],
        conda:
            "../envs/plink1_9.yml"
        log:
            "logs/{refGenome}/populations/{prefix}_gone_plink_{k}_{seed}.log",
        shell:
            """
            exec > {log} 2>&1
            plink --vcf {input.vcf} --out {params.out} --double-id --make-bed --allow-extra-chr
            plink --bfile {params.out} --recode --out {params.out} --allow-extra-chr
            mv {params.out}.map {output.map_orig}
            """

    rule gone_rename_map:
        """GONE requires chromosomes numbered 1..n."""
        input:
            map_temp="results/{refGenome}/populations/{prefix}_GONE/plink/{prefix}_pop_{k}_seed_{seed}.orig.map",
        output:
            mapp="results/{refGenome}/populations/{prefix}_GONE/plink/{prefix}_pop_{k}_seed_{seed}.map",
        log:
            "logs/{refGenome}/populations/{prefix}_gone_rename_map_{k}_{seed}.log",
        script:
            "../scripts/rename_map.py"

    rule gone:
        """Run GONE in a shadow directory.

        GONE insists on running in its own working directory next to
        PROGRAMMES/, so inputs are symlinked in and the parameter file is
        patched for thread count rather than editing the installation.
        """
        input:
            mapp="results/{refGenome}/populations/{prefix}_GONE/plink/{prefix}_pop_{k}_seed_{seed}.map",
            ped="results/{refGenome}/populations/{prefix}_GONE/plink/{prefix}_pop_{k}_seed_{seed}.ped",
            gone_script=f"{gone_dir()}/script_GONE.sh",
            gone_params=f"{gone_dir()}/INPUT_PARAMETERS_FILE",
        output:
            d2="results/{refGenome}/populations/{prefix}_GONE/run_{k}_seed_{seed}/Output_d2_{prefix}_pop_{k}_seed_{seed}",
            ne="results/{refGenome}/populations/{prefix}_GONE/run_{k}_seed_{seed}/Output_Ne_{prefix}_pop_{k}_seed_{seed}",
            out="results/{refGenome}/populations/{prefix}_GONE/run_{k}_seed_{seed}/OUTPUT_{prefix}_pop_{k}_seed_{seed}",
        params:
            run_id="{prefix}_pop_{k}_seed_{seed}",
            gone_base=lambda wc, input: str(Path(input.gone_script).parent.resolve()),
        threads: 10
        shadow:
            "shallow"
        log:
            "logs/{refGenome}/populations/{prefix}_gone_{k}_{seed}.log",
        shell:
            """
            exec > {log} 2>&1
            ln -s $(readlink -f {input.mapp}) {params.run_id}.map
            ln -s $(readlink -f {input.ped}) {params.run_id}.ped

            sed 's/^threads=-99/threads={threads}/' {params.gone_base}/INPUT_PARAMETERS_FILE \
              > INPUT_PARAMETERS_FILE

            # Pass the thread count through to GONEparallel.sh, which otherwise
            # ignores it. A no-op on versions that already do.
            sed 's|\\./PROGRAMMES/GONEparallel\\.sh -hc \\$hc outfileLD \\$REPS|{params.gone_base}/PROGRAMMES/GONEparallel.sh -hc $hc outfileLD $REPS $threads|' \
                {params.gone_base}/script_GONE.sh > script_GONE.sh
            chmod +x script_GONE.sh
            ln -s {params.gone_base}/PROGRAMMES PROGRAMMES

            bash script_GONE.sh {params.run_id}

            mv Output_d2_{params.run_id} {output.d2}
            mv Output_Ne_{params.run_id} {output.ne}
            mv OUTPUT_{params.run_id} {output.out}
            """

    rule gone_summarize:
        """Combine GONE replicates for one population into a mean trajectory."""
        input:
            GONE_out_Ne=expand(
                "results/{{refGenome}}/populations/{{prefix}}_GONE/run_{{k}}_seed_{seed}/Output_Ne_{{prefix}}_pop_{{k}}_seed_{seed}",
                seed=GONE_SEEDS,
            ),
        output:
            GONE_plot="results/{refGenome}/populations/{prefix}_GONE/plots/{prefix}_population_{k}_all_seeds.png",
            GONE_avg_plot="results/{refGenome}/populations/{prefix}_GONE/plots/{prefix}_population_{k}_mean.png",
            GONE_summary_stats="results/{refGenome}/populations/{prefix}_GONE/plots/{prefix}_summary_stats_population_{k}.csv",
            GONE_Ne_change="results/{refGenome}/populations/{prefix}_GONE/plots/{prefix}_Ne_change_population_{k}.csv",
        params:
            project_id="{prefix}",
            population_id="{k}",
            seeds=GONE_SEEDS,
        conda:
            "../envs/rplots.yml"
        log:
            "logs/{refGenome}/populations/{prefix}_gone_summarize_{k}.log",
        script:
            "../scripts/plot_GONE_all_seeds.R"


# -------------------------------------------------------------------- gather
def per_population_targets(wc):
    """Everything that has to exist, once per surviving population."""
    ks = active_k(wc)
    targets = [
        *expand(
            "results/{refGenome}/populations/{prefix}_population_pi/pi_population_{k}.csv",
            refGenome=wc.refGenome, prefix=wc.prefix, k=ks,
        ),
        *expand(
            "results/{refGenome}/populations/{prefix}_population_pi/het_population_{k}.csv",
            refGenome=wc.refGenome, prefix=wc.prefix, k=ks,
        ),
        *expand(
            "results/{refGenome}/populations/{prefix}_merged_roh/{prefix}_population_{k}_top.froh",
            refGenome=wc.refGenome, prefix=wc.prefix, k=ks,
        ),
        *expand(
            "results/{refGenome}/populations/{prefix}_merged_roh/pop_{k}_all_indiv_roh.pdf",
            refGenome=wc.refGenome, prefix=wc.prefix, k=ks,
        ),
    ]
    if RUN_GONE:
        targets += expand(
            "results/{refGenome}/populations/{prefix}_GONE/plots/{prefix}_Ne_change_population_{k}.csv",
            refGenome=wc.refGenome, prefix=wc.prefix, k=ks,
        )
    return targets


rule population_stats:
    """Sentinel: every per-population statistic has been computed."""
    input:
        per_population_targets,
    output:
        touch("results/{refGenome}/populations/{prefix}_stats.done"),
