"""Ancestry inference: PCA, ADMIXTURE, and optionally TESS.

ADMIXTURE runs for K = 1..admixture_kmax across admixture_seeds replicates; the
first seed is the one plotted and used for population assignment. TESS is the
spatially-explicit counterpart and only exists when coordinates are available.

The composite PDF is the figure to look at before choosing manual_k_assignment.

Outputs, under {outdir}/structure/:
  {prefix}.eigenvec / .eigenval                PCA
  admixture/Q_files/{prefix}.{K}.{seed}.Q      ancestry proportions
  {prefix}_best_K.txt                          CV scores, best K, Evanno dK
  {prefix}_TESS_qmatrix.csv                     TESS ancestry (coords only)
  {prefix}_TESS_pieplot.pdf                     TESS pies on a map (coords only)
  {prefix}_structure_composite.pdf              everything, in one figure
"""


rule pca:
    input:
        vcf=out("prepare/{prefix}_subsampled.vcf.gz"),
    output:
        eigenvec=out("structure/{prefix}.eigenvec"),
        eigenval=out("structure/{prefix}.eigenval"),
    params:
        out=lambda wc, output: output.eigenvec[: -len(".eigenvec")],
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000,
    conda:
        "../envs/plink.yml"
    log:
        logp("structure/{prefix}_pca.log"),
    shell:
        """
        plink2 --vcf {input.vcf} --pca 10 --out {params.out} \
          --allow-extra-chr --autosome-num 95 --const-fid --bad-freqs > {log} 2>&1
        """


rule admixture_plink_input:
    """plink binary set for ADMIXTURE.

    The .bim is written aside as .bim.orig because ADMIXTURE needs integer
    chromosome names, supplied by the next rule. Keeping the two files separate
    means neither rule mutates its own input.
    """
    input:
        vcf=out("prepare/{prefix}_subsampled.vcf.gz"),
    output:
        bed=out("structure/admixture/{prefix}.bed"),
        bim_orig=out("structure/admixture/{prefix}.bim.orig"),
        fam=out("structure/admixture/{prefix}.fam"),
    params:
        out=lambda wc, output: output.bed[: -len(".bed")],
    conda:
        "../envs/admixture.yml"
    log:
        logp("structure/{prefix}_admixture_input.log"),
    shell:
        """
        plink --vcf {input.vcf} --make-bed --out {params.out} \
          --allow-extra-chr --autosome-num 95 --const-fid > {log} 2>&1
        mv {params.out}.bim {output.bim_orig}
        """


rule admixture_integer_contigs:
    """Renumber contigs 1..n in the .bim; ADMIXTURE rejects names."""
    input:
        bim=out("structure/admixture/{prefix}.bim.orig"),
        fai=CONTIGS,
    output:
        bim=out("structure/admixture/{prefix}.bim"),
    log:
        logp("structure/{prefix}_admixture_bim.log"),
    script:
        "../scripts/contigs4admixture.py"


rule admixture:
    input:
        bed=out("structure/admixture/{prefix}.bed"),
        bim=out("structure/admixture/{prefix}.bim"),
        fam=out("structure/admixture/{prefix}.fam"),
    output:
        Q=out("structure/admixture/Q_files/{prefix}.{K}.{seed}.Q"),
        # Allele frequencies; nothing downstream reads them.
        P=temp(out("structure/admixture/P_files/{prefix}.{K}.{seed}.P")),
    log:
        # Parsed by admixture_cloud.R for CV scores, so the _K{K}_seed{seed}.log
        # shape matters.
        out("structure/admixture/logs/{prefix}_admix_K{K}_seed{seed}.log"),
    threads: 4
    conda:
        "../envs/admixture.yml"
    shell:
        # ADMIXTURE always writes <bed basename>.<K>.{P,Q} into the working
        # directory, so each replicate runs in its own temp dir; otherwise
        # concurrent K/seed jobs would collide on those names.
        """
        BED=$(readlink -f {input.bed})
        LOG=$(realpath -m {log})
        QOUT=$(realpath -m {output.Q})
        POUT=$(realpath -m {output.P})
        mkdir -p $(dirname $LOG) $(dirname $QOUT) $(dirname $POUT)

        TMP=$(mktemp -d)
        trap 'rm -rf $TMP' EXIT
        cd $TMP
        admixture --cv=5 -j{threads} --seed={wildcards.seed} $BED {wildcards.K} > $LOG 2>&1

        BASE=$(basename $BED .bed)
        mv $BASE.{wildcards.K}.P $POUT
        mv $BASE.{wildcards.K}.Q $QOUT
        """


rule admixture_best_k:
    """Cross-validation scores across all K and seeds, plus Evanno's delta K."""
    input:
        logs=expand(
            out("structure/admixture/logs/{prefix}_admix_K{K}_seed{seed}.log"),
            K=ADMIX_K,
            seed=ADMIX_SEEDS,
        ),
    output:
        best_k=out("structure/admixture/logs/{prefix}_best_K.txt"),
        cv_plot=out("structure/admixture/logs/{prefix}_admixture_cv.png"),
        deltak_plot=out("structure/admixture/logs/{prefix}_admixture_deltak.png"),
    params:
        species=PREFIX,
        # admixture_cloud.R globs *.log in this directory; it must end in "/".
        output_path=out("structure/admixture/logs/"),
    conda:
        "../envs/pie.yml"
    log:
        logp("structure/{prefix}_best_k.log"),
    script:
        "../scripts/admixture_cloud.R"


if RUN_TESS:

    rule tess_env:
        """Install tess3r, the one package TESS needs that is not on conda.

        It is a small Rcpp package, so this is a quick compile. With a shared
        --conda-prefix it happens once per environment rather than per project;
        the check makes later projects a no-op.
        """
        output:
            done=out("structure/{prefix}_tess_env.txt"),
        params:
            ref=config["tess3r_ref"],
        conda:
            "../envs/tess.yml"
        log:
            logp("structure/{prefix}_tess_env.log"),
        script:
            "../scripts/setup_tess_env.R"

    rule tess:
        """Spatially-explicit ancestry (TESS3) over tess_kvals, with the
        cross-entropy K selection that algatr's tess_ktest(K_selection="auto")
        used."""
        input:
            vcf=out("prepare/{prefix}_subsampled.vcf.gz"),
            coords=out("prepare/{prefix}.coords.txt"),
            env=out("structure/{prefix}_tess_env.txt"),
        output:
            qmatrix=out("structure/{prefix}_TESS_qmatrix.csv"),
            xval=out("structure/{prefix}_TESS_xval.csv"),
            xval_plot=out("structure/{prefix}_TESS_xval.png"),
        params:
            species=PREFIX,
            kvals=config["tess_kvals"],
            seed=config["tess_seed"],
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000,
        conda:
            "../envs/tess.yml"
        log:
            logp("structure/{prefix}_tess.log"),
        benchmark:
            logp("structure/{prefix}_tess.benchmark.tsv")
        script:
            "../scripts/run_tess.R"

    rule tess_pieplot:
        input:
            tess=out("structure/{prefix}_TESS_qmatrix.csv"),
            coords=out("prepare/{prefix}.coords.txt"),
        output:
            plot=out("structure/{prefix}_TESS_pieplot.pdf"),
        params:
            species=PREFIX,
        conda:
            "../envs/pie.yml"
        log:
            logp("structure/{prefix}_tess_pieplot.log"),
        script:
            "../scripts/pie_plot.R"


rule structure_composite:
    """One figure summarising every ancestry result, for choosing K.

    Panels present depend on the mode: ADMIXTURE barplot, CV/deltaK and PCA are
    always drawn; maps need coordinates; the TESS panels and the
    ADMIXTURE-vs-TESS agreement panels need TESS.
    """
    input:
        cv_errors=out("structure/admixture/logs/{prefix}_best_K.txt"),
        cv_plot=out("structure/admixture/logs/{prefix}_admixture_cv.png"),
        deltak_plot=out("structure/admixture/logs/{prefix}_admixture_deltak.png"),
        eigenvec=out("structure/{prefix}.eigenvec"),
        eigenval=out("structure/{prefix}.eigenval"),
        fam=out("structure/admixture/{prefix}.fam"),
        coords=out("prepare/{prefix}.coords.txt"),
        admixture_logs=expand(
            out("structure/admixture/logs/{prefix}_admix_K{K}_seed{seed}.log"),
            K=ADMIX_K, seed=[PLOT_SEED],
        ),
        admixture_q=expand(
            out("structure/admixture/Q_files/{prefix}.{K}.{seed}.Q"),
            K=ADMIX_K, seed=[PLOT_SEED],
        ),
        tess=[out("structure/{prefix}_TESS_qmatrix.csv")] if RUN_TESS else [],
        xval=[out("structure/{prefix}_TESS_xval.csv")] if RUN_TESS else [],
    output:
        pdf=out("structure/{prefix}_structure_composite.pdf"),
        merged=out("structure/{prefix}_admixture_tess_merged.csv"),
    params:
        admixture_path=out("structure/admixture/Q_files/{prefix}"),
        plot_seed=PLOT_SEED,
        use_coords=USE_COORDS,
        run_tess=RUN_TESS,
        kmax=config["admixture_kmax"],
    conda:
        "../envs/pie.yml"
    log:
        logp("structure/{prefix}_structure_composite.log"),
    script:
        "../scripts/structure_composite.R"
