"""Ancestry inference: PCA, ADMIXTURE, and optionally TESS.

ADMIXTURE runs for K = 1..admixture_kmax across admixture_seeds replicates; the
first seed is the one plotted and used for population assignment. TESS is the
spatially-explicit counterpart and only exists when coordinates are available.

The composite PDF is the figure to look at before choosing manual_k_assignment.

Outputs:
  {prefix}.eigenvec / .eigenval                PCA
  admixture/Q_files/{prefix}.{K}.{seed}.Q      ancestry proportions
  {prefix}_best_K.txt                          CV scores, best K, Evanno dK
  {prefix}_TESS_qmatrix.csv                     TESS ancestry (coords only)
  {prefix}_TESS_pieplot.pdf                     TESS pies on a map (coords only)
  {prefix}_structure_composite.pdf              everything, in one figure
"""


rule pca:
    input:
        vcf="results/{refGenome}/prepare/{prefix}_subsampled.vcf.gz",
    output:
        eigenvec="results/{refGenome}/structure/{prefix}.eigenvec",
        eigenval="results/{refGenome}/structure/{prefix}.eigenval",
    params:
        out=lambda wc, output: output.eigenvec[: -len(".eigenvec")],
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000,
    conda:
        "../envs/plink.yml"
    log:
        "logs/{refGenome}/structure/{prefix}_pca.log",
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
        vcf="results/{refGenome}/prepare/{prefix}_subsampled.vcf.gz",
    output:
        bed="results/{refGenome}/structure/admixture/{prefix}.bed",
        # NOT "{prefix}.orig.bim": that filename is itself matched by the
        # "{prefix}.bim" output of admixture_integer_contigs, which makes the two
        # rules mutually ambiguous and sends DAG resolution into recursion.
        bim_orig="results/{refGenome}/structure/admixture/{prefix}.bim.orig",
        fam="results/{refGenome}/structure/admixture/{prefix}.fam",
    params:
        out=lambda wc, output: output.bed[: -len(".bed")],
    conda:
        "../envs/admixture.yml"
    log:
        "logs/{refGenome}/structure/{prefix}_admixture_input.log",
    shell:
        """
        plink --vcf {input.vcf} --make-bed --out {params.out} \
          --allow-extra-chr --autosome-num 95 --const-fid > {log} 2>&1
        mv {params.out}.bim {output.bim_orig}
        """


rule admixture_integer_contigs:
    """Renumber contigs 1..n in the .bim; ADMIXTURE rejects names."""
    input:
        bim="results/{refGenome}/structure/admixture/{prefix}.bim.orig",
        fai="results/{refGenome}/data/genome/{refGenome}.fna.fai",
    output:
        bim="results/{refGenome}/structure/admixture/{prefix}.bim",
    log:
        "logs/{refGenome}/structure/{prefix}_admixture_bim.log",
    script:
        "../scripts/contigs4admixture.py"


rule admixture:
    input:
        bed="results/{refGenome}/structure/admixture/{prefix}.bed",
        bim="results/{refGenome}/structure/admixture/{prefix}.bim",
        fam="results/{refGenome}/structure/admixture/{prefix}.fam",
    output:
        Q="results/{refGenome}/structure/admixture/Q_files/{prefix}.{K}.{seed}.Q",
        P="results/{refGenome}/structure/admixture/P_files/{prefix}.{K}.{seed}.P",
    log:
        # Parsed by admixture_cloud.R for CV scores, so the _K{K}_seed{seed}.log
        # shape matters.
        "results/{refGenome}/structure/admixture/logs/{prefix}_admix_K{K}_seed{seed}.log",
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
        admixture --cv=5 --seed={wildcards.seed} $BED {wildcards.K} > $LOG 2>&1

        BASE=$(basename $BED .bed)
        mv $BASE.{wildcards.K}.P $POUT
        mv $BASE.{wildcards.K}.Q $QOUT
        """


rule admixture_best_k:
    """Cross-validation scores across all K and seeds, plus Evanno's delta K."""
    input:
        logs=expand(
            "results/{{refGenome}}/structure/admixture/logs/{{prefix}}_admix_K{K}_seed{seed}.log",
            K=ADMIX_K,
            seed=ADMIX_SEEDS,
        ),
    output:
        best_k="results/{refGenome}/structure/admixture/logs/{prefix}_best_K.txt",
        cv_plot="results/{refGenome}/structure/admixture/logs/{prefix}_admixture_cv.png",
        deltak_plot="results/{refGenome}/structure/admixture/logs/{prefix}_admixture_deltak.png",
    params:
        species=lambda wc: wc.prefix,
        # admixture_cloud.R globs *.log in this directory; it must end in "/".
        output_path=lambda wc: f"results/{wc.refGenome}/structure/admixture/logs/",
    conda:
        "../envs/pie.yml"
    log:
        "logs/{refGenome}/structure/{prefix}_best_k.log",
    script:
        "../scripts/admixture_cloud.R"


if RUN_TESS:

    rule algatr_env:
        """Install the R packages that are not on conda (algatr, wingen, tess3r).

        Kept as its own rule because the install is slow and easy to get wrong;
        doing it once up front makes failures obvious.
        """
        input:
            vcf="results/{refGenome}/prepare/{prefix}_subsampled.vcf.gz",
        output:
            done="results/{refGenome}/structure/{prefix}_algatr_env.txt",
        conda:
            "../envs/algatr.yml"
        log:
            "logs/{refGenome}/structure/{prefix}_algatr_env.log",
        script:
            "../scripts/setup_algatr_env.R"

    rule tess:
        """Spatially-explicit ancestry (TESS3), over tess_kvals."""
        input:
            vcf="results/{refGenome}/prepare/{prefix}_subsampled.vcf.gz",
            vcf_idx="results/{refGenome}/prepare/{prefix}_subsampled.vcf.gz.csi",
            eigenvec="results/{refGenome}/structure/{prefix}.eigenvec",
            env="results/{refGenome}/structure/{prefix}_algatr_env.txt",
            coords="results/{refGenome}/prepare/{prefix}.coords.txt",
        output:
            qmatrix="results/{refGenome}/structure/{prefix}_TESS_qmatrix.csv",
            xval="results/{refGenome}/structure/{prefix}_TESS_xval.csv",
            peakram="results/{refGenome}/structure/{prefix}_TESS_peakRAM.csv",
            xval_plot="results/{refGenome}/structure/{prefix}_TESS_xval.png",
        params:
            # positional, in this order: species, data_path, rmislands, kvals,
            # output_path, incl_env
            species=lambda wc: wc.prefix,
            data_path=lambda wc: f"results/{wc.refGenome}/",
            rmislands="false",
            kvals=config["tess_kvals"],
            output_path=lambda wc: f"results/{wc.refGenome}/structure/",
            incl_env="false",
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000,
        conda:
            "../envs/algatr.yml"
        log:
            "logs/{refGenome}/structure/{prefix}_tess.log",
        script:
            "../scripts/TESS_cloud.R"

    rule tess_pieplot:
        input:
            tess="results/{refGenome}/structure/{prefix}_TESS_qmatrix.csv",
            coords="results/{refGenome}/prepare/{prefix}.coords.txt",
        output:
            plot="results/{refGenome}/structure/{prefix}_TESS_pieplot.pdf",
        params:
            species=lambda wc: wc.prefix,
        conda:
            "../envs/pie.yml"
        log:
            "logs/{refGenome}/structure/{prefix}_tess_pieplot.log",
        script:
            "../scripts/pie_plot.R"


rule structure_composite:
    """One figure summarising every ancestry result, for choosing K.

    Panels present depend on the mode: ADMIXTURE barplot, CV/deltaK and PCA are
    always drawn; maps need coordinates; the TESS panels and the
    ADMIXTURE-vs-TESS agreement panels need TESS.
    """
    input:
        unpack(
            lambda wc: {
                "cv_errors": f"results/{wc.refGenome}/structure/admixture/logs/{wc.prefix}_best_K.txt",
                "cv_plot": f"results/{wc.refGenome}/structure/admixture/logs/{wc.prefix}_admixture_cv.png",
                "deltak_plot": f"results/{wc.refGenome}/structure/admixture/logs/{wc.prefix}_admixture_deltak.png",
                "eigenvec": f"results/{wc.refGenome}/structure/{wc.prefix}.eigenvec",
                "eigenval": f"results/{wc.refGenome}/structure/{wc.prefix}.eigenval",
                "fam": f"results/{wc.refGenome}/structure/admixture/{wc.prefix}.fam",
                "coords": f"results/{wc.refGenome}/prepare/{wc.prefix}.coords.txt",
                "admixture_logs": expand(
                    "results/{refGenome}/structure/admixture/logs/{prefix}_admix_K{K}_seed{seed}.log",
                    refGenome=wc.refGenome, prefix=wc.prefix, K=ADMIX_K, seed=[PLOT_SEED],
                ),
                "admixture_q": expand(
                    "results/{refGenome}/structure/admixture/Q_files/{prefix}.{K}.{seed}.Q",
                    refGenome=wc.refGenome, prefix=wc.prefix, K=ADMIX_K, seed=[PLOT_SEED],
                ),
                **(
                    {
                        "tess": f"results/{wc.refGenome}/structure/{wc.prefix}_TESS_qmatrix.csv",
                        "xval": f"results/{wc.refGenome}/structure/{wc.prefix}_TESS_xval.csv",
                    }
                    if RUN_TESS
                    else {}
                ),
            }
        ),
    output:
        pdf="results/{refGenome}/structure/{prefix}_structure_composite.pdf",
        merged="results/{refGenome}/structure/{prefix}_admixture_tess_merged.csv",
    params:
        admixture_path=lambda wc: f"results/{wc.refGenome}/structure/admixture/Q_files/{wc.prefix}",
        plot_seed=PLOT_SEED,
        use_coords=USE_COORDS,
        run_tess=RUN_TESS,
        kmax=config["admixture_kmax"],
    conda:
        "../envs/pie.yml"
    log:
        "logs/{refGenome}/structure/{prefix}_structure_composite.log",
    script:
        "../scripts/structure_composite.R"
