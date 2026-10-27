"""clam: callable-site-aware pi, dxy, Fst and ROH-corrected heterozygosity.

Everything here reads per-sample GVCFs rather than the joint VCF, which is what
makes the statistics comparable across projects: a site that was not callable in
a sample is excluded from that sample's denominator instead of being silently
treated as homozygous-reference.

  clam collect   per-sample per-base depth from the GVCFs -> depths zarr
  clam loci      depth -> per-sample callable masks, using per-chromosome
                 Poisson-SD depth bounds
  clam stat      windowed pi / dxy / Fst, and heterozygosity in and out of ROH

Population membership flows from the assignment checkpoint into samples.tsv, and
from there through collect -> loci -> stat via zarr metadata.

The depths zarr is big: roughly n_samples * genome_length * 4 bytes before
compression. It is deleted once callable loci are computed unless
clam_keep_depths_zarr is set; point clam_depths_zarr_dir at bulk storage for
large projects.

Outputs, under {outdir}/clam/:
  {prefix}_pi.tsv, {prefix}_heterozygosity.tsv      per window
  {prefix}_dxy.tsv, {prefix}_fst.tsv                per window, K >= 2 only
  {prefix}_summary_per_pop.csv                      genome-wide, per population
  {prefix}_summary_pairwise.csv                     genome-wide, per population pair
  {prefix}_fst_dxy_manhattan.png                    K >= 2 only
"""


if RUN_CLAM:

    DEPTHS_ZARR = depths_zarr_dir()

    def _depths_zarr_output(path):
        """temp() unless the user asked to keep the depth store."""
        d = directory(path)
        return d if as_bool(config["clam_keep_depths_zarr"]) else temp(d)

    def clam_gvcf_inputs(wc):
        """Every GVCF (and index) clam will read, as declared inputs.

        Listing them individually rather than depending on the directory means a
        replaced GVCF triggers a rerun, and a missing one fails during DAG
        construction with a useful message rather than deep inside clam.
        """
        active_k()  # gate on the assignment checkpoint
        pops = population_files()
        samples = [s for ss in pops.values() for s in ss]
        directory_ = gvcf_dir()
        files = []
        missing = []
        for sample in samples:
            gvcf = os.path.join(directory_, f"{sample}.g.vcf.gz")
            files.append(gvcf)
            if not os.path.exists(gvcf):
                missing.append(gvcf)
            tbi = f"{gvcf}.tbi"
            if os.path.exists(tbi):
                files.append(tbi)
            else:
                missing.append(tbi)
        if missing:
            shown = "\n  ".join(missing[:10])
            raise WorkflowError(
                f"{len(missing)} GVCF file(s) or index(es) missing from "
                f"gvcf_dir ({directory_}):\n  {shown}"
                + (f"\n  ... and {len(missing) - 10} more" if len(missing) > 10 else "")
                + "\nclam expects one bgzipped, tabix-indexed GVCF per sample named "
                "{sample}.g.vcf.gz, with sample names matching the VCF."
            )
        return files

    rule clam_samples_tsv:
        """clam's --samples table: sample_name, file_path, population."""
        input:
            pops_done=out("populations/{prefix}_populations-done.txt"),
        output:
            samples_tsv=out("clam/samples.tsv"),
            n_pops=out("clam/n_populations.txt"),
        run:
            pops = population_files()
            directory_ = gvcf_dir()
            rows = ["sample_name\tfile_path\tpopulation"]
            for pop, samples in sorted(pops.items(), key=lambda kv: int(kv[0])):
                for sample in samples:
                    path = os.path.join(directory_, f"{sample}.g.vcf.gz")
                    rows.append(f"{sample}\t{path}\tpop{pop}")
            Path(output.samples_tsv).write_text("\n".join(rows) + "\n")
            Path(output.n_pops).write_text(f"{len(pops)}\n")

    rule clam_collect:
        """Per-sample per-base depth from the GVCFs into a zarr store."""
        input:
            gvcfs=clam_gvcf_inputs,
            samples_tsv=out("clam/samples.tsv"),
            exclude_file=out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt"),
        output:
            depths_zarr=_depths_zarr_output(DEPTHS_ZARR),
        params:
            min_gq=config["clam_min_gq"],
        threads: config["clam_threads"]
        conda:
            "../envs/clam.yml"
        log:
            logp("clam/{prefix}_collect.log"),
        shell:
            """
            rm -rf {output.depths_zarr}
            EXCLUDE_ARG=""
            if [ -s {input.exclude_file} ]; then
              EXCLUDE_ARG="--exclude-file {input.exclude_file}"
            fi
            clam collect \
              -o {output.depths_zarr} \
              --min-gq {params.min_gq} \
              --samples {input.samples_tsv} \
              -t {threads} \
              $EXCLUDE_ARG > {log} 2>&1
            """

    rule clam_thresholds:
        """Per-chromosome callable-depth bounds, mean +/- N*sqrt(mean).

        Poisson SD, measured across all samples and positions in the depth
        store. This mirrors the per-chromosome coverage filter snpArcher
        applies, except that clam applies the bounds per individual.
        """
        input:
            depths_zarr=DEPTHS_ZARR,
            exclude_file=out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt"),
        output:
            thresholds=out("clam/thresholds.tsv"),
        params:
            sd_scale=config["clam_mean_depth_sd_scale"],
        threads: 4
        conda:
            "../envs/clam.yml"
        log:
            logp("clam/{prefix}_thresholds.log"),
        script:
            "../scripts/compute_thresholds.py"

    rule clam_sample_depth:
        """Per-sample mean coverage, genome-wide and per contig.

        Deliberately not wired into any target: it consumes the depth store,
        which is temporary, so requesting it after a completed run would force
        `clam collect` to run again. Ask for it explicitly on a fresh run:

          snakemake <outdir>/clam/<prefix>_sample_depth.csv
        """
        input:
            depths_zarr=DEPTHS_ZARR,
        output:
            summary=out("clam/{prefix}_sample_depth.csv"),
            per_contig=out("clam/{prefix}_sample_depth_per_contig.csv"),
        threads: 4
        conda:
            "../envs/clam.yml"
        log:
            logp("clam/{prefix}_sample_depth.log"),
        script:
            "../scripts/compute_sample_depth.py"

    rule clam_loci:
        """Per-sample callable masks from the depth store.

        --per-sample is required for accurate per-sample
        heterozygosity-outside-ROH; without it clam has to approximate the
        callable-not-in-ROH denominator. It costs roughly n_samples times more
        disk in callable.zarr.
        """
        input:
            depths_zarr=DEPTHS_ZARR,
            samples_tsv=out("clam/samples.tsv"),
            exclude_file=out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt"),
            thresholds=out("clam/thresholds.tsv"),
        output:
            done=touch(out("clam/callable.zarr.done")),
        params:
            zarr=out("clam/callable.zarr"),
        threads: config["clam_threads"]
        conda:
            "../envs/clam.yml"
        log:
            logp("clam/{prefix}_loci.log"),
        shell:
            """
            rm -rf {params.zarr}
            EXCLUDE_ARG=""
            if [ -s {input.exclude_file} ]; then
              EXCLUDE_ARG="--exclude-file {input.exclude_file}"
            fi
            clam loci \
              -o {params.zarr} \
              --per-sample \
              --thresholds-file {input.thresholds} \
              --samples {input.samples_tsv} \
              -t {threads} \
              $EXCLUDE_ARG \
              {input.depths_zarr} > {log} 2>&1
            """

    rule clam_prepare_vcf:
        """clam wants a tabix index; the analysis VCF carries a .csi."""
        input:
            vcf=out("prepare/{prefix}_analysis_set.vcf.gz"),
        output:
            vcf=out("clam/variants.vcf.gz"),
            tbi=out("clam/variants.vcf.gz.tbi"),
        conda:
            "../envs/clam.yml"
        log:
            logp("clam/{prefix}_prepare_vcf.log"),
        shell:
            """
            exec 2> {log}
            ln -sf $(readlink -f {input.vcf}) {output.vcf}
            tabix -f -p vcf {output.vcf}
            """

    def clam_roh_inputs(wc):
        return expand(out("populations/{prefix}_merged_roh/population_{k}.rg.roh"), k=active_k())

    rule clam_combine_roh:
        """All populations' merged ROH as one sorted, indexed BED.

        Runs shorter than roh_min_size are dropped so that clam's
        heterozygosity-outside-ROH uses the same ROH definition as F_ROH.
        """
        input:
            pops_done=out("populations/{prefix}_populations-done.txt"),
            roh_files=clam_roh_inputs,
        output:
            bed_gz=out("clam/combined_roh.bed.gz"),
            bed_tbi=out("clam/combined_roh.bed.gz.tbi"),
        params:
            min_roh_size=config["roh_min_size"],
        conda:
            "../envs/clam.yml"
        log:
            logp("clam/{prefix}_combine_roh.log"),
        shell:
            """
            exec 2> {log}
            cat {input.roh_files} \
              | awk -v m={params.min_roh_size} 'BEGIN{{OFS="\\t"}} NF>=4 && ($3-$2)>=m {{print $1,$2,$3,$4}}' \
              | sort -k1,1 -k2,2n \
              | bgzip -c > {output.bed_gz}
            tabix -p bed {output.bed_gz}
            """

    rule clam_stat:
        """Windowed pi and heterozygosity, plus dxy/Fst when there are >1 populations.

        clam writes fixed file names into -o, so they are renamed to carry the
        project prefix. dxy/fst only exist for K >= 2, hence the guarded moves.
        """
        input:
            zarr_done=out("clam/callable.zarr.done"),
            vcf=out("clam/variants.vcf.gz"),
            tbi=out("clam/variants.vcf.gz.tbi"),
            roh_bed=out("clam/combined_roh.bed.gz"),
            roh_tbi=out("clam/combined_roh.bed.gz.tbi"),
            exclude_file=out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt"),
        output:
            pi=out("clam/{prefix}_pi.tsv"),
            het=out("clam/{prefix}_heterozygosity.tsv"),
            done=touch(out("clam/{prefix}_clam_stat.done")),
        params:
            zarr=out("clam/callable.zarr"),
            outdir=out("clam"),
            dxy=out("clam/{prefix}_dxy.tsv"),
            fst=out("clam/{prefix}_fst.tsv"),
            window=config["clam_window"],
        threads: config["clam_threads"]
        conda:
            "../envs/clam.yml"
        log:
            logp("clam/{prefix}_stat.log"),
        shell:
            """
            EXCLUDE_ARG=""
            if [ -s {input.exclude_file} ]; then
              EXCLUDE_ARG="--exclude-file {input.exclude_file}"
            fi
            clam stat \
              -o {params.outdir} \
              -w {params.window} \
              -c {params.zarr} \
              -r {input.roh_bed} \
              --force-samples \
              -t {threads} \
              $EXCLUDE_ARG \
              {input.vcf} > {log} 2>&1

            mv {params.outdir}/pi.tsv             {output.pi}
            mv {params.outdir}/heterozygosity.tsv {output.het}
            [ -f {params.outdir}/dxy.tsv ] && mv {params.outdir}/dxy.tsv {params.dxy} || true
            [ -f {params.outdir}/fst.tsv ] && mv {params.outdir}/fst.tsv {params.fst} || true
            """

    rule clam_summary:
        """Genome-wide summaries from the per-window tables.

        Per population: pi, heterozygosity and heterozygosity-outside-ROH,
        restricted to scaffolds >= min_scaffold_size so that short scaffolds do
        not dominate. Pairwise: weighted dxy and Hudson's Fst from genome-wide
        pi and dxy, over all scaffolds. Pairwise is header-only when K = 1.
        """
        input:
            pi=out("clam/{prefix}_pi.tsv"),
            het=out("clam/{prefix}_heterozygosity.tsv"),
            n_pops=out("clam/n_populations.txt"),
            samples_tsv=out("clam/samples.tsv"),
            stat_done=out("clam/{prefix}_clam_stat.done"),
            fai=CONTIGS,
        output:
            per_pop=out("clam/{prefix}_summary_per_pop.csv"),
            pairwise=out("clam/{prefix}_summary_pairwise.csv"),
        params:
            dxy=out("clam/{prefix}_dxy.tsv"),
            fst=out("clam/{prefix}_fst.tsv"),
            min_scaffold_size=config["min_scaffold_size"],
        conda:
            "../envs/clam.yml"
        log:
            logp("clam/{prefix}_summary.log"),
        script:
            "../scripts/summarize_clam.py"

    rule clam_manhattan:
        """Genome-wide Fst and dxy per population pair.

        Only meaningful for K >= 2; single-population runs just touch the
        sentinel. fst/dxy are params, not inputs, because they do not exist
        when K = 1.
        """
        input:
            stat_done=out("clam/{prefix}_clam_stat.done"),
            per_pop=out("clam/{prefix}_summary_per_pop.csv"),
            fai=CONTIGS,
        output:
            done=touch(out("clam/{prefix}_manhattan.done")),
        params:
            prefix=PREFIX,
            outdir=out("clam"),
            dxy=out("clam/{prefix}_dxy.tsv"),
            fst=out("clam/{prefix}_fst.tsv"),
        conda:
            "../envs/plot.yml"
        log:
            logp("clam/{prefix}_manhattan.log"),
        script:
            "../scripts/plot_clam_manhattan.py"
