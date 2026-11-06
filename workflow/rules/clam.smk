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
  {prefix}_summary_per_sample.csv                   genome-wide, per sample
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

    GVCF_SOURCE = (config["gvcf_source"] or "").replace("{prefix}", PREFIX)

    def _longest_contig(vcf):
        """Longest ##contig length in a bgzipped VCF's header (0 if none)."""
        import gzip
        import re
        longest = 0
        with gzip.open(vcf, "rt") as fh:
            for line in fh:
                if not line.startswith("##"):
                    break
                if line.startswith("##contig="):
                    m = re.search(r"length=(\d+)", line)
                    if m:
                        longest = max(longest, int(m.group(1)))
        return longest

    # A .tbi cannot address positions past 2^29 - 1, so references with a
    # longer contig (58-Anaxyrus, 58-Pseudacris) use .csi for the GVCFs, the
    # variants and the ROH BEDs: CLAM_IDX is the extension, TABIX_CSI the flag.
    # Reading .csi needs a clam with CSI support (cademirch/clam#64); point
    # clam_bin_dir at one until it is released.
    CLAM_IDX = "csi" if _longest_contig(INPUT_VCF) > 2**29 - 1 else "tbi"
    TABIX_CSI = "-C" if CLAM_IDX == "csi" else ""
    CLAM_PATH = f"export PATH={config['clam_bin_dir']}:$PATH;" if config.get("clam_bin_dir") else ""

    def clam_all_samples():
        """Every analysed sample, assigned to a population or not (NA).

        clam collect and loci cover all of them, so the per-sample callable
        masks do not depend on K: after a change of K only clam stat, which
        takes the populations from clam/samples.tsv, needs rerunning. The
        .fam is the analysis set's sample list and exists once the assignment
        checkpoint has run.
        """
        fam = Path(out("structure/admixture/{prefix}.fam"))
        # Gate on the assignment checkpoint only when the .fam does not exist
        # yet (a fresh run). Otherwise collect would depend on the checkpoint,
        # and rerunning the assignment (a change of K or exclude_samples) would
        # redo collect, i.e. download every GVCF again.
        if not fam.exists():
            active_k()
        return [line.split()[1] for line in fam.read_text().splitlines() if line.strip()]

    def clam_gvcf_inputs(wc):
        """Every GVCF (and index) clam will read, as declared inputs.

        Listing them individually rather than depending on the directory means a
        replaced GVCF triggers a rerun, and a missing one fails during DAG
        construction with a useful message rather than deep inside clam. With
        gvcf_source set they are produced by clam_fetch_gvcf instead, so they
        are not expected to exist yet.
        """
        samples = clam_all_samples()
        directory_ = gvcf_dir()
        if GVCF_SOURCE:
            return [
                os.path.join(directory_, f"{sample}.g.vcf.gz{ext}")
                for sample in samples
                for ext in ("", f".{CLAM_IDX}")
            ]
        files = []
        missing = []
        for sample in samples:
            gvcf = os.path.join(directory_, f"{sample}.g.vcf.gz")
            files.append(gvcf)
            if not os.path.exists(gvcf):
                missing.append(gvcf)
            tbi = f"{gvcf}.{CLAM_IDX}"
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

    if GVCF_SOURCE:

        rule clam_fetch_gvcf:
            """Download one sample's GVCF and index from gvcf_source.

            temp(): GVCFs are the largest inputs by far, and clam collect is
            their only reader, so each pair is deleted as soon as the depth
            store is built. Concurrency is capped by the gvcf_fetch resource;
            set it with --resources gvcf_fetch=N.
            """
            output:
                gvcf=temp(os.path.join(gvcf_dir(), "{sample}.g.vcf.gz")),
                tbi=temp(os.path.join(gvcf_dir(), "{sample}.g.vcf.gz." + CLAM_IDX)),
            wildcard_constraints:
                sample=r"[^/]+",
            params:
                source=GVCF_SOURCE.rstrip("/"),
                idx=CLAM_IDX,
            resources:
                gvcf_fetch=1,
            log:
                logp("clam/fetch/{sample}.log"),
            shell:
                """
                rsync -a --partial --timeout=600 \
                  {params.source}/{wildcards.sample}.g.vcf.gz \
                  {params.source}/{wildcards.sample}.g.vcf.gz.{params.idx} \
                  $(dirname {output.gvcf})/ > {log} 2>&1
                """

    rule clam_collect:
        """Per-sample per-base depth from the GVCFs into a zarr store."""
        input:
            gvcfs=clam_gvcf_inputs,
            # ancient(): a later change to the exclude list must not redo the
            # depth store (the GVCFs are gone by then). clam stat applies the
            # current list, and thresholds are per chromosome, so newly excluded
            # scaffolds drop out exactly; a scaffold no longer excluded is absent
            # until clam_collect is forced.
            exclude_file=ancient(out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt")),
        output:
            depths_zarr=_depths_zarr_output(DEPTHS_ZARR),
        params:
            min_gq=config["clam_min_gq"],
            # All analysed samples, NA included (see clam_all_samples). Written
            # to a table in the job rather than taken as an input file: a new
            # input file would make every finished project collect again.
            samples=lambda wc: " ".join(clam_all_samples()),
            gvcf_dir=gvcf_dir(),
            clam_path=CLAM_PATH,
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
            SAMPLES=$(mktemp)
            trap 'rm -f $SAMPLES' EXIT
            printf 'sample_name\tfile_path\n' > $SAMPLES
            for s in {params.samples}; do
              printf '%s\t%s/%s.g.vcf.gz\n' $s {params.gvcf_dir} $s >> $SAMPLES
            done
            {params.clam_path}
            clam collect \
              -o {output.depths_zarr} \
              --min-gq {params.min_gq} \
              --samples $SAMPLES \
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
            exclude_file=ancient(out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt")),  # see clam_collect
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
            exclude_file=ancient(out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt")),  # see clam_collect
            thresholds=out("clam/thresholds.tsv"),
        output:
            done=touch(out("clam/callable.zarr.done")),
        params:
            zarr=out("clam/callable.zarr"),
            clam_path=CLAM_PATH,
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
            {params.clam_path}
            clam loci \
              -o {params.zarr} \
              --per-sample \
              --thresholds-file {input.thresholds} \
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
            tbi=out("clam/variants.vcf.gz." + CLAM_IDX),
        params:
            csi=TABIX_CSI,
        conda:
            "../envs/clam.yml"
        log:
            logp("clam/{prefix}_prepare_vcf.log"),
        shell:
            """
            exec 2> {log}
            ln -sf $(readlink -f {input.vcf}) {output.vcf}
            tabix -f {params.csi} -p vcf {output.vcf}
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
            bed_tbi=out("clam/combined_roh.bed.gz." + CLAM_IDX),
        params:
            min_roh_size=config["roh_min_size"],
            csi=TABIX_CSI,
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
            tabix {params.csi} -p bed {output.bed_gz}
            """

    rule clam_stat:
        """Windowed pi and heterozygosity, plus dxy/Fst when there are >1 populations.

        clam writes fixed file names into -o, so they are renamed to carry the
        project prefix. dxy/fst only exist for K >= 2, hence the guarded moves.
        """
        input:
            zarr_done=out("clam/callable.zarr.done"),
            vcf=out("clam/variants.vcf.gz"),
            tbi=out("clam/variants.vcf.gz." + CLAM_IDX),
            roh_bed=out("clam/combined_roh.bed.gz"),
            roh_tbi=out("clam/combined_roh.bed.gz." + CLAM_IDX),
            exclude_file=out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt"),
            # Populations come from here, not from callable.zarr's metadata,
            # which keeps the ones it was built with: after a change of K only
            # clam stat then needs rerunning (the per-sample masks do not
            # depend on populations).
            samples_tsv=out("clam/samples.tsv"),
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
            clam_path=CLAM_PATH,
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
            {params.clam_path}
            clam stat \
              -o {params.outdir} \
              -w {params.window} \
              -c {params.zarr} \
              -r {input.roh_bed} \
              --samples {input.samples_tsv} \
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

        Per sample: heterozygosity and heterozygosity-outside-ROH. Per
        population: pi, pooled heterozygosity, and the mean of the samples'
        heterozygosity. Pairwise: dxy and Hudson's Fst as a ratio of averages.
        Everything is restricted to scaffolds >= min_scaffold_size, so pi, dxy
        and Fst share one set of sequence. Pairwise is header-only when K = 1.
        """
        input:
            pi=out("clam/{prefix}_pi.tsv"),
            het=out("clam/{prefix}_heterozygosity.tsv"),
            n_pops=out("clam/n_populations.txt"),
            samples_tsv=out("clam/samples.tsv"),
            stat_done=out("clam/{prefix}_clam_stat.done"),
            fai=CONTIGS,
        output:
            per_sample=out("clam/{prefix}_summary_per_sample.csv"),
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

    # --- heterozygosity outside ROH at other ROH cutoffs ----------------------
    # New rules rather than wildcards on clam_combine_roh/clam_stat, so adding
    # them changed nothing for runs that already have those outputs.

    rule clam_combine_roh_extra:
        """As clam_combine_roh, keeping runs of at least {size} bp."""
        input:
            pops_done=out("populations/{prefix}_populations-done.txt"),
            roh_files=clam_roh_inputs,
        output:
            bed_gz=out("clam/combined_roh_{size}.bed.gz"),
            bed_tbi=out("clam/combined_roh_{size}.bed.gz." + CLAM_IDX),
        wildcard_constraints:
            size=r"\d+",
        params:
            csi=TABIX_CSI,
        conda:
            "../envs/clam.yml"
        log:
            logp("clam/{prefix}_combine_roh_{size}.log"),
        shell:
            """
            exec 2> {log}
            cat {input.roh_files} \
              | awk -v m={wildcards.size} 'BEGIN{{OFS="\\t"}} NF>=4 && ($3-$2)>=m {{print $1,$2,$3,$4}}' \
              | sort -k1,1 -k2,2n \
              | bgzip -c > {output.bed_gz}
            tabix {params.csi} -p bed {output.bed_gz}
            """

    rule clam_het_roh_extra:
        """Per-window heterozygosity in and out of ROH >= {size} bp.

        `clam stat` again on the kept callable store and variants, with the
        other ROH BED. It also recomputes pi (and dxy/Fst), identical to
        clam_stat's; those are discarded.
        """
        input:
            zarr_done=out("clam/callable.zarr.done"),
            vcf=out("clam/variants.vcf.gz"),
            tbi=out("clam/variants.vcf.gz." + CLAM_IDX),
            roh_bed=out("clam/combined_roh_{size}.bed.gz"),
            roh_tbi=out("clam/combined_roh_{size}.bed.gz." + CLAM_IDX),
            exclude_file=out("populations/{prefix}_sex_id/{prefix}_exclude_scaffolds.txt"),
            samples_tsv=out("clam/samples.tsv"),  # populations; see clam_stat
        output:
            het=out("clam/{prefix}_heterozygosity_roh{size}.tsv"),
        wildcard_constraints:
            size=r"\d+",
        params:
            zarr=out("clam/callable.zarr"),
            outdir=out("clam/stat_roh{size}"),
            window=config["clam_window"],
            clam_path=CLAM_PATH,
        threads: config["clam_threads"]
        conda:
            "../envs/clam.yml"
        log:
            logp("clam/{prefix}_stat_roh{size}.log"),
        shell:
            """
            EXCLUDE_ARG=""
            if [ -s {input.exclude_file} ]; then
              EXCLUDE_ARG="--exclude-file {input.exclude_file}"
            fi
            rm -rf {params.outdir}; mkdir -p {params.outdir}
            {params.clam_path}
            clam stat \
              -o {params.outdir} \
              -w {params.window} \
              -c {params.zarr} \
              -r {input.roh_bed} \
              --samples {input.samples_tsv} \
              --force-samples \
              -t {threads} \
              $EXCLUDE_ARG \
              {input.vcf} > {log} 2>&1
            mv {params.outdir}/heterozygosity.tsv {output.het}
            rm -rf {params.outdir}
            """

    rule clam_summary_het_roh_extra:
        """Genome-wide heterozygosity outside ROH >= {size} bp, per sample and population."""
        input:
            het=out("clam/{prefix}_heterozygosity_roh{size}.tsv"),
            fai=CONTIGS,
        output:
            summary=out("clam/{prefix}_summary_het_roh{size}.csv"),
        wildcard_constraints:
            size=r"\d+",
        params:
            min_scaffold_size=config["min_scaffold_size"],
        conda:
            "../envs/clam.yml"
        log:
            logp("clam/{prefix}_summary_het_roh{size}.log"),
        script:
            "../scripts/summarize_clam_het.py"
