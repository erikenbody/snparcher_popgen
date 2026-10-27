"""Clean SNPs -> the analysis sample set and the SNP subsets used downstream.

This is where the workflow begins: the `vcf` named in the config is an input,
not something produced here. Site filtering belongs to the variant caller.

The one branch here is coordinates. With use_coords: true, samples without
coordinates are dropped (drop_samples_without_coords) so that every downstream
analysis — PCA, ADMIXTURE, TESS, ROH, pi, clam — sees exactly the same samples.
With use_coords: false the coordinate files are written empty and the analysis
set is the full clean SNP set.

Outputs, under {outdir}/prepare/:
  {prefix}.contigs.tsv           contig <TAB> length, in reference order
  {prefix}_annotated.vcf.gz      clean SNPs with CHROM_POS_REF_ALT variant IDs
  {prefix}.coords.txt            sample <TAB> longitude <TAB> latitude (no header)
  {prefix}_analysis_set.vcf.gz   the sample set every analysis uses
  {prefix}_pruned_ld.vcf.gz      LD-pruned (ld_r2), for PCA/ADMIXTURE/TESS
  {prefix}_subsampled.vcf.gz     ~1M SNPs, stricter MAF, for ADMIXTURE/TESS
  {prefix}_dist.dist(.id)        plink identity-by-state distance matrices
"""


rule contig_lengths:
    """Contig names and lengths, in reference order.

    Taken from the .fai when `fai` is set, otherwise from the VCF's ##contig
    header lines — which GATK, Sentieon and bcftools write for every reference
    contig, including ones left with no variants after filtering. Either way
    the result is the first two .fai columns, which is all anything downstream
    reads. Fails if the header has no contig lengths, or if a contig carrying
    variants is missing from it, rather than silently shrinking a denominator.
    """
    input:
        vcf=INPUT_VCF,
        fai=[INPUT_FAI] if INPUT_FAI else [],
    output:
        contigs=CONTIGS,
    params:
        from_fai=bool(INPUT_FAI),
    conda:
        "../envs/bcftools.yml"
    log:
        logp("prepare/{prefix}_contig_lengths.log"),
    shell:
        """
        exec 2> {log}
        if [ "{params.from_fai}" = "True" ]; then
            cut -f1,2 {input.fai} > {output.contigs}
        else
            bcftools view -h {input.vcf} \
              | awk '/^##contig=/ {{
                    id = ""; len = ""
                    if (match($0, /ID=[^,>]+/))      id  = substr($0, RSTART + 3, RLENGTH - 3)
                    if (match($0, /length=[0-9]+/))  len = substr($0, RSTART + 7, RLENGTH - 7)
                    if (id == "" || len == "") {{ print "contig header line without ID or length: " $0 > "/dev/stderr"; bad = 1; exit 1 }}
                    print id "\\t" len
                }}
                END {{ exit bad }}' > {output.contigs}
            if [ ! -s {output.contigs} ]; then
                echo "{input.vcf} has no ##contig lines; set fai in the config" >&2
                exit 1
            fi
        fi

        # Every contig that carries a variant must have a length.
        MISSING=$(bcftools query -f '%CHROM\\n' {input.vcf} | uniq | sort -u \
                    | comm -23 - <(cut -f1 {output.contigs} | sort -u) | head -5)
        if [ -n "$MISSING" ]; then
            echo "contigs with variants but no length: $MISSING" >&2
            exit 1
        fi
        echo "$(wc -l < {output.contigs}) contigs"
        """


rule annotate_variant_ids:
    """Give every SNP a stable ID so LD pruning can be applied by ID later."""
    input:
        vcf=INPUT_VCF,
    output:
        vcf=out("prepare/{prefix}_annotated.vcf.gz"),
        idx=out("prepare/{prefix}_annotated.vcf.gz.csi"),
    conda:
        "../envs/bcftools.yml"
    log:
        logp("prepare/{prefix}_annotate.log"),
    shell:
        """
        exec 2> {log}
        bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \
          -O z -o {output.vcf} {input.vcf}
        bcftools index {output.vcf}
        """


rule vcf_sample_list:
    input:
        vcf=out("prepare/{prefix}_annotated.vcf.gz"),
    output:
        samples=out("prepare/{prefix}_vcf_samples.txt"),
    conda:
        "../envs/bcftools.yml"
    log:
        logp("prepare/{prefix}_sample_list.log"),
    shell:
        "bcftools query -l {input.vcf} > {output.samples} 2> {log}"


rule prepare_coords:
    """Validate the user coordinate sheet against the VCF sample list.

    Writes the internal coordinate file consumed by TESS and every map panel,
    plus the list of samples that have no coordinates. Both are empty when
    use_coords is false. See scripts/prepare_coords.py for the accepted header
    spellings and the validation it applies.
    """
    input:
        samples=out("prepare/{prefix}_vcf_samples.txt"),
        sheet=[config["coords_sheet"]] if USE_COORDS else [],
    output:
        coords=out("prepare/{prefix}.coords.txt"),
        no_coords=out("prepare/{prefix}.no_coords.txt"),
    params:
        coords_sheet=config["coords_sheet"],
        use_coords=USE_COORDS,
        drop_samples_without_coords=config["drop_samples_without_coords"],
    log:
        logp("prepare/{prefix}_prepare_coords.log"),
    script:
        "../scripts/prepare_coords.py"


rule analysis_set:
    """The VCF every downstream analysis uses.

    The only thing dropped here is samples with no coordinates. Alleles are
    re-trimmed and AF recomputed so that sites which became invariant after
    dropping samples do not linger.

    Sample and site exclusion upstream of this is snpArcher's job. Sex-linked
    and anomalous scaffolds are detected and excluded later, in populations.

    When there is nothing to drop the annotated VCF is symlinked through rather
    than copied, which keeps the rest of the workflow free of coordinate-mode
    branching without storing the largest file twice.
    """
    input:
        vcf=out("prepare/{prefix}_annotated.vcf.gz"),
        idx=out("prepare/{prefix}_annotated.vcf.gz.csi"),
        no_coords=out("prepare/{prefix}.no_coords.txt"),
    output:
        vcf=out("prepare/{prefix}_analysis_set.vcf.gz"),
        idx=out("prepare/{prefix}_analysis_set.vcf.gz.csi"),
    conda:
        "../envs/bcftools.yml"
    log:
        logp("prepare/{prefix}_analysis_set.log"),
    shell:
        """
        exec > {log} 2>&1
        if [ ! -s {input.no_coords} ]; then
            echo "no samples to drop; linking the annotated VCF through"
            ln -sf $(basename {input.vcf}) {output.vcf}
            ln -sf $(basename {input.idx}) {output.idx}
        else
            echo "dropping $(wc -l < {input.no_coords}) sample(s) without coordinates:"
            sed 's/^/  /' {input.no_coords}
            bcftools view -S ^{input.no_coords} -f .,PASS {input.vcf} -a -U -O u \
              | bcftools +fill-tags -Ou \
              | bcftools view -m2 -e 'AF==0 | ref="N" | ALT="."' -O z -o {output.vcf}
            bcftools index -f {output.vcf}
        fi
        """


rule ld_prune_list:
    input:
        vcf=out("prepare/{prefix}_analysis_set.vcf.gz"),
    output:
        prune_in=out("prepare/{prefix}_ld.prune.in"),
        prune_out=out("prepare/{prefix}_ld.prune.out"),
    params:
        out_prefix=lambda wc, output: output.prune_in[: -len(".prune.in")],
        r2=config["ld_r2"],
    threads: 4
    resources:
        mem_mb=16000,
    conda:
        "../envs/plink.yml"
    log:
        logp("prepare/{prefix}_ld_prune_list.log"),
    shell:
        """
        plink2 --vcf {input.vcf} --threads {threads} \
          --indep-pairwise 50kb {params.r2} \
          --out {params.out_prefix} \
          --allow-extra-chr --autosome-num 95 --const-fid --bad-freqs > {log} 2>&1
        """


rule ld_pruned_vcf:
    input:
        vcf=out("prepare/{prefix}_analysis_set.vcf.gz"),
        prune_in=out("prepare/{prefix}_ld.prune.in"),
    output:
        vcf=out("prepare/{prefix}_pruned_ld.vcf.gz"),
        idx=out("prepare/{prefix}_pruned_ld.vcf.gz.csi"),
    conda:
        "../envs/bcftools.yml"
    log:
        logp("prepare/{prefix}_ld_pruned_vcf.log"),
    shell:
        """
        exec 2> {log}
        bcftools view -i 'ID=@{input.prune_in}' -O z -o {output.vcf} {input.vcf}
        bcftools index {output.vcf}
        """


rule distance_matrix:
    """plink1.9 identity-by-state distances, on the full and LD-pruned sets."""
    input:
        vcf=out("prepare/{prefix}_analysis_set.vcf.gz"),
        pruned=out("prepare/{prefix}_pruned_ld.vcf.gz"),
    output:
        dist=out("prepare/{prefix}_dist.dist"),
        dist_id=out("prepare/{prefix}_dist.dist.id"),
        dist_ld=out("prepare/{prefix}_pruned_ld.dist"),
        dist_ld_id=out("prepare/{prefix}_pruned_ld.dist.id"),
    params:
        out=lambda wc, output: output.dist[: -len(".dist")],
        out_ld=lambda wc, output: output.dist_ld[: -len(".dist")],
    resources:
        mem_mb=32000,
    conda:
        "../envs/plink.yml"
    log:
        logp("prepare/{prefix}_distance.log"),
    shell:
        """
        exec > {log} 2>&1
        plink --vcf {input.vcf} --out {params.out} \
          --allow-extra-chr --autosome-num 95 --distance square --const-fid
        plink --vcf {input.pruned} --out {params.out_ld} \
          --allow-extra-chr --autosome-num 95 --distance square --const-fid
        """


rule subsample_snps:
    """Thin the LD-pruned set to roughly 1M SNPs for ADMIXTURE and TESS.

    A stricter MAF (subsample_maf) is applied first: rare variants carry little
    information for ancestry estimation and inflate runtime. If fewer than 1M
    SNPs survive, all of them are kept.
    """
    input:
        vcf=out("prepare/{prefix}_pruned_ld.vcf.gz"),
    output:
        filtered=temp(out("prepare/{prefix}_pruned_ld_maffilt.vcf.gz")),
        filtered_idx=temp(out("prepare/{prefix}_pruned_ld_maffilt.vcf.gz.csi")),
        vcf=out("prepare/{prefix}_subsampled.vcf.gz"),
        idx=out("prepare/{prefix}_subsampled.vcf.gz.csi"),
    params:
        maf=config["subsample_maf"],
        miss=config["missingness"],
    conda:
        "../envs/bcftools.yml"
    log:
        logp("prepare/{prefix}_subsample.log"),
    shell:
        """
        exec > {log} 2>&1
        bcftools view \
          -e 'ALT="*" | F_MISSING > {params.miss} | ref="N" | AF < {params.maf} | AF > 1-{params.maf}' \
          {input.vcf} -O z -o {output.filtered}
        bcftools index {output.filtered}

        ALLSITES=$(bcftools query -f '%CHROM\\t%POS\\n' {output.filtered} | wc -l)
        WINDOW=$(( ALLSITES / 1000000 ))
        echo "sites after MAF filter: $ALLSITES; thinning window: $WINDOW"

        if [ "$WINDOW" -gt 1 ]; then
            bcftools +prune -w $WINDOW -n 1 -N rand -O z -o {output.vcf} {output.filtered}
        else
            bcftools view -O z -o {output.vcf} {output.filtered}
        fi
        bcftools index {output.vcf}
        """
