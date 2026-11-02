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
  {prefix}.coords.txt            sample <TAB> longitude <TAB> latitude (no header)
  {prefix}_analysis_set.vcf.gz   the sample set every analysis uses; a link to
                                 the input VCF when no samples are dropped
  {prefix}_pruned_ld.vcf.gz      LD-pruned (ld_r2), for PCA/ADMIXTURE/TESS, with
                                 CHROM_POS_REF_ALT variant IDs
  {prefix}_subsampled.vcf.gz     ~1M SNPs, stricter MAF, for ADMIXTURE/TESS (temp)
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


rule vcf_sample_list:
    input:
        vcf=INPUT_VCF,
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

    The only thing dropped is samples with no coordinates; when any are,
    alleles are re-trimmed and AF recomputed so that sites which became
    invariant do not linger. When none are, this is a link to the input VCF
    plus an index: rewriting it cost hours on large projects (58-Sceloporus:
    79 GB) and stored a second copy, only to add variant IDs that the LD-pruning
    steps now make themselves.

    Sample and site exclusion upstream of this is snpArcher's job. Sex-linked
    and anomalous scaffolds are detected and excluded later, in populations.
    """
    input:
        vcf=INPUT_VCF,
        no_coords=out("prepare/{prefix}.no_coords.txt"),
    output:
        vcf=out("prepare/{prefix}_analysis_set.vcf.gz"),
        idx=out("prepare/{prefix}_analysis_set.vcf.gz.csi"),
    # Compression is ~80% of a rewrite; 4 threads made it ~3x faster on
    # 58-Sceloporus, and more barely helped.
    threads: 4
    conda:
        "../envs/bcftools.yml"
    log:
        logp("prepare/{prefix}_analysis_set.log"),
    shell:
        """
        exec > {log} 2>&1
        if [ ! -s {input.no_coords} ]; then
            echo "no samples to drop; linking the input VCF"
            ln -sf $(readlink -f {input.vcf}) {output.vcf}
        else
            echo "dropping $(wc -l < {input.no_coords}) sample(s) without coordinates:"
            sed 's/^/  /' {input.no_coords}
            bcftools view --threads {threads} -S ^{input.no_coords} -f .,PASS -a -U -O u {input.vcf} \
              | bcftools +fill-tags -Ou \
              | bcftools view --threads {threads} -m2 -e 'AF==0 | ref="N" | ALT="."' -O z -o {output.vcf}
        fi
        bcftools index --threads {threads} -f {output.vcf}
        """


rule ld_prune_list:
    input:
        vcf=out("prepare/{prefix}_analysis_set.vcf.gz"),
        contigs=CONTIGS,
    output:
        prune_in=temp(out("prepare/{prefix}_ld.prune.in")),
        prune_out=temp(out("prepare/{prefix}_ld.prune.out")),
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
        exec > {log} 2>&1
        # IDs as CHROM_POS_REF_ALT, matching what ld_pruned_vcf sets. plink2
        # rewrites contig names it reads as chromosome codes (1, chr1, X, MT...);
        # ld_pruned_vcf checks that every listed ID is found.
        # --output-chr keeps the VCF's style for the names plink2 recodes: chr
        # prefix or not, M or MT. Other contig names pass through unchanged.
        OUTPUT_CHR=$(awk '$1 ~ /^chr([0-9]+|X|Y|XY|M|MT)$/ {{p = "chr"}} $1 ~ /^(chr)?M$/ {{m = "M"}}
                          END {{print p (m ? m : "MT")}}' {input.contigs})
        plink2 --vcf {input.vcf} --threads {threads} \
          --set-all-var-ids '@_#_$r_$a' --new-id-max-allele-len 1000 \
          --output-chr $OUTPUT_CHR \
          --indep-pairwise 50kb {params.r2} \
          --out {params.out_prefix} \
          --allow-extra-chr --autosome-num 95 --const-fid --bad-freqs
        """


rule ld_pruned_vcf:
    input:
        vcf=out("prepare/{prefix}_analysis_set.vcf.gz"),
        prune_in=out("prepare/{prefix}_ld.prune.in"),
    output:
        vcf=out("prepare/{prefix}_pruned_ld.vcf.gz"),
        idx=out("prepare/{prefix}_pruned_ld.vcf.gz.csi"),
    threads: 4
    conda:
        "../envs/bcftools.yml"
    log:
        logp("prepare/{prefix}_ld_pruned_vcf.log"),
    shell:
        """
        exec 2> {log}
        bcftools annotate --threads {threads} --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -O u {input.vcf} \
          | bcftools view --threads {threads} -i 'ID=@{input.prune_in}' -O z -o {output.vcf}
        bcftools index --threads {threads} {output.vcf}
        # Every SNP plink2 kept must be found. Fewer means the IDs did not match
        # (plink2 renamed a contig), which would silently thin the data.
        KEPT=$(bcftools index -n {output.vcf})
        LISTED=$(wc -l < {input.prune_in})
        echo "plink2 kept $LISTED SNPs; $KEPT found in the VCF" >&2
        if [ "$KEPT" -lt "$LISTED" ]; then
            echo "ERROR: $(( LISTED - KEPT )) pruned SNP IDs not found; see ld_prune_list's contig names" >&2
            exit 1
        fi
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
        # Only PCA, the ADMIXTURE plink set and TESS read this.
        vcf=temp(out("prepare/{prefix}_subsampled.vcf.gz")),
        idx=temp(out("prepare/{prefix}_subsampled.vcf.gz.csi")),
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
