"""Clean SNPs -> the analysis sample set and the SNP subsets used downstream.

This is where the workflow begins: results/{refGenome}/{prefix}_clean_snps.vcf.gz
is an input, not something produced here. Site filtering belongs to the variant
caller.

The one branch here is coordinates. With use_coords: true, samples without
coordinates are dropped (drop_samples_without_coords) so that every downstream
analysis — PCA, ADMIXTURE, TESS, ROH, pi, clam — sees exactly the same samples.
With use_coords: false the coordinate files are written empty and the analysis
set is the full clean SNP set.

Outputs:
  {prefix}_annotated.vcf.gz      clean SNPs with CHROM_POS_REF_ALT variant IDs
  {prefix}.coords.txt            sample <TAB> longitude <TAB> latitude (no header)
  {prefix}_analysis_set.vcf.gz   the sample set every analysis uses
  {prefix}_pruned_ld.vcf.gz      LD-pruned (ld_r2), for PCA/ADMIXTURE/TESS
  {prefix}_subsampled.vcf.gz     ~1M SNPs, stricter MAF, for ADMIXTURE/TESS
  {prefix}_dist.dist(.id)        plink identity-by-state distance matrices
"""


rule annotate_variant_ids:
    """Give every SNP a stable ID so LD pruning can be applied by ID later."""
    input:
        vcf="results/{refGenome}/{prefix}_clean_snps.vcf.gz",
    output:
        vcf="results/{refGenome}/prepare/{prefix}_annotated.vcf.gz",
        idx="results/{refGenome}/prepare/{prefix}_annotated.vcf.gz.csi",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/{refGenome}/prepare/{prefix}_annotate.log",
    shell:
        """
        exec 2> {log}
        bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \
          -O z -o {output.vcf} {input.vcf}
        bcftools index {output.vcf}
        """


rule vcf_sample_list:
    input:
        vcf="results/{refGenome}/prepare/{prefix}_annotated.vcf.gz",
    output:
        samples="results/{refGenome}/prepare/{prefix}_vcf_samples.txt",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/{refGenome}/prepare/{prefix}_sample_list.log",
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
        samples="results/{refGenome}/prepare/{prefix}_vcf_samples.txt",
        sheet=[config["coords_sheet"]] if USE_COORDS else [],
    output:
        coords="results/{refGenome}/prepare/{prefix}.coords.txt",
        no_coords="results/{refGenome}/prepare/{prefix}.no_coords.txt",
    params:
        coords_sheet=config["coords_sheet"],
        use_coords=USE_COORDS,
        drop_samples_without_coords=config["drop_samples_without_coords"],
    log:
        "logs/{refGenome}/prepare/{prefix}_prepare_coords.log",
    script:
        "../scripts/prepare_coords.py"


rule analysis_set:
    """The VCF every downstream analysis uses.

    The only thing dropped here is samples with no coordinates. Alleles are
    re-trimmed and AF recomputed so that sites which became invariant after
    dropping samples do not linger.

    Sample and site exclusion upstream of this is snpArcher's job. Sex-linked
    and anomalous scaffolds are detected and excluded later, in populations.

    When there is nothing to drop this is a straight copy, which keeps the rest
    of the workflow free of coordinate-mode branching.
    """
    input:
        vcf="results/{refGenome}/prepare/{prefix}_annotated.vcf.gz",
        idx="results/{refGenome}/prepare/{prefix}_annotated.vcf.gz.csi",
        no_coords="results/{refGenome}/prepare/{prefix}.no_coords.txt",
    output:
        vcf="results/{refGenome}/prepare/{prefix}_analysis_set.vcf.gz",
        idx="results/{refGenome}/prepare/{prefix}_analysis_set.vcf.gz.csi",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/{refGenome}/prepare/{prefix}_analysis_set.log",
    shell:
        """
        exec > {log} 2>&1
        if [ ! -s {input.no_coords} ]; then
            echo "no samples to drop; copying the annotated VCF through"
            cp {input.vcf} {output.vcf}
            cp {input.idx} {output.idx}
        else
            echo "dropping $(wc -l < {input.no_coords}) sample(s) without coordinates:"
            sed 's/^/  /' {input.no_coords}
            bcftools view -S ^{input.no_coords} -f .,PASS {input.vcf} -a -U -O u \
              | bcftools +fill-tags -Ou \
              | bcftools view -m2 -e 'AF==0 | ref="N" | ALT="."' -O z -o {output.vcf}
            bcftools index {output.vcf}
        fi
        """


rule ld_prune_list:
    input:
        vcf="results/{refGenome}/prepare/{prefix}_analysis_set.vcf.gz",
    output:
        prune_in="results/{refGenome}/prepare/{prefix}_ld.prune.in",
        prune_out="results/{refGenome}/prepare/{prefix}_ld.prune.out",
    params:
        out_prefix=lambda wc, output: output.prune_in[: -len(".prune.in")],
        r2=config["ld_r2"],
    resources:
        mem_mb=16000,
    conda:
        "../envs/plink.yml"
    log:
        "logs/{refGenome}/prepare/{prefix}_ld_prune_list.log",
    shell:
        """
        plink2 --vcf {input.vcf} --make-bed \
          --indep-pairwise 50kb {params.r2} \
          --out {params.out_prefix} \
          --allow-extra-chr --autosome-num 95 --const-fid --bad-freqs > {log} 2>&1
        """


rule ld_pruned_vcf:
    input:
        vcf="results/{refGenome}/prepare/{prefix}_analysis_set.vcf.gz",
        prune_in="results/{refGenome}/prepare/{prefix}_ld.prune.in",
    output:
        vcf="results/{refGenome}/prepare/{prefix}_pruned_ld.vcf.gz",
        idx="results/{refGenome}/prepare/{prefix}_pruned_ld.vcf.gz.csi",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/{refGenome}/prepare/{prefix}_ld_pruned_vcf.log",
    shell:
        """
        exec 2> {log}
        bcftools view -i 'ID=@{input.prune_in}' -O z -o {output.vcf} {input.vcf}
        bcftools index {output.vcf}
        """


rule distance_matrix:
    """plink1.9 identity-by-state distances, on the full and LD-pruned sets."""
    input:
        vcf="results/{refGenome}/prepare/{prefix}_analysis_set.vcf.gz",
        pruned="results/{refGenome}/prepare/{prefix}_pruned_ld.vcf.gz",
    output:
        dist="results/{refGenome}/prepare/{prefix}_dist.dist",
        dist_id="results/{refGenome}/prepare/{prefix}_dist.dist.id",
        dist_ld="results/{refGenome}/prepare/{prefix}_pruned_ld.dist",
        dist_ld_id="results/{refGenome}/prepare/{prefix}_pruned_ld.dist.id",
    params:
        out=lambda wc, output: output.dist[: -len(".dist")],
        out_ld=lambda wc, output: output.dist_ld[: -len(".dist")],
    resources:
        mem_mb=32000,
    conda:
        "../envs/plink.yml"
    log:
        "logs/{refGenome}/prepare/{prefix}_distance.log",
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
        vcf="results/{refGenome}/prepare/{prefix}_pruned_ld.vcf.gz",
    output:
        filtered=temp("results/{refGenome}/prepare/{prefix}_pruned_ld_maffilt.vcf.gz"),
        filtered_idx=temp("results/{refGenome}/prepare/{prefix}_pruned_ld_maffilt.vcf.gz.csi"),
        vcf="results/{refGenome}/prepare/{prefix}_subsampled.vcf.gz",
        idx="results/{refGenome}/prepare/{prefix}_subsampled.vcf.gz.csi",
    params:
        maf=config["subsample_maf"],
        miss=config["missingness"],
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/{refGenome}/prepare/{prefix}_subsample.log",
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
