"""UCSC Genome Browser track hub: windowed diversity, depth, allele frequency,
non-callable regions, and the SNPs themselves. Included only when run_trackhub
is true.

Ported from snpArcher's v1 trackhub module, with its later fixes (visibility
full/hide, the hub.txt email line, the published citation, no scratch files in
the working directory). The hub uses `useOneFile on` and points at UCSC's
GenArk assembly for `refgenome`, so refgenome must be the assembly accession
(e.g. GCA_024610735.1). Everything the hub references is written into one
directory, {outdir}/trackhub/, which can be copied to a web server as is.
"""

TRACKHUB_WINDOWS = [1000, 10000, 100000]
# Don't change this unless you add rules for more statistics.
TRACKHUB_FILE_TYPES = ["Tajima", "SNP-Density", "Pi"]


wildcard_constraints:
    window=r"\d+",
    file_type="|".join(TRACKHUB_FILE_TYPES),


rule trackhub_write_hub:
    input:
        window_stat_files=expand(
            out("trackhub/{file_type}_{window}.bw"),
            file_type=TRACKHUB_FILE_TYPES, window=TRACKHUB_WINDOWS,
        ),
        vcf=out("trackhub/{prefix}_snps.vcf.gz"),
        tbi=out("trackhub/{prefix}_snps.vcf.gz.tbi"),
        callable_sites=out("trackhub/non_callable_sites.bb"),
        allele_freq=out("trackhub/allele_freq.bw"),
        depth=out("trackhub/depth.bw"),
    output:
        trackhub_file=out("trackhub/hub.txt"),
        html=out("trackhub/index.html"),
    params:
        refGenome=config["refgenome"],
        file_types=TRACKHUB_FILE_TYPES,
        windows=TRACKHUB_WINDOWS,
        email=config["trackhub_email"],
    log:
        logp("trackhub/{prefix}_write_hub.log"),
    script:
        "../scripts/write_hub_files.py"


rule trackhub_strip_vcf:
    """The VCF track: SNPs with INFO and all FORMAT fields but GT removed, to keep
    the file the browser fetches small."""
    input:
        vcf=INPUT_VCF,
    output:
        vcf=out("trackhub/{prefix}_snps.vcf.gz"),
        tbi=out("trackhub/{prefix}_snps.vcf.gz.tbi"),
    conda:
        "../envs/trackhub.yml"
    log:
        logp("trackhub/{prefix}_strip_vcf.log"),
    shell:
        """
        exec 2> {log}
        bcftools annotate -x INFO,^FORMAT/GT {input.vcf} -O z -o {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule trackhub_window_stat:
    """vcftools windowed Tajima's D, SNP density or pi."""
    input:
        vcf=INPUT_VCF,
    output:
        temp(out("trackhub/{window}.{file_type}")),
    params:
        flag=lambda wc: {
            "Tajima": "--TajimaD",
            "SNP-Density": "--SNPdensity",
            "Pi": "--window-pi",
        }[wc.file_type],
    conda:
        "../envs/trackhub.yml"
    log:
        logp("trackhub/{prefix}_{file_type}_{window}.log"),
    shell:
        "vcftools --gzvcf {input.vcf} {params.flag} {wildcards.window} --stdout > {output} 2> {log}"


rule trackhub_chrom_sizes:
    input:
        CONTIGS,
    output:
        out("trackhub/chrom.sizes"),
    shell:
        "cp {input} {output}"


rule trackhub_depth:
    input:
        vcf=INPUT_VCF,
        chrom_sizes=out("trackhub/chrom.sizes"),
    output:
        bg=temp(out("trackhub/depth.bg")),
        bw=out("trackhub/depth.bw"),
    conda:
        "../envs/trackhub.yml"
    log:
        logp("trackhub/{prefix}_depth.log"),
    shell:
        """
        exec 2> {log}
        bcftools query -f '%CHROM\\t%POS\\t%POS\\t%DP\\n' {input.vcf} \
          | awk -v OFS='\\t' '{{print $1,$2-1,$2,$4}}' > {output.bg}
        bedGraphToBigWig {output.bg} {input.chrom_sizes} {output.bw}
        """


rule trackhub_allele_freq:
    input:
        vcf=INPUT_VCF,
        chrom_sizes=out("trackhub/chrom.sizes"),
    output:
        bg=temp(out("trackhub/allele_freq.bg")),
        bw=out("trackhub/allele_freq.bw"),
    conda:
        "../envs/trackhub.yml"
    log:
        logp("trackhub/{prefix}_allele_freq.log"),
    shell:
        """
        exec 2> {log}
        vcftools --gzvcf {input.vcf} --freq2 --stdout | cut -f 1,2,6 | tail -n +2 \
          | awk -v OFS='\\t' '{{print $1,$2-1,$2,$3}}' > {output.bg}
        bedGraphToBigWig {output.bg} {input.chrom_sizes} {output.bw}
        """


rule trackhub_to_bedgraph:
    input:
        stat_file=out("trackhub/{window}.{file_type}"),
        chrom_sizes=out("trackhub/chrom.sizes"),
    output:
        temp(out("trackhub/{file_type}_{window}.bg")),
    log:
        logp("trackhub/{prefix}_{file_type}_{window}_bg.log"),
    script:
        "../scripts/vcftools_out_to_bg.py"


rule trackhub_to_bigwig:
    input:
        bg=out("trackhub/{file_type}_{window}.bg"),
        chrom_sizes=out("trackhub/chrom.sizes"),
    output:
        out("trackhub/{file_type}_{window}.bw"),
    conda:
        "../envs/trackhub.yml"
    log:
        logp("trackhub/{prefix}_{file_type}_{window}_bw.log"),
    shell:
        "bedGraphToBigWig {input.bg} {input.chrom_sizes} {output} 2> {log}"


rule trackhub_non_callable:
    """Complement of the callable sites, which is easier to read in a browser."""
    input:
        callable_sites=CALLABLE_BED,
        chrom_sizes=out("trackhub/chrom.sizes"),
    output:
        bb=out("trackhub/non_callable_sites.bb"),
    conda:
        "../envs/trackhub.yml"
    log:
        logp("trackhub/{prefix}_non_callable.log"),
    shadow:
        "minimal"
    shell:
        """
        exec 2> {log}
        sort -k1,1V {input.chrom_sizes} > sorted.chrom.sizes
        # zcat -f reads the BED whether or not it is gzipped
        zcat -f {input.callable_sites} | sort -k1,1V -k2,2n > sorted_callable_sites.bed
        bedtools complement -i sorted_callable_sites.bed -g sorted.chrom.sizes > non_callable.bed
        bedSort non_callable.bed non_callable.sorted.bed
        bedToBigBed non_callable.sorted.bed sorted.chrom.sizes {output.bb}
        """
