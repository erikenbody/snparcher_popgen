"""Raw VCF -> clean biallelic SNP set.

Only included when start_from is "raw". If you already have
results/{refGenome}/{prefix}_clean_snps.vcf.gz (for example straight out of
snpArcher), set start_from: clean_snps and none of these rules exist.

Two ways to remove unreliable sites, chosen by whether callable_sites_bed is set:
  unset  Depth outliers are found in the VCF itself: sites whose total depth is
         more than cov_threshold_stdev SDs from the mean are dropped. Needs no
         extra input files.
  set    Sites are restricted to the callable-sites BED, optionally also
         dropping short contigs (min_contig_size).
"""


if START_FROM == "raw":

    rule basic_filter:
        """Drop excluded samples and hard-fail sites, then recompute AF.

        -a -U trims alleles no longer seen and drops sites that become
        uncalled once samples are removed; +fill-tags recomputes AF/F_MISSING
        on the reduced sample set so the thresholds downstream mean what they
        say.
        """
        input:
            vcf="results/{refGenome}/{prefix}_raw.vcf.gz",
            exclude=[config["exclude_samples"]] if config["exclude_samples"] else [],
        output:
            vcf=temp("results/{refGenome}/filter/{prefix}_sample_filtered.vcf.gz"),
            idx=temp("results/{refGenome}/filter/{prefix}_sample_filtered.vcf.gz.csi"),
        params:
            drop=lambda wc, input: f"-S ^{input.exclude[0]}" if input.exclude else "",
        conda:
            "../envs/filter.yml"
        log:
            "logs/{refGenome}/filter/{prefix}_basic_filter.log",
        shell:
            """
            exec 2> {log}
            bcftools view {params.drop} -f .,PASS {input.vcf} -a -U -O u \
              | bcftools +fill-tags -Ou \
              | bcftools view -m2 -e 'AF==0 | ref="N" | ALT="."' -O z -o {output.vcf}
            bcftools index {output.vcf}
            """

    if not USE_BED:

        rule site_depth:
            input:
                vcf="results/{refGenome}/filter/{prefix}_sample_filtered.vcf.gz",
                idx="results/{refGenome}/filter/{prefix}_sample_filtered.vcf.gz.csi",
            output:
                depth=temp("results/{refGenome}/filter/{prefix}_depth.txt"),
            conda:
                "../envs/filter.yml"
            log:
                "logs/{refGenome}/filter/{prefix}_site_depth.log",
            shell:
                "bcftools query -f '%CHROM\\t%POS\\t%DP\\n' {input.vcf} > {output.depth} 2> {log}"

        rule depth_outlier_sites:
            """List sites whose total depth is an outlier; excluded below."""
            input:
                depth="results/{refGenome}/filter/{prefix}_depth.txt",
            output:
                pos="results/{refGenome}/filter/{prefix}_excluded_positions.txt",
            params:
                cov_threshold_stdev=config["cov_threshold_stdev"],
            conda:
                "../envs/filter.yml"
            log:
                "logs/{refGenome}/filter/{prefix}_depth_outliers.log",
            script:
                "../scripts/get_snp_depth_pos.py"

        rule strict_filter:
            input:
                pos="results/{refGenome}/filter/{prefix}_excluded_positions.txt",
                vcf="results/{refGenome}/filter/{prefix}_sample_filtered.vcf.gz",
                idx="results/{refGenome}/filter/{prefix}_sample_filtered.vcf.gz.csi",
            output:
                vcf=temp("results/{refGenome}/filter/{prefix}_strict.vcf.gz"),
                idx=temp("results/{refGenome}/filter/{prefix}_strict.vcf.gz.csi"),
                gzpos="results/{refGenome}/{prefix}_excluded_positions.tsv.gz",
            params:
                miss=config["missingness"],
                maf=config["maf"],
                chr_ex=config["scaffolds_to_exclude"],
            conda:
                "../envs/filter.yml"
            log:
                "logs/{refGenome}/filter/{prefix}_strict_filter.log",
            shell:
                """
                exec 2> {log}
                bgzip -c {input.pos} > {output.gzpos}
                tabix -s1 -b2 -e2 {output.gzpos}

                if [ -z "{params.chr_ex}" ]; then
                    bcftools view -T ^{output.gzpos} -m2 -M2 \
                      -e 'F_MISSING > {params.miss} | AF<{params.maf}' \
                      {input.vcf} -O z -o {output.vcf}
                else
                    bcftools view -t ^{params.chr_ex} -m2 -M2 \
                      -e 'F_MISSING > {params.miss} | AF<{params.maf}' \
                      -O u {input.vcf} \
                      | bcftools view -T ^{output.gzpos} -O z -o {output.vcf}
                fi
                bcftools index {output.vcf}
                """

    else:

        rule callable_regions:
            """Callable sites to keep, dropping contigs below min_contig_size."""
            input:
                bed=CALLABLE_BED,
                fai="results/{refGenome}/data/genome/{refGenome}.fna.fai",
            output:
                bed="results/{refGenome}/filter/{prefix}_include_sites.bed",
            params:
                min_contig_size=config["min_contig_size"],
            conda:
                "../envs/bed.yml"
            log:
                "logs/{refGenome}/filter/{prefix}_callable_regions.log",
            shell:
                """
                exec 2> {log}
                awk -v m={params.min_contig_size} 'BEGIN{{OFS="\\t"}} $2>=m {{print $1,0,$2}}' \
                  {input.fai} | bedtools sort -i - > {output.bed}.contigs
                bedtools intersect -a {input.bed} -b {output.bed}.contigs \
                  | bedtools sort -i - | bedtools merge -i - > {output.bed}
                rm -f {output.bed}.contigs
                """

        rule strict_filter:
            input:
                bed="results/{refGenome}/filter/{prefix}_include_sites.bed",
                vcf="results/{refGenome}/filter/{prefix}_sample_filtered.vcf.gz",
                idx="results/{refGenome}/filter/{prefix}_sample_filtered.vcf.gz.csi",
            output:
                vcf=temp("results/{refGenome}/filter/{prefix}_strict.vcf.gz"),
                idx=temp("results/{refGenome}/filter/{prefix}_strict.vcf.gz.csi"),
            params:
                miss=config["missingness"],
                maf=config["maf"],
                chr_ex=config["scaffolds_to_exclude"],
            conda:
                "../envs/filter.yml"
            log:
                "logs/{refGenome}/filter/{prefix}_strict_filter.log",
            shell:
                """
                exec 2> {log}
                if [ -z "{params.chr_ex}" ]; then
                    bcftools view -R {input.bed} -m2 -M2 \
                      -e 'F_MISSING > {params.miss} | AF<{params.maf}' \
                      {input.vcf} -O z -o {output.vcf}
                else
                    bcftools view -t ^{params.chr_ex} -R {input.bed} -m2 -M2 \
                      -e 'F_MISSING > {params.miss} | AF<{params.maf}' \
                      {input.vcf} -O z -o {output.vcf}
                fi
                bcftools index {output.vcf}
                """

    rule subset_snps:
        """Clean, biallelic SNPs. This is the entry point when start_from: clean_snps."""
        input:
            vcf="results/{refGenome}/filter/{prefix}_strict.vcf.gz",
            idx="results/{refGenome}/filter/{prefix}_strict.vcf.gz.csi",
        output:
            vcf="results/{refGenome}/{prefix}_clean_snps.vcf.gz",
            idx="results/{refGenome}/{prefix}_clean_snps.vcf.gz.tbi",
        conda:
            "../envs/filter.yml"
        log:
            "logs/{refGenome}/filter/{prefix}_subset_snps.log",
        shell:
            """
            exec 2> {log}
            bcftools view -v snps -e 'TYPE ~ "indel"' -O z -o {output.vcf} {input.vcf}
            bcftools index -t {output.vcf}
            """

    rule subset_indels:
        """Clean indels. Not used downstream, but cheap and usually wanted."""
        input:
            vcf="results/{refGenome}/filter/{prefix}_strict.vcf.gz",
            idx="results/{refGenome}/filter/{prefix}_strict.vcf.gz.csi",
        output:
            vcf="results/{refGenome}/{prefix}_clean_indels.vcf.gz",
            idx="results/{refGenome}/{prefix}_clean_indels.vcf.gz.tbi",
        conda:
            "../envs/filter.yml"
        log:
            "logs/{refGenome}/filter/{prefix}_subset_indels.log",
        shell:
            """
            exec 2> {log}
            bcftools view -v indels -O z -o {output.vcf} {input.vcf}
            bcftools index -t {output.vcf}
            """
