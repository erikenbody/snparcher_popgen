"""Rewrite a plink .bim so chromosome names are integers.

ADMIXTURE refuses non-numeric chromosome names. The mapping is taken from the
reference .fai, in file order, so contig N becomes N.
"""

# Send stdout and stderr to the rule's log file. Snakemake does this for R
# scripts (via sink) but not for python ones, which otherwise leaves the log
# empty and the traceback buried in the master log.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

def generate_mapping(input_file, bim_file, output_file):

    # fai column 1 is the contig name; number contigs 1..n in fai order.
    conversion_dict = {}
    with open(input_file, 'r') as f:
        for i, line in enumerate(f, start=1):
            contig = line.split('\t')[0]
            conversion_dict[contig] = str(i)

    # read bim file and replace the scaffold names with numbering 1:n (n = number of scaffolds)
    updated_lines = []
    with open(bim_file, 'r') as f:
        for line in f:
            elements = line.strip().split('\t')
            scaffold = elements[0]
            if scaffold in conversion_dict:
                elements[0] = conversion_dict[scaffold]
            updated_lines.append('\t'.join(elements))

    with open(output_file, 'w') as f:
        for line in updated_lines:
            f.write(line + '\n')

input_file = snakemake.input.fai
bim_file = snakemake.input.bim
output_file = snakemake.output.bim
generate_mapping(input_file, bim_file, output_file)
