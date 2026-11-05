# Per-sample genotype counts in sampled windows, for the scaffold QC (rule
# sex_id_extract).
#
# Input:  bcftools query -f '%CHROM\t%POS[\t%GT:%DP]\n', restricted (-R) to the
#         windows.
# -v step=<window spacing>: each window sits inside one step-sized block, and the
#         block's 0-based start is the window's id.
# Output: one line per window and sample with any data: chrom, block, sample
#         index (0-based, in `bcftools query -l` order), called genotypes,
#         heterozygous genotypes, summed DP, genotypes with a DP. A window can come
#         out on more than one line (its sites reached in separate pieces); the
#         reader sums them.
# Written for speed (substr/index, no regex or split per genotype: ~10x faster);
# genotypes that are not one character / separator / one character (multi-digit
# alleles, haploid calls) take the slower general path.
BEGIN { FS = OFS = "\t" }

function flush(   i) {
    if (key == "") return
    for (i = 3; i <= nf; i++)
        if (called[i] || dpn[i])
            print chrom, block, i - 3, called[i] + 0, het[i] + 0, dps[i] + 0, dpn[i] + 0
    delete called; delete het; delete dps; delete dpn
}

{
    b = int(($2 - 1) / step) * step
    k = $1 SUBSEP b
    if (k != key) { flush(); key = k; chrom = $1; block = b; nf = NF }
    for (i = 3; i <= NF; i++) {
        g = $i
        p = index(g, ":")
        if (p == 4) {
            a1 = substr(g, 1, 1)
            if (a1 != ".") {
                called[i]++
                a2 = substr(g, 3, 1)
                if (a1 != a2 && a2 != ".") het[i]++
            }
        } else {
            gt = p ? substr(g, 1, p - 1) : g
            if (gt !~ /^\.([\/|]\.)?$/) {
                called[i]++
                if (split(gt, a, /[\/|]/) == 2 && a[1] != a[2] && a[1] != "." && a[2] != ".") het[i]++
            }
        }
        if (p) {
            v = substr(g, p + 1)
            if (v != ".") { dps[i] += v; dpn[i]++ }
        }
    }
}

END { flush() }
