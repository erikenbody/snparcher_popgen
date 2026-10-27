## TESS3 ancestry over a range of K, with cross-entropy K selection.
##
## This calls tess3r directly. It replaces algatr::tess_ktest(K_selection =
## "auto") and algatr::pops_helper, and wingen::vcf_to_dosage, reproducing their
## logic line for line (algatr R/TESS.R, commit c0841bd) so that results match
## the earlier algatr-based runs. The output files and their columns are
## unchanged (bar the peak-RAM table, superseded by the rule's benchmark): pie_plot.R, structure_composite.R and ancestry_assignment.R read
## them as before.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(tess3r)
})

if (exists("snakemake") && length(snakemake@log) != 0) {
  log <- file(snakemake@log[1][[1]], open = "wt")
  sink(log, append = TRUE)
  sink(log, append = TRUE, type = "message")
}

species <- snakemake@params[["species"]]
kvals <- eval(parse(text = snakemake@params[["kvals"]]))
set.seed(as.integer(snakemake@params[["seed"]]))

# --- input ------------------------------------------------------------------
# Rows are individuals, columns SNPs, values the ALT allele count (0/1/2/NA):
# what vcfR2genlight + as.matrix gives, as in wingen/algatr vcf_to_dosage.
read_dosage <- function(vcf_path) {
  vcf <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
  as.matrix(vcfR::vcfR2genlight(vcf))
}

# Coordinates in genotype order. Every genotyped sample must have coordinates
# (prepare already dropped the ones that do not); extra coordinate rows are
# ignored.
read_coords <- function(coords_path, sample_ids) {
  coords <- read_tsv(coords_path, col_names = c("INDV", "x", "y"),
                     show_col_types = FALSE)
  missing <- setdiff(sample_ids, coords$INDV)
  if (length(missing)) {
    stop("Missing coordinate data for: ", paste(missing, collapse = ", "))
  }
  coords <- coords[match(sample_ids, coords$INDV), ]
  m <- as.matrix(coords[, c("x", "y")])
  colnames(m) <- c("x", "y")
  m
}

gen <- read_dosage(snakemake@input[["vcf"]])
coords <- read_coords(snakemake@input[["coords"]], rownames(gen))

# --- TESS -------------------------------------------------------------------
# algatr::bestK: the smallest K whose cross-entropy slope is within the upper
# quartile of slopes. Indexes into K, so it assumes kvals is 1..n — as did
# algatr.
best_k <- function(tess3_obj, kvals) {
  ce <- list()
  for (k in kvals) ce[[k]] <- tess3_obj[[k]]$crossentropy
  ce.K <- c()
  for (k in kvals) ce.K[k] <- min(ce[[k]])
  diff <- ce.K[-1] - ce.K[-max(kvals)]
  slope <- exp(-diff) - 1
  min(which(slope <= quantile(slope)[4]))
}

# algatr::pops_helper: Q matrix for K with columns K1..KK and the argmax.
pops_helper <- function(gen, tess3_obj, K) {
  qmat <- as.data.frame(tess3r::qmatrix(tess3_obj, K = K))
  colnames(qmat) <- gsub("V", "K", colnames(qmat))
  data.frame(individual = rownames(gen), qmat) %>%
    rowwise() %>%
    mutate(pop_assignment = which.max(c_across(-individual))) %>%
    ungroup()
}

tess3_obj <- tess3r::tess3(X = gen, coord = coords, K = kvals,
                           method = "projected.ls", lambda = 1, ploidy = 2)
bestK <- if (length(kvals) > 1) best_k(tess3_obj, kvals) else kvals

ce.K <- c()
for (k in kvals) ce.K[k] <- min(tess3_obj[[k]]$crossentropy)
minK <- kvals[which.min(ce.K[kvals])]

# --- export -----------------------------------------------------------------
output_qmatrix <- snakemake@output[["qmatrix"]]
if (length(kvals) > 1) {
  qvals <- lapply(seq_along(kvals), function(x) {
    pops_helper(gen, tess3_obj, K = x) %>%
      pivot_longer(cols = "K1":paste0("K", x), names_to = "kval",
                   values_to = "qvalue") %>%
      mutate(total_K = paste0(x))
  }) %>%
    bind_rows()
  qvals$best_k <- bestK
  qvals$min_k <- minK
} else {
  qvals <- pops_helper(gen, tess3_obj, K = kvals)
}
write_csv(qvals, output_qmatrix, col_names = TRUE)

xval <- data.frame(
  K = sapply(kvals, function(k) tess3_obj[[k]]$K),
  xent = sapply(kvals, function(k) tess3_obj[[k]]$crossentropy),
  rmse = sapply(kvals, function(k) tess3_obj[[k]]$rmse)
)
write_csv(xval, snakemake@output[["xval"]], col_names = TRUE)

p <- ggplot(xval, aes(x = K, y = xent)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ggtitle(species) +
  ylab("Cross-entropy value") +
  scale_x_continuous(breaks = kvals)
if (length(kvals) > 1) {
  p <- p +
    geom_vline(xintercept = minK, color = "red", linetype = "dashed") +
    annotate("text", x = minK - 0.2, y = mean(xval$xent), label = "Minimum K",
             size = 4, color = "red", angle = 90) +
    geom_vline(xintercept = bestK, color = "blue", linetype = "dashed") +
    annotate("text", x = bestK - 0.2, y = mean(xval$xent), label = "Best K",
             size = 4, color = "blue", angle = 90)
}
ggsave(snakemake@output[["xval_plot"]], p, width = 20, height = 15, units = "cm")
