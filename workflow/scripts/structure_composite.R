#' Composite ancestry figure: everything you need to choose K.
#'
#' Panels are assembled according to what the run actually produced:
#'
#'   always        ADMIXTURE barplot at best K, CV error, Evanno dK, PCA,
#'                 ADMIXTURE Fst between inferred populations
#'   coordinates   ADMIXTURE ancestry pies on a map
#'   TESS          TESS pies on a map, TESS cross-entropy, and the
#'                 ADMIXTURE-vs-TESS correlation heatmap
#'
#' Also writes a CSV of per-sample ancestry coefficients: ADMIXTURE and TESS
#' side by side with components matched up by correlation when TESS ran,
#' ADMIXTURE alone otherwise.

suppressMessages({
  # maps must load before tidyverse: maps::map would otherwise mask purrr::map.
  library(maps)
  library(tidyverse)
  library(patchwork)
  library(scatterpie)
})

log_smk <- function() {
  if (exists("snakemake") && length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

source(file.path(snakemake@scriptdir, "basemap.R"))

# ---------------------------------------------------------------------- inputs
cv_errors_file <- snakemake@input[["cv_errors"]]
eigenvec_file <- snakemake@input[["eigenvec"]]
eigenval_file <- snakemake@input[["eigenval"]]
fam_file <- snakemake@input[["fam"]]
coords_file <- snakemake@input[["coords"]]

admixture_path <- snakemake@params[["admixture_path"]]
plot_seed <- snakemake@params[["plot_seed"]]
use_coords <- as.logical(snakemake@params[["use_coords"]])
run_tess <- as.logical(snakemake@params[["run_tess"]])
kmax <- as.integer(snakemake@params[["kmax"]])

output_pdf <- snakemake@output[["pdf"]]
output_merged <- snakemake@output[["merged"]]

fam <- read_delim(fam_file, col_names = FALSE, delim = " ", show_col_types = FALSE)
sample_ids <- fam$X2

coords <- NULL
if (use_coords) {
  coords <- read_tsv(coords_file,
                     col_names = c("Sample", "Longitude", "Latitude"),
                     show_col_types = FALSE)
}

cv_scores <- read_csv(cv_errors_file, show_col_types = FALSE)
admix_best_k <- unique(cv_scores$best_k)
admix_min_k <- unique(cv_scores$min_k)
admix_evanno_k <- unique(cv_scores$evanno_k)

# ------------------------------------------------------------- ADMIXTURE Q data
#' Read the Q matrix for one K from the plotted seed.
read_q <- function(k) {
  path <- paste0(admixture_path, ".", k, ".", plot_seed, ".Q")
  q <- read_delim(path, col_names = paste0("K", seq_len(k)), delim = " ",
                  show_col_types = FALSE)
  q$Sample <- sample_ids
  q$samp_num <- seq_len(nrow(q))
  q$K_value <- k
  q
}

all_admix <- map_dfr(seq_len(kmax), read_q)
admix_best <- all_admix %>% filter(K_value == admix_best_k)
k_columns <- paste0("K", seq_len(admix_best_k))

#' Sample ordering that groups samples by their dominant component, then sorts
#' within group by coefficient. Applied consistently across panels.
assignment_order <- function(q, kcols) {
  q %>%
    mutate(likely_assignment = kcols[apply(select(., all_of(kcols)), 1, which.max)]) %>%
    arrange(likely_assignment, across(all_of(kcols), desc)) %>%
    mutate(sample_order = row_number()) %>%
    select(Sample, likely_assignment, sample_order)
}

order_admix <- assignment_order(admix_best, k_columns)
admix_best <- left_join(admix_best, order_admix, by = "Sample")

admix_long <- admix_best %>%
  pivot_longer(all_of(k_columns), names_to = "component", values_to = "prob")

pal <- ancestry_palette(admix_best_k)

# --------------------------------------------------------------- always panels
p_bar <- admix_long %>%
  ggplot(aes(x = factor(sample_order), y = prob, fill = component)) +
  geom_col() +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  labs(x = NULL, y = "Ancestry",
       title = paste0("ADMIXTURE, K = ", admix_best_k, " (seed ", plot_seed, ")")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, size = 4),
    legend.position = "none"
  )

#' Dashed marker plus label at a chosen K. Silently omitted when that K is
#' undefined — Evanno's delta K needs at least three values of K, so short K
#' ranges legitimately have none.
k_marker <- function(k, label = NULL, colour = "black", y = NULL) {
  if (length(k) == 0 || all(is.na(k))) return(NULL)
  layers <- list(geom_vline(xintercept = k, color = colour, linetype = "dashed"))
  if (!is.null(label) && !is.null(y)) {
    layers <- c(layers, list(
      annotate("text", x = k - 0.2, y = y, label = label, size = 3,
               color = colour, angle = 90)
    ))
  }
  layers
}

p_cv <- cv_scores %>%
  group_by(kval) %>%
  mutate(mean = mean(cv), sd = sd(cv)) %>%
  ungroup() %>%
  ggplot(aes(x = kval, y = mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, color = "grey") +
  geom_point() +
  geom_line(aes(group = 1)) +
  k_marker(admix_min_k, "Minimum K", "red", mean(cv_scores$cv)) +
  k_marker(admix_best_k, "Best K", "blue", mean(cv_scores$cv)) +
  k_marker(admix_evanno_k, "\u0394K", "purple", max(cv_scores$cv)) +
  scale_x_continuous(breaks = unique(cv_scores$kval)) +
  labs(title = "ADMIXTURE cross-validation", x = "K", y = "CV error") +
  theme_bw()

p_deltak <- cv_scores %>%
  arrange(seed, kval) %>%
  mutate(secondrate = abs((lead(cv) - cv) - (cv - lag(cv)))) %>%
  group_by(kval) %>%
  mutate(delta_k = mean(secondrate) / sd(cv)) %>%
  ungroup() %>%
  select(kval, delta_k) %>%
  distinct() %>%
  ggplot(aes(x = kval, y = delta_k)) +
  k_marker(admix_evanno_k, colour = "purple") +
  geom_point() +
  geom_line(aes(group = 1)) +
  scale_x_continuous(breaks = unique(cv_scores$kval)) +
  labs(title = "Evanno ΔK", x = "K", y = "ΔK") +
  theme_bw()

# PCA, coloured by the ADMIXTURE assignment so the two views can be compared.
eigenvec <- read_tsv(eigenvec_file, col_names = FALSE, skip = 1, show_col_types = FALSE)
colnames(eigenvec) <- c("Sample", paste0("PC", seq_len(ncol(eigenvec) - 1)))
eigenval <- read_tsv(eigenval_file, col_names = FALSE, show_col_types = FALSE)
eigenvec <- left_join(eigenvec, order_admix, by = "Sample")

pc_lab <- function(i) {
  paste0("PC", i, " (", round(100 * eigenval[i, 1] / sum(eigenval), 2), "%)")
}

p_pca <- ggplot(eigenvec, aes(x = PC1, y = PC2, color = likely_assignment)) +
  geom_point() +
  scale_color_manual(values = pal, na.value = "grey70") +
  labs(title = "PCA, coloured by ADMIXTURE assignment",
       x = pc_lab(1), y = pc_lab(2), color = NULL) +
  theme_bw()

# Fst between inferred populations, as reported in the ADMIXTURE logs.
extract_fst <- function(file) {
  lines <- readLines(file, warn = FALSE)
  k <- as.numeric(str_match(file, "_K([0-9]+)_")[, 2])
  start <- grep("Fst divergences between estimated populations:", lines)
  if (length(start) == 0 || k < 2) return(NULL)
  block <- lines[(start + 2):(start + k)]
  rows <- block %>%
    str_split("\t+") %>%
    purrr::map(~ .x[.x != ""]) %>%
    purrr::map(~ suppressWarnings(as.numeric(.x[-1])))
  tibble(Fst = unlist(rows), K = k) %>% filter(!is.na(Fst))
}

fst_df <- map_dfr(snakemake@input[["admixture_logs"]], extract_fst)

p_fst <- if (nrow(fst_df) > 0) {
  ggplot(fst_df, aes(x = K, y = Fst)) +
    geom_boxplot(aes(group = K)) +
    geom_point(size = 0.8) +
    scale_x_continuous(breaks = unique(fst_df$K)) +
    labs(title = "Fst between ADMIXTURE populations") +
    theme_bw()
} else {
  ggplot() + geom_blank() + labs(title = "No Fst reported (K = 1 only)") + theme_bw()
}

panels <- list(p_bar, p_cv, p_deltak, p_pca, p_fst)

# ----------------------------------------------------------- coordinate panels
if (use_coords) {
  admix_map <- inner_join(coords, admix_best, by = "Sample")
  bm <- make_basemap(admix_map)
  p_admix_map <- ggplot() +
    bm$layers +
    geom_scatterpie(data = admix_map,
                    aes(x = Longitude, y = Latitude, r = bm$r),
                    cols = k_columns, color = NA) +
    scale_fill_manual(values = pal) +
    labs(title = paste0("ADMIXTURE ancestry, K = ", admix_best_k),
         x = NULL, y = NULL, fill = NULL) +
    theme_minimal()
  panels <- c(panels, list(p_admix_map))
}

# ----------------------------------------------------------------- TESS panels
merged_out <- admix_best %>% select(Sample, samp_num, all_of(k_columns), sample_order)

if (run_tess) {
  tess_in <- read_csv(snakemake@input[["tess"]], show_col_types = FALSE)
  names(tess_in) <- c("Sample", "pop_assignment", "kval", "qvalue",
                      "total_K", "best_k", "min_k")
  tess_best_k <- unique(tess_in$best_k)
  tess_min_k <- unique(tess_in$min_k)

  tess_wide_all <- tess_in %>%
    pivot_wider(id_cols = c(Sample, pop_assignment, total_K),
                names_from = kval, values_from = qvalue) %>%
    rename(K_value = total_K) %>%
    mutate(K_value = as.numeric(K_value))

  tess_best <- tess_wide_all %>% filter(K_value == tess_best_k)
  tess_k_cols <- paste0("K", seq_len(tess_best_k))

  if (use_coords) {
    tess_map <- inner_join(coords, tess_best, by = "Sample")
    bm_t <- make_basemap(tess_map)
    p_tess_map <- ggplot() +
      bm_t$layers +
      geom_scatterpie(data = tess_map,
                      aes(x = Longitude, y = Latitude, r = bm_t$r),
                      cols = tess_k_cols, color = NA) +
      scale_fill_manual(values = ancestry_palette(tess_best_k)) +
      labs(title = paste0("TESS ancestry, K = ", tess_best_k),
           x = NULL, y = NULL, fill = NULL) +
      theme_minimal()
    panels <- c(panels, list(p_tess_map))
  }

  xval <- read_csv(snakemake@input[["xval"]], show_col_types = FALSE)
  p_xent <- ggplot(xval, aes(x = K, y = xent)) +
    geom_point() +
    geom_line() +
    k_marker(tess_min_k, "Minimum K", "red", mean(xval$xent)) +
    k_marker(tess_best_k, "Best K", "blue", mean(xval$xent)) +
    scale_x_continuous(breaks = xval$K) +
    labs(title = "TESS cross-entropy", x = "K", y = "Cross-entropy") +
    theme_bw()
  panels <- c(panels, list(p_xent))

  # ADMIXTURE and TESS label their components arbitrarily, so before the two can
  # be compared the components have to be matched up. For each K, correlate every
  # ADMIXTURE component against every TESS component and take the best match.
  # If two ADMIXTURE components claim the same TESS component the matching is
  # ambiguous and that K is reported as NA rather than guessed at.
  match_k <- function(k) {
    a <- all_admix %>%
      filter(K_value == k) %>%
      select(Sample, all_of(paste0("K", seq_len(k)))) %>%
      rename_with(~ paste0(.x, ".admix"), -Sample)
    t <- tess_wide_all %>%
      filter(K_value == k) %>%
      select(Sample, any_of(paste0("K", seq_len(k)))) %>%
      rename_with(~ paste0(.x, ".tess"), -Sample)

    joined <- inner_join(a, t, by = "Sample") %>% select(where(~ !any(is.na(.))))
    ax <- joined %>% select(ends_with(".admix"))
    tx <- joined %>% select(ends_with(".tess"))
    if (nrow(joined) == 0 || ncol(ax) == 0 || ncol(tx) == 0) return(NULL)

    cor_matrix <- suppressWarnings(cor(ax, tx))
    if (k == 1) cor_matrix[1] <- 1

    best_match <- colnames(tx)[apply(cor_matrix, 1, which.max)]
    ambiguous <- anyDuplicated(best_match) > 0

    diag_cor <- if (ambiguous) rep(NA_real_, ncol(ax)) else
      vapply(seq_len(ncol(ax)), function(i) cor_matrix[i, best_match[i]], numeric(1))

    list(
      summary = tibble(K_value = k, component = paste0("K", seq_along(diag_cor)),
                       correlation = diag_cor),
      samples = if (ambiguous) NULL else joined %>% mutate(K_value = k)
    )
  }

  matches <- map(seq_len(min(kmax, max(tess_wide_all$K_value))), match_k)
  matches <- matches[!map_lgl(matches, is.null)]

  heatmap_data <- map_dfr(matches, "summary") %>%
    mutate(component = factor(as.numeric(str_remove(component, "K"))))

  p_heat <- ggplot(heatmap_data,
                   aes(x = component, y = as.factor(K_value), fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = ifelse(is.na(correlation), "", sprintf("%.2f", correlation))),
              size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, na.value = "grey85", limits = c(-1, 1)) +
    labs(title = "ADMIXTURE vs TESS component correlation",
         x = "Component", y = "K", fill = "r") +
    theme_minimal()
  panels <- c(panels, list(p_heat))

  matched_samples <- map_dfr(matches, "samples")
  if (nrow(matched_samples) > 0) merged_out <- matched_samples
}

write_csv(merged_out, output_merged)

# ------------------------------------------------------------------- assemble
n <- length(panels)
combined <- wrap_plots(panels, ncol = 2) +
  plot_annotation(
    tag_levels = "A",
    title = paste0(basename(admixture_path), " — ancestry summary"),
    subtitle = paste0(
      "ADMIXTURE best K = ", admix_best_k,
      "; coordinates: ", ifelse(use_coords, "yes", "no"),
      "; TESS: ", ifelse(run_tess, "yes", "no")
    )
  )

ggsave(output_pdf, combined, width = 11, height = 3.2 * ceiling(n / 2),
       units = "in", limitsize = FALSE)
