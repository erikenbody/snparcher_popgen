#' Assign samples to populations from ancestry coefficients.
#'
#' This is the checkpoint that everything downstream depends on: it writes one
#' population_<k>.txt per population, and those files determine which ROH, pi,
#' F_ROH, GONE and clam jobs exist.
#'
#' Which coefficients are used is set by `method`:
#'   admixture  the plotted ADMIXTURE replicate at K
#'   tess       TESS at K (requires coordinates)
#' K is `manual_k_assignment` when set, otherwise the automatically chosen best K
#' for that method.
#'
#' A sample whose top coefficient is below qvalue_threshold is treated as
#' admixed: it lands in population_NA.txt, which no downstream rule reads, so it
#' is excluded from every per-population statistic. That exclusion is recorded in
#' the populations-done file and in the all-qvalues table.
#'
#' Panels: ancestry barplot always; ancestry pies on a map when coordinates are
#' available; the same pair for the method not chosen, for comparison, when both
#' are available.

suppressMessages({
  # maps must load before tidyverse: maps::map would otherwise mask purrr::map.
  library(maps)
  library(tidyverse)
  library(patchwork)
  library(scatterpie)
  library(ggrepel)
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
fam_file <- snakemake@input[["fam"]]
coords_file <- snakemake@input[["coords"]]
cv_errors_file <- snakemake@input[["cv_errors"]]

admixture_path <- snakemake@params[["admixture_path"]]
plot_seed <- snakemake@params[["plot_seed"]]
method <- snakemake@params[["method"]]
manual_k <- snakemake@params[["manual_k_assignment"]]
threshold <- as.numeric(snakemake@params[["threshold_qval"]])
use_coords <- as.logical(snakemake@params[["use_coords"]])
run_tess <- as.logical(snakemake@params[["run_tess"]])
project_id <- snakemake@params[["project_id"]]

out_pdf <- snakemake@output[["plot"]]
out_dir <- snakemake@output[["outdir"]]
out_qvalues <- snakemake@output[["all_qvalues"]]
out_done <- snakemake@output[["pops_done"]]

fam <- read_delim(fam_file, col_names = FALSE, delim = " ", show_col_types = FALSE)
sample_ids <- fam$X2

coords <- NULL
if (use_coords) {
  coords <- read_tsv(coords_file,
                     col_names = c("Sample", "Longitude", "Latitude"),
                     show_col_types = FALSE)
}

# ------------------------------------------------------------------ choose K
cv_scores <- read_csv(cv_errors_file, show_col_types = FALSE)
admix_best_k <- unique(cv_scores$best_k)

tess_wide <- NULL
tess_best_k <- NULL
if (run_tess) {
  tess_in <- read_csv(snakemake@input[["tess"]], show_col_types = FALSE)
  names(tess_in) <- c("Sample", "pop_assignment", "kval", "qvalue",
                      "total_K", "best_k", "min_k")
  tess_best_k <- unique(tess_in$best_k)
  tess_wide <- tess_in %>%
    mutate(total_K = as.numeric(total_K)) %>%
    pivot_wider(id_cols = c(Sample, total_K), names_from = kval, values_from = qvalue)
}

auto_k <- if (method == "tess") tess_best_k else admix_best_k

set_k <- if (is.null(manual_k) || identical(as.character(manual_k), "NA")) {
  message("manual_k_assignment not set — using the automatically chosen K = ", auto_k,
          " (", method, ")")
  as.integer(auto_k)
} else {
  as.integer(manual_k)
}

message("assigning populations with method=", method, ", K=", set_k,
        ", qvalue_threshold=", threshold)

k_columns <- paste0("K", seq_len(set_k))

# --------------------------------------------------------- coefficient tables
read_admixture_q <- function(k) {
  path <- paste0(admixture_path, ".", k, ".", plot_seed, ".Q")
  q <- read_delim(path, col_names = paste0("K", seq_len(k)), delim = " ",
                  show_col_types = FALSE)
  q$Sample <- sample_ids
  q$samp_num <- seq_len(nrow(q))
  q
}

admix_q <- read_admixture_q(set_k)

tess_q <- NULL
if (run_tess) {
  tess_q <- tess_wide %>%
    filter(total_K == set_k) %>%
    select(Sample, any_of(k_columns))
  if (nrow(tess_q) == 0) {
    stop("TESS produced no results for K = ", set_k,
         ". Either widen tess_kvals or pick a K that TESS ran.")
  }
  tess_q$samp_num <- match(tess_q$Sample, sample_ids)
}

#' Add the dominant component, its coefficient, and a plotting order.
annotate_assignment <- function(q, kcols) {
  mat <- as.matrix(select(q, all_of(kcols)))
  q %>%
    mutate(
      top_component = kcols[apply(mat, 1, which.max)],
      top_qvalue = apply(mat, 1, max)
    ) %>%
    arrange(top_component, desc(top_qvalue)) %>%
    mutate(sample_order = row_number())
}

admix_assign <- annotate_assignment(admix_q, k_columns)
tess_assign <- if (run_tess) annotate_assignment(tess_q, k_columns) else NULL

chosen <- if (method == "tess") tess_assign else admix_assign

pal <- ancestry_palette(set_k)

# ------------------------------------------------------------------- plotting
barplot_panel <- function(assign, kcols, title) {
  assign %>%
    pivot_longer(all_of(kcols), names_to = "component", values_to = "prob") %>%
    ggplot(aes(x = factor(sample_order), y = prob, fill = component)) +
    geom_col() +
    scale_fill_manual(values = pal) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 1)) +
    labs(x = NULL, y = "Ancestry", title = title) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, size = 4),
      legend.position = "none"
    )
}

map_panels <- function(assign, kcols, title) {
  dat <- inner_join(coords, assign, by = "Sample")
  if (nrow(dat) == 0) return(NULL)
  bm <- make_basemap(dat)
  base <- ggplot() +
    bm$layers +
    geom_scatterpie(data = dat, aes(x = Longitude, y = Latitude, r = bm$r),
                    cols = kcols, color = NA) +
    scale_fill_manual(values = pal) +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    theme_minimal()
  labelled <- base +
    geom_text_repel(data = dat,
                    aes(x = Longitude, y = Latitude, label = sample_order),
                    size = 2, max.overlaps = 10, force = 2, nudge_y = 0.3,
                    segment.color = "black") +
    labs(title = paste0(title, " (sample labels)"))
  list(base, labelled)
}

panels <- list()

primary_title <- paste0(project_id, ": ", toupper(method), " K = ", set_k)
if (use_coords) {
  panels <- c(panels, map_panels(chosen, k_columns, primary_title))
}
panels <- c(panels, list(barplot_panel(chosen, k_columns, primary_title)))

# The method that was not used for assignment, shown for comparison.
if (run_tess) {
  other <- if (method == "tess") admix_assign else tess_assign
  other_name <- if (method == "tess") "ADMIXTURE" else "TESS"
  other_title <- paste0(other_name, " K = ", set_k, " (not used for assignment)")
  if (use_coords) {
    panels <- c(panels, map_panels(other, k_columns, other_title)[1])
  }
  panels <- c(panels, list(barplot_panel(other, k_columns, other_title)))
}

panels <- panels[!map_lgl(panels, is.null)]

combined <- wrap_plots(panels, ncol = 2) + plot_annotation(tag_levels = "A")
ggsave(out_pdf, combined, width = 11,
       height = 3.4 * ceiling(length(panels) / 2), units = "in",
       limitsize = FALSE)

# ------------------------------------------------------------------- assignment
assigned <- chosen %>%
  mutate(population = if_else(top_qvalue < threshold, NA_character_,
                              str_remove(top_component, "^K")))

n_admixed <- sum(is.na(assigned$population))
message(nrow(assigned) - n_admixed, " samples assigned; ",
        n_admixed, " below the qvalue threshold and excluded")

write_tsv(assigned, out_qvalues)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

assigned %>%
  filter(!is.na(population)) %>%
  group_split(population) %>%
  walk(function(df) {
    pop <- unique(df$population)
    write_tsv(df %>% select(Sample),
              file.path(out_dir, paste0("population_", pop, ".txt")),
              col_names = FALSE)
  })

excluded <- assigned %>% filter(is.na(population)) %>% pull(Sample)
if (length(excluded) > 0) {
  write_tsv(tibble(Sample = excluded),
            file.path(out_dir, "population_NA.txt"), col_names = FALSE)
}

# A human-readable record of what was decided, and the marker every downstream
# rule waits on.
con <- file(out_done, open = "wt")
writeLines(c(
  paste0("project: ", project_id),
  paste0("method: ", method),
  paste0("K: ", set_k, if (is.null(manual_k) || identical(as.character(manual_k), "NA")) " (automatic)" else " (manual)"),
  paste0("qvalue_threshold: ", threshold),
  paste0("samples assigned: ", nrow(assigned) - n_admixed),
  paste0("samples excluded as admixed: ", n_admixed),
  "",
  "excluded samples:"
), con)
writeLines(excluded, con)
writeLines("", con)
for (pop in sort(unique(na.omit(assigned$population)))) {
  writeLines(paste0("@", pop), con)
  writeLines(assigned$Sample[!is.na(assigned$population) & assigned$population == pop], con)
  writeLines("", con)
}
close(con)
