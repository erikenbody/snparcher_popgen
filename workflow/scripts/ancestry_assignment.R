#' Assign samples to populations and plot the assignment (final_admixture_plot.R
#' from CCGP-analysis).
#'
#' This is the checkpoint that everything downstream depends on: it writes one
#' population_<k>.txt per population, and those files determine which ROH, pi,
#' F_ROH, GONE and clam jobs exist.
#'
#' `method` picks the coefficients: the best-likelihood ADMIXTURE replicate at K,
#' or TESS at K (requires coordinates). K is `manual_k_assignment` when set,
#' otherwise the automatically chosen best K for that method. The PDF shows the
#' chosen method only: pies on California, the same with sample labels, and the
#' barplot (barplot alone without coordinates).
#'
#' A sample whose top coefficient is below qvalue_threshold is treated as
#' admixed: it lands in population_NA.txt, which no downstream rule reads.

suppressMessages({
  library(tidyverse)
  library(patchwork)
  library(maps)
  library(scatterpie)
  library(ggrepel)
})


#set up log file writing
log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

source(file.path(snakemake@scriptdir, "plot_longitude.R"))

#' Pies + labelled pies over the barplot, or the barplot alone without
#' coordinates.
save_combined <- function(p1, p1_labs, p2, use_coords, output_path) {
  if (use_coords) {
    layout <- "
            AABB
            AABB
            AABB
            CCCC
            "

    combined_plot <- p1 + p1_labs + p2 + plot_annotation(tag_levels = 'A') +
      plot_layout(design = layout)
    height <- 6
  } else {
    combined_plot <- p2
    height <- 2
  }
  ggplot2::ggsave(output_path, plot = combined_plot, width = 8.5, height = height, units = "in", device = cairo_pdf)
}

create_plots_tess <- function(coords_file, tess_file, fam, set_k, output_path) {

  coords <- read_tsv(coords_file, col_names = c("Sample", "Longitude", "Latitude"))
  coords$Longitude <- plot_longitude(coords$Longitude)

  # Load the Tess data
  tess.in <- read_csv(tess_file)
  names(tess.in) <- c("Sample", "pop_assignment", "kval", "qvalue", "total_K", "best_k", "min_k")

  tess <- tess.in %>% filter(total_K == set_k)
  if (nrow(tess) == 0) {
    stop("TESS produced no results for K = ", set_k,
         ". Either widen tess_kvals or pick a K that TESS ran.")
  }
  df_tess <- left_join(tess, coords, by = "Sample")

  wide_data <- df_tess %>%
    pivot_wider(names_from = kval, values_from = qvalue)
  wide_data <- wide_data %>% rowwise() %>%
    mutate(pop_name = names(select(., starts_with("K")))[which.max(c_across(starts_with("K")))])

  #get samp num from fam file
  df_fam <- read_delim(fam, col_names = FALSE)
  df_fam$samp_num <- row.names(df_fam)
  df_fam2 <- df_fam %>% select(X2, samp_num)
  wide_data <- left_join(wide_data, df_fam2, by = c("Sample" = "X2"))
  df_tess <- left_join(df_tess, df_fam2, by = c("Sample" = "X2"))

  # Plot 1: Outline of California with pie plots for each individual
  california_map <- map_data("state", region = "california")

  tess_assign <- df_tess %>%
    select(Sample, samp_num, kval, qvalue) %>%
    group_by(Sample) %>%
    mutate(likely_assignment = kval[which.max(qvalue)],
          assingment_prob = max(qvalue)) %>%
    pivot_wider(id_cols = c("likely_assignment", "Sample", "samp_num"), names_from = "kval", values_from = "qvalue") %>%
    arrange(likely_assignment, across(all_of(unique(df_tess$kval)), desc)) %>%
    ungroup() %>%
    mutate(sample_order = row_number())

  sample.order <- tess_assign %>% select(Sample, sample_order) %>% distinct(.keep_all = T)
  df_tess <- left_join(df_tess, sample.order, by = c("Sample"))
  wide_data <- left_join(wide_data, sample.order, by = c("Sample"))

  # plots -------------------------------------------------------------------

  # Plot 3: TESS Admixture Proportions
  kvals <- unique(df_tess$kval)

  p1 <- ggplot() +
    geom_polygon(data = california_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
    geom_scatterpie(data = with_pie_r(wide_data, california_map), aes(x = Longitude, y = Latitude, r = pie_r), cols = kvals, color = NA) +
    coord_fixed() +
    labs(title = "TESS Admixture Proportions", x = NULL, y = NULL, fill = "pop") +
    theme_minimal() +
    theme(legend.position = c(0.8, 0.8)) +
    scale_fill_manual(values = selected_colors)

  p1_labs <- p1 +  geom_text_repel(data = wide_data, aes(x = Longitude, y = Latitude, label = sample_order),
                    size = 2, color = "black",
                    max.overlaps = 10,
                    force = 2,
                    nudge_y = 0.3,
                    segment.length = 0.8,
                    segment.color = "black") +
                    labs(title = "Sample Labels")

   p2 <- df_tess %>%
      ggplot(aes(x = factor(sample_order), y = qvalue, fill = factor(kval), text = sample_order)) +
      geom_col(aes(fill = factor(kval)), size = 0.1) +
      theme_minimal() +
      labs(x = NULL, y = "Ancestry") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_discrete(expand = expansion(add = 1)) +
      theme(
        panel.spacing.x = unit(0.0, "lines"),
        axis.text.x = element_text(angle = 90, size = 4),
        panel.grid = element_blank(),
        strip.text.x = element_text(angle = 0, size = 6),
        legend.position = "none"
      ) +
      xlab(NULL)  + theme(legend.position = "none") +
      scale_fill_manual(values = selected_colors)

  save_combined(p1, p1_labs, p2, TRUE, output_path)

  return(tess_assign)

}

create_plots_admixture <- function(coords_file, use_coords, admixture_path, fam, set_k, output_path) {
  # Load the admixture data for the best K value
  num_columns <- set_k
  admix_k <- set_k
  admix_best_k <- set_k

  column_names <- paste0("K", 1:num_columns)  # Generate column names from K1 to K{num_columns}

  #for simplicity, just load admixture_data as the relavent Q. Easier to plot
  admixture_data <- read_delim(paste0(admixture_path, ".", admix_best_k, ".Q"),
                               col_names = column_names,
                               delim = " ")
  df_fam <- read_delim(fam, col_names = FALSE)
  df_fam$samp_num <- row.names(df_fam)

  admixture_data$Sample <- df_fam$X2
  admixture_data$samp_num <- df_fam$samp_num
  if (use_coords) {
    coords <- read_tsv(coords_file, col_names = c("Sample", "Longitude", "Latitude"))
    coords$Longitude <- plot_longitude(coords$Longitude)
    admixture_data <- inner_join(coords, admixture_data, by = "Sample")
  } else {
    admixture_data <- admixture_data %>% mutate(Longitude = NA_real_, Latitude = NA_real_)
  }

  # Plot 1: Outline of California with pie plots for each individual
  california_map <- map_data("state", region = "california")

  k_columns <- names(admixture_data)[startsWith(names(admixture_data), "K")]

  edit_admix <- admixture_data
  names(edit_admix) <- gsub("K", "pop", names(edit_admix)) #rename ancestral pops
  numpops <- grep("pop", names(edit_admix), value = T)

  cat_admx <- edit_admix %>%
      pivot_longer(names_to = "popGroup", values_to = "prob",
        cols = -c(Sample, samp_num, Longitude, Latitude)) %>%
      mutate(popGroup = gsub("pop", "K", popGroup))

  #includes pop assignment if you want to use
  cat_admx.wide <- cat_admx %>%
    select(Sample, samp_num, popGroup, prob) %>%
    group_by(Sample) %>%
    mutate(likely_assignment = popGroup[which.max(prob)],
    assingment_prob = max(prob)) %>%
    pivot_wider(id_cols = c("likely_assignment", "Sample", "samp_num"), names_from = "popGroup", values_from = "prob") %>%
    arrange(likely_assignment, across(all_of(unique(cat_admx$popGroup)), desc)) %>%
    ungroup() %>%
    mutate(sample_order = row_number())

  sample.order <- cat_admx.wide %>% select(Sample, sample_order) %>% distinct(.keep_all = T)
  cat_admx <- left_join(cat_admx, sample.order, by = c("Sample"))
  admixture_data <- left_join(admixture_data, sample.order, by = c("Sample"))

  p1 <- ggplot() +
    geom_polygon(data = california_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
    geom_scatterpie(data = with_pie_r(admixture_data, california_map), aes(x = Longitude, y = Latitude, r = pie_r), cols = k_columns, color = NA) +
    coord_fixed() +
    labs(title = "Admixture Proportions", x = NULL, y = NULL, fill = "pop") +
    theme_minimal() +
    theme(legend.position = c(0.8, 0.8)) +
    scale_fill_manual(values = selected_colors)

  p1_labs <- p1 +  geom_text_repel(data = admixture_data, aes(x = Longitude, y = Latitude, label = sample_order),
                    size = 2, color = "black",
                    max.overlaps = 10,
                    force = 2,
                    nudge_y = 0.3,
                    segment.length = 0.8,
                    segment.color = "black") +
                    labs(title = "Sample Labels")

  p2 <- cat_admx %>%
      ggplot(aes(x = factor(sample_order), y = prob, fill = factor(popGroup), text = sample_order)) +
      geom_col(aes(fill = factor(popGroup)), size = 0.1) +
      theme_minimal() +
      labs(x = NULL, y = "Ancestry") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_discrete(expand = expansion(add = 1)) +
      theme(
        panel.spacing.x = unit(0.0, "lines"),
        axis.text.x = element_text(angle = 90, size = 4),
        panel.grid = element_blank(),
        strip.text.x = element_text(angle = 0, size = 6),
        legend.position = "none"
      ) +
      xlab(NULL)  + theme(legend.position = "none") +
      scale_fill_manual(values = selected_colors)

  save_combined(p1, p1_labs, p2, use_coords, output_path)

  return(cat_admx.wide)
}

# ------------------------------------------------------------------ choose K
method <- snakemake@params[["method"]]
manual_k <- snakemake@params[["manual_k_assignment"]]
threshold <- as.numeric(snakemake@params[["threshold_qval"]])
use_coords <- as.logical(snakemake@params[["use_coords"]])
project_id <- snakemake@params[["project_id"]]
k_path <- snakemake@output[["outdir"]]

auto_k <- if (method == "tess") {
  unique(read_csv(snakemake@input[["tess"]], show_col_types = FALSE)$best_k)
} else {
  unique(read_csv(snakemake@input[["cv_errors"]], show_col_types = FALSE)$best_k)
}

automatic <- is.null(manual_k) || identical(as.character(manual_k), "NA")
set_k <- if (automatic) {
  message("manual_k_assignment not set — using the automatically chosen K = ", auto_k,
          " (", method, ")")
  as.integer(auto_k)
} else {
  as.integer(manual_k)
}

message("assigning populations with method=", method, ", K=", set_k,
        ", qvalue_threshold=", threshold)

#select color palette: calecopal::cal_palette("chaparral1"), reordered
colors <- c("#DCC27A", "#B0B9BE", "#63605F", "#985E5C", "#AEBFA8", "#F19B34")
selected_colors <- colors[c(1, 4, 5, 2, 3, 6)]
if (set_k > length(selected_colors)) {
  # The CCGP palette has six colours; beyond that, recycle rather than fail.
  selected_colors <- rep(selected_colors, length.out = set_k)
}

if (method == "tess") {
  df_out <- create_plots_tess(
    coords_file = snakemake@input[["coords"]],
    tess_file = snakemake@input[["tess"]],
    fam = snakemake@input[["fam"]],
    set_k = set_k,
    output_path = snakemake@output[["plot"]]
  )
} else {
  df_out <- create_plots_admixture(
    coords_file = snakemake@input[["coords"]],
    use_coords = use_coords,
    admixture_path = snakemake@params[["admixture_path"]],
    fam = snakemake@input[["fam"]],
    set_k = set_k,
    output_path = snakemake@output[["plot"]]
  )
}

# ------------------------------------------------------------------- assignment
df_out <- df_out %>%
  rowwise() %>%
  mutate(likely_assignment = if_else(get(likely_assignment) < threshold, NA_character_, likely_assignment)) %>%
  ungroup()

write_tsv(df_out, snakemake@output[["all_qvalues"]])

dir.create(k_path, recursive = TRUE, showWarnings = FALSE)

df_out %>%
  filter(!is.na(likely_assignment)) %>%
  mutate(likely_assignment = gsub("K","", likely_assignment)) %>%
  split(.$likely_assignment) %>%
  walk2(., names(.), ~ write_tsv(.x %>% select(Sample), paste0(k_path, "/population_", .y, ".txt"), col_names = FALSE))

excluded <- df_out %>% filter(is.na(likely_assignment)) %>% pull(Sample)
if (length(excluded) > 0) {
  write_tsv(tibble(Sample = excluded),
            file.path(k_path, "population_NA.txt"), col_names = FALSE)
}
message(nrow(df_out) - length(excluded), " samples assigned; ",
        length(excluded), " below the qvalue threshold and excluded")

# A human-readable record of what was decided, and the marker every downstream
# rule waits on.
assigned <- df_out %>% mutate(population = gsub("K", "", likely_assignment))
con <- file(snakemake@output[["pops_done"]], open = "wt")
writeLines(c(
  paste0("project: ", project_id),
  paste0("method: ", method),
  paste0("K: ", set_k, if (automatic) " (automatic)" else " (manual)"),
  paste0("qvalue_threshold: ", threshold),
  paste0("samples assigned: ", nrow(df_out) - length(excluded)),
  paste0("samples excluded as admixed: ", length(excluded)),
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
