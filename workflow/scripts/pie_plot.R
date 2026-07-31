#' TESS ancestry proportions as pie charts at each sampling locality.
#'
#' The map extent follows the coordinates (see basemap.R), so this works for any
#' region rather than a fixed one.

suppressMessages({
  # maps must load before tidyverse: maps::map would otherwise mask purrr::map.
  library(maps)
  library(tidyverse)
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

species <- snakemake@params[["species"]]
tess_path <- snakemake@input[["tess"]]
coords_path <- snakemake@input[["coords"]]
output_path <- snakemake@output[["plot"]]

df_tess <- read_csv(tess_path, show_col_types = FALSE)
df_coord <- read_tsv(coords_path,
                     col_names = c("Sample", "Longitude", "Latitude"),
                     show_col_types = FALSE)

df_tess <- df_tess %>%
  left_join(df_coord, by = c("individual" = "Sample")) %>%
  filter(total_K == best_k)   # only the selected K

wide_data <- df_tess %>%
  pivot_wider(names_from = kval, values_from = qvalue)

kvals <- unique(df_tess$kval)
bm <- make_basemap(wide_data)

ggplot() +
  bm$layers +
  geom_scatterpie(
    data = wide_data,
    aes(x = Longitude, y = Latitude, group = individual, r = bm$r),
    cols = kvals, color = NA
  ) +
  scale_fill_manual(values = ancestry_palette(length(kvals)), name = "Population") +
  theme_void() +
  ggtitle(paste0(species, ": TESS ancestry (K = ", unique(df_tess$best_k), ")"))

ggsave(output_path, width = 8, height = 6)
