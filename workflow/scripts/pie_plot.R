library(tidyverse)
library(scatterpie)

#set up log file writing
log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

species <- snakemake@params[[1]]
input_path <- snakemake@input[[1]]
input_path2 <- snakemake@input[[2]]
output_path <- snakemake@output[[1]]

df_tess <- read_csv(input_path)
df_coord <- read_table(input_path2, col_names = F)
source(file.path(snakemake@scriptdir, "plot_longitude.R"))
df_coord$X2 <- plot_longitude(df_coord$X2)

df_tess <- left_join(df_tess, df_coord, by = c("individual"="X1")) %>%
  rename(lon = X2, lat = X3) %>%
  filter(total_K == best_k) #return only best K

ca_map <- map_data("state") %>%
  filter(region == "california")

wide_data <- df_tess %>%
  pivot_wider(names_from = kval, values_from = qvalue)

kvals <- unique(df_tess$kval)

ggplot() +
  geom_polygon(data = ca_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
  geom_scatterpie(data = wide_data, aes(x = lon, y = lat, group = individual, r = 0.15),
                  cols = kvals, color = NA) +
  coord_fixed(ratio = 1) +
  theme_void() +
  scale_fill_discrete(name = "Population") +
  ggtitle(species)

ggsave(output_path, width = 8, height = 6)
