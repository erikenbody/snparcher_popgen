suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
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

# example call: `Rscript admixture_cloud.R "59-Ursus" "outputs/admixture/"`

#!/usr/bin/env Rscript # leave line commented

species = snakemake@params[[1]]
output_path = snakemake@params[[2]]


#' Marker for Evanno's delta K, omitted when it is undefined.
evanno_layers <- function(cv_scores, y = NULL) {
  k <- unique(cv_scores$evanno_k)
  if (length(k) == 0 || all(is.na(k))) return(NULL)
  layers <- list(ggplot2::geom_vline(xintercept = k, color = "purple", linetype = "dashed"))
  if (!is.null(y)) {
    layers <- c(layers, list(ggplot2::annotate(
      "text", x = k - 0.2, y = y, label = "\u0394K", size = 4, color = "purple", angle = 90)))
  }
  layers
}


# Get CV scores -----------------------------------------------------------

admix_cv <- function(output_path) {
  files <- list.files(output_path, pattern = ".log", full.names = TRUE)
  shortfiles <- list.files(output_path, pattern = ".log", full.names = FALSE)
  
  cv_scores <-
    1:length(files) %>% 
    lapply(function(x) {
      cv <- readLines(con = files[[x]])
      cv <- cv[grepl("^CV error", cv)]
      cv <- as.data.frame(cv) %>% dplyr::mutate(filename = shortfiles[[x]])
      return(cv)
    }) %>% 
    dplyr::bind_rows() %>%
    tidyr::separate(cv, sep = ": ", into = c("temp", "cv")) %>%
    # K and seed come from the log file name. Matched by pattern rather than by
    # splitting on "_", so that project prefixes containing underscores work.
    dplyr::mutate(kval = as.numeric(stringr::str_match(filename, "_K([0-9]+)_")[, 2]),
                  seed = as.numeric(stringr::str_match(filename, "_seed([0-9]+)\\.log$")[, 2])) %>%
    dplyr::select(-c(temp, filename)) %>%
    dplyr::mutate(cv = as.numeric(cv),
                  kval = as.numeric(kval), 
                  seed = as.numeric(seed)) %>% 
    dplyr::group_by(kval) %>% 
    dplyr::mutate(min_cv_seed = min(cv, na.rm = T)) %>% #EDE added na.rm = T here. This allows seeds that fail to includes, seems fine though
    dplyr::ungroup()

  cv_scores_seeds <- cv_scores %>% 
    dplyr::arrange(kval) %>% 
    dplyr::select(kval, min_cv_seed) %>% 
    distinct()
  
  # Get min_k value across all Ks
  min_k <- cv_scores_seeds %>% 
    dplyr::filter(min_cv_seed == min(min_cv_seed)) %>%
    dplyr::pull(kval)
  
  # Get best_k value
  kvals <- cv_scores_seeds$kval
  ce.K <- cv_scores_seeds$min_cv_seed
  # index by length, not by max(kval): the two only coincide when K starts at 1
  diff <- ce.K[-1] - ce.K[-length(ce.K)]
  slope <- exp(-diff) - 1
  best_k <- if (length(slope) == 0) kvals[1] else min(which(slope <= quantile(slope)[4]))

  # Evanno's delta K value; adapted from code provided by Katherine Drotos
  delta_k <- cv_scores %>%
    dplyr::arrange(seed, kval) %>%
    dplyr::mutate(secondrate = abs((lead(cv) - cv) - (cv - lag(cv)))) %>%
    dplyr::group_by(kval) %>%
    dplyr::mutate(delta_k = mean(secondrate) / sd(cv)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(kval, delta_k)
  
  # Evanno's delta K is a second difference, so it needs at least three K
  # values. With fewer, or when every delta is NaN, report it as missing rather
  # than failing: min_k and best_k are still perfectly usable.
  evanno_k <- delta_k$kval[which.max(delta_k$delta_k)]
  if (length(evanno_k) == 0) {
    evanno_k <- NA_real_
    message("Evanno delta K is undefined for this K range (needs >= 3 values of K)")
  }
  
  cv_scores <- cv_scores %>% 
    dplyr::mutate(min_k = min_k,
                  best_k = best_k,
                  evanno_k = evanno_k)
  
  return(cv_scores)
}

# Run above function
cv_scores <- admix_cv(output_path)


# Export results ----------------------------------------------------------

export_admix <- function(output_path, species, cv_scores) {
  readr::write_csv(cv_scores,
                   file = paste0(output_path, species, "_best_K.txt"),
                   col_names = TRUE)
  
  # Export plot of CV error
  cv_scores %>%
    dplyr::group_by(kval) %>% 
    dplyr::mutate(mean = mean(cv)) %>% 
    dplyr::mutate(sd = sd(cv)) %>% 
    dplyr::ungroup() %>% 
    ggplot2::ggplot(aes(x = kval, y = mean)) +
    ggplot2::geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2, color = "grey") +
    ggplot2::geom_point() +
    ggplot2::geom_line(aes(group = 1)) +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(xintercept = unique(cv_scores$min_k),
                        color = "red",
                        linetype = "dashed") +
    ggplot2::annotate("text", x = unique(cv_scores$min_k-0.2), y = mean(cv_scores$cv), label = "Minimum K", size = 4, color = "red", angle = 90) +
    ggplot2::geom_vline(xintercept = unique(cv_scores$best_k),
                        color = "blue",
                        linetype = "dashed") +
    ggplot2::annotate("text", x = unique(cv_scores$best_k-0.2), y = mean(cv_scores$cv), label = "Best K", size = 4, color = "blue", angle = 90) +
    evanno_layers(cv_scores, y = max(cv_scores$cv)) +
    ggplot2::ggtitle(paste0(species, " admixture results")) +
    ggplot2::ylab("Cross-validation score") +
    ggplot2::scale_x_continuous(breaks = unique(cv_scores$kval))

  ggplot2::ggsave(paste0(output_path, species, "_admixture_cv.png"),
                  width = 20, height = 15, units = "cm")
  
  # Export Evanno's delta K plot; adapted from code provided by Katherine Drotos
  cv_scores %>%
    dplyr::arrange(seed, kval) %>%
    dplyr::mutate(secondrate = abs((lead(cv) - cv) - (cv - lag(cv)))) %>%
    dplyr::group_by(kval) %>%
    dplyr::mutate(delta_k = mean(secondrate) / sd(cv)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(kval, delta_k) %>% 
    ggplot(aes(x = kval, y = delta_k)) +
    evanno_layers(cv_scores) +
    geom_point() + 
    geom_line(aes(group = 1)) +
    scale_x_continuous(breaks = unique(cv_scores$kval)) +
    ylab("Δk") +
    ggtitle(paste0(species, " delta K"))
  
  ggplot2::ggsave(paste0(output_path, species, "_admixture_deltak.png"),
                  width = 20, height = 15, units = "cm")
}

export_admix(output_path, species, cv_scores)