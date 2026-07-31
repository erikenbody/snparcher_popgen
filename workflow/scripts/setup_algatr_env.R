## Install the R packages that are not available on conda.
##
## algatr, wingen and tess3r are only on GitHub/CRAN, so they are compiled into
## the rule's conda environment on first use. This takes a long time (dozens of
## source dependencies) and people interrupt it, which is why the stale-lock
## cleanup below exists: R leaves a 00LOCK-<pkg> directory behind when an install
## is killed, and every later attempt then fails with "failed to lock directory"
## for that package -- and for everything that depends on it.

suppressMessages({
  library(ggplot2)
  library(terra)
  library(raster)
  library(fields)
  library(rworldmap)
  library(automap)
  library(cowplot)
  library(here)
  library(LEA)
  library(dplyr)
  library(rgdal)
  library(peakRAM)
})

output <- snakemake@output[[1]]

lib <- .libPaths()[1]
stale <- list.files(lib, pattern = "^00LOCK", full.names = TRUE)
if (length(stale)) {
  message("removing ", length(stale), " stale install lock(s) left by an ",
          "interrupted run: ", paste(basename(stale), collapse = ", "))
  unlink(stale, recursive = TRUE, force = TRUE)
}

install_if_missing <- function(pkg, installer) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(pkg, ": already installed")
    return(invisible(TRUE))
  }
  message(pkg, ": installing")
  try(installer(), silent = FALSE)
  # Checked rather than assumed: install.packages() and install_github() only
  # warn on failure, so without this the rule would create its sentinel and
  # report success, and TESS would fail later with an unrelated-looking error.
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " could not be installed. This step compiles many packages from ",
         "source and needs network access; the messages above say which ",
         "dependency failed. If an earlier attempt was interrupted, look for ",
         "leftover 00LOCK-* directories in ", lib, call. = FALSE)
  }
  message(pkg, ": installed")
}

# upgrade = "never" is load-bearing, not tidiness. The default ("ask", which is
# non-interactive here) rebuilds already-satisfied dependencies from source, and
# the heavy geospatial ones (s2, sf, units, terra) are supplied as conda binaries
# by envs/algatr.yml. A source rebuild of sf that fails on s2 does not fail
# harmlessly: it removes the working conda-installed sf first, so the environment
# ends up worse than before and every later run fails with "no package called
# 'sf'".
install_if_missing("algatr", function() {
  devtools::install_github("TheWangLab/algatr", quiet = TRUE, upgrade = "never")
})
install_if_missing("wingen", function() {
  install.packages("wingen", repos = "http://cran.us.r-project.org")
})
install_if_missing("tess3r", function() {
  devtools::install_github("bcm-uga/TESS3_encho_sen", quiet = TRUE,
                           upgrade = "never")
})

file.create(output)
