## Install tess3r into the rule's conda environment, if it is not there yet.
##
## tess3r is only on GitHub. Its dependencies (Rcpp, RcppEigen) and a C++
## toolchain come from envs/tess.yml, so this compiles one small package and
## nothing else. upgrade = "never" stops remotes from rebuilding conda-supplied
## packages from source, which can leave the environment broken.

if (exists("snakemake") && length(snakemake@log) != 0) {
  log <- file(snakemake@log[1][[1]], open = "wt")
  sink(log, append = TRUE)
  sink(log, append = TRUE, type = "message")
}

lib <- .libPaths()[1]
# An interrupted install leaves 00LOCK-<pkg> behind, and every later attempt
# then fails with "failed to lock directory".
stale <- list.files(lib, pattern = "^00LOCK", full.names = TRUE)
if (length(stale)) {
  message("removing stale install lock(s): ", paste(basename(stale), collapse = ", "))
  unlink(stale, recursive = TRUE, force = TRUE)
}

ref <- snakemake@params[["ref"]]
if (requireNamespace("tess3r", quietly = TRUE)) {
  message("tess3r ", as.character(packageVersion("tess3r")), ": already installed")
} else {
  message("tess3r: installing bcm-uga/TESS3_encho_sen@", ref)
  remotes::install_github(paste0("bcm-uga/TESS3_encho_sen@", ref),
                          upgrade = "never", dependencies = FALSE, quiet = FALSE)
  if (!requireNamespace("tess3r", quietly = TRUE)) {
    stop("tess3r could not be installed; the messages above say why. This step ",
         "needs network access.", call. = FALSE)
  }
}

file.create(snakemake@output[["done"]])
