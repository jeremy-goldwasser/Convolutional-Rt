## ============================================================================
## install_packages.R
##
## Install all R package dependencies for the girt Rt pipeline (code/analyses).
## Run once after cloning:
##
##   Rscript code/install_packages.R
##
## (EpiNow2 brings rstan; the EpiNow2 fits themselves are run on a cluster.)
##
## ============================================================================

cat("Installing R package dependencies...\n\n")

# Helper: install from CRAN if not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s from CRAN...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(sprintf("  %s already installed (v%s)\n", pkg, packageVersion(pkg)))
  }
}

# Helper: install from GitHub if not already installed
install_github_if_missing <- function(repo, pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s from GitHub (%s)...\n", pkg, repo))
    if (!requireNamespace("remotes", quietly = TRUE))
      install.packages("remotes", repos = "https://cloud.r-project.org")
    remotes::install_github(repo)
  } else {
    cat(sprintf("  %s already installed (v%s)\n", pkg, packageVersion(pkg)))
  }
}

# ---- CRAN packages ----

cran_pkgs <- c(
  "dplyr",        # data manipulation
  "tidyr",        # data manipulation
  "tibble",       # data manipulation
  "ggplot2",      # plotting
  "EpiEstim",     # Rt estimation (renewal equation) baseline
  "rtestim",      # Rt estimation (trend-filtered Poisson) baseline
  "EpiLPS",       # Rt estimation (Bayesian P-splines) baseline
  "splines",      # B-spline basis (base R, always available)
  "EpiNow2",      # Bayesian Rt estimation (Stan-based); pulls rstan
  "rstan",        # Stan interface for R
  "remotes"       # for GitHub installs
)

# splines is a base package, skip install check
cran_pkgs <- setdiff(cran_pkgs, "splines")

cat("--- CRAN packages ---\n")
for (pkg in cran_pkgs) {
  install_if_missing(pkg)
}

# ---- GitHub packages ----

cat("\n--- GitHub packages ---\n")
install_github_if_missing("covid-19-Re/estimateR", "estimateR")

# ---- Verify all packages load ----

cat("\n--- Verification ---\n")
all_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "EpiEstim", "rtestim",
              "EpiLPS", "splines", "EpiNow2", "estimateR", "rstan")

ok <- TRUE
for (pkg in all_pkgs) {
  loaded <- tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    TRUE
  }, error = function(e) FALSE)
  status <- if (loaded) "OK" else "FAILED"
  cat(sprintf("  %-12s  %s\n", pkg, status))
  if (!loaded) ok <- FALSE
}

if (ok) {
  cat("\nAll packages installed and loadable.\n")
} else {
  cat("\nWARNING: Some packages failed to load. Check errors above.\n")
}
