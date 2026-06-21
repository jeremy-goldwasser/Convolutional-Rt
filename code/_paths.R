# =============================================================================
# Code/_paths.R
#
# Defines repo-relative paths that scripts use instead of hard-coded absolute
# paths. Source this file at the top of any analysis script:
#
#   source("<path-to-Code>/_paths.R")
#
# After sourcing you have:
#   repo_root, code_dir, mechrt_dir, data_dir, results_dir, figures_dir
#
# Paths under each of the data / results / figures roots:
#   <X>/real/covidestim/      (real COVID data: covidestim outputs)
#   <X>/real/flu/             (real flu data: HHS hospitalizations)
#   <X>/sim/flu_based/  (flu-based GP-Rt simulation: tuning_sim_v2)
#   <X>/sim/rtestim/    (rtestim-paper SEIR simulations)
# =============================================================================

.mechrt_find_repo_root <- function() {
  # 1. Try to locate this file's path
  args <- commandArgs(trailingOnly = FALSE)
  this_arg <- args[grep("^--file=", args)]
  this <- if (length(this_arg)) sub("^--file=", "", this_arg[1]) else NULL
  if (is.null(this) && exists("sys.frame", mode = "function")) {
    sf <- try(sys.frame(1)$ofile, silent = TRUE)
    if (!inherits(sf, "try-error") && !is.null(sf)) this <- sf
  }
  start <- if (is.null(this) || identical(this, "")) getwd()
           else dirname(normalizePath(this, mustWork = FALSE))
  # 2. Walk upward looking for the four canonical top-level folders
  d <- start
  for (.. in 1:8) {
    needed <- c("code", "data", "results", "figures")
    if (all(needed %in% list.files(d))) return(normalizePath(d))
    p <- dirname(d)
    if (p == d) break
    d <- p
  }
  # 3. Last resort: MECHRT_REPO_ROOT env var (must point at a valid root).
  #    No silent getwd() fallback — that's how we ended up writing
  #    figures/results/sim/... when getwd() happened to be figures/.
  env <- Sys.getenv("MECHRT_REPO_ROOT", unset = "")
  if (nzchar(env)) {
    env_norm <- normalizePath(env, mustWork = FALSE)
    if (all(c("code", "data", "results", "figures") %in% list.files(env_norm))) {
      return(env_norm)
    }
    stop("MECHRT_REPO_ROOT='", env, "' does not contain code/data/results/figures.")
  }
  stop("Could not locate repo root from this=", if (is.null(this)) "<NULL>" else this,
       " (start=", start, "); set MECHRT_REPO_ROOT or run from inside the repo.")
}

repo_root   <- .mechrt_find_repo_root()
stopifnot(all(c("code", "data", "results", "figures") %in% list.files(repo_root)))
code_dir    <- file.path(repo_root, "code")
mechrt_dir  <- file.path(code_dir,  "mechrt")
data_dir    <- file.path(repo_root, "data")
results_dir <- file.path(repo_root, "results")
figures_dir <- file.path(repo_root, "figures")

# Convenience: the canonical MechRt loader file path.
mechrt_loader <- file.path(mechrt_dir, "mechrt.R")
