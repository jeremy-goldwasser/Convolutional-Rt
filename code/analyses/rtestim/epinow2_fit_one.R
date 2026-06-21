# =============================================================================
# rtestim_paper_epinow2_fit_one.R
#
# Fit EpiNow2 GP (tight MCMC) on ONE scenario of the rtestim paper retro
# benchmark (SEIR variant only).  Parallels the other methods (Ours /
# estimateR / EpiEstim / rtestim) fit in rtestim_paper_retro_benchmark.r; this
# script is intended to run on the cluster (each EpiNow2 fit takes 25-45 min
# with tight settings, way too slow to do inline in the benchmark script).
#
# Observation channel: y = sim$seir_reports (infections from SEIR simulator,
# convolved forward through E->Y gamma).
#
# Scenarios (4, sit = "fake" for all):
#   1 = piecewise constant
#   2 = piecewise exponential
#   3 = piecewise linear
#   4 = periodic
#
# Priors (match rtestim_paper_retro_benchmark.r conventions):
#   GI            = paper "fake" SI: Gamma(mean = 8.4, sd = 3.8)
#   Delay (E->Y)  = from outputs/rtestim/data/seir_params.rds (mean_EY, sd_EY)
#   Rt prior      = LogNormal(meanlog = 0, sdlog = 1.0)
#                   -- sdlog bumped from the standard 0.5 because these
#                      scenarios have Rt reaching ~3 and dipping below 0.5;
#                      0.5 is too tight.
#   Obs family    = Poisson, no week_effect (no DoW in the simulator).
#
# MCMC: EpiNow2 defaults (4 chains x 500 warmup x 2000 samples, adapt 0.95,
#       treedepth 12).  cores = 4 -> chains run in parallel when the sbatch
#       allocates --cpus-per-task=4.
#
# Usage:
#   Rscript rtestim_paper_epinow2_fit_one.R <scenario_idx>
# e.g.
#   Rscript rtestim_paper_epinow2_fit_one.R 1
#
# Output:
#   outputs/rtestim/seir/data/epinow2/
#     fit_epinow2_<tag>.rds           (full list incl. $fit if save_full=TRUE)
#     fit_epinow2_<tag>__summary.rds  (summary: meta + timings + diag + rt_df)
# where <tag> = <sit>__<scenario_slug> (matches the existing benchmark's tag
# convention, so the retro benchmark loader can find them by scenario).
# =============================================================================

suppressPackageStartupMessages({
  library(EpiNow2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})
options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L)
  stop("Usage: Rscript rtestim_paper_epinow2_fit_one.R <scenario_idx>")
variant      <- "seir"
scenario_idx <- as.integer(args[1])
stopifnot(is.finite(scenario_idx), scenario_idx %in% 1:4)
obs_channel  <- "seir_reports"

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
script_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  this <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(this) == 0L) {
    sf <- try(sys.frame(1)$ofile, silent = TRUE)
    if (!inherits(sf, "try-error") && !is.null(sf)) this <- sf
  }
  if (length(this) == 0L || identical(this, "")) return(getwd())
  normalizePath(dirname(this), mustWork = FALSE)
})()
source(file.path(script_dir, "..", "..", "_paths.R"))
rtestim_data_dir    <- file.path(data_dir,    "sim", "rtestim")
rtestim_results_dir <- file.path(results_dir, "sim", "rtestim")
rtestim_figures_dir <- file.path(figures_dir, "sim", "rtestim")

# Reuse the flu EpiNow2 lib for extract_rt_df / extract_diagnostics /
# save_fit_result (writes summary-only by default, per our earlier cleanup).
source(file.path(code_dir, "analyses", "flu", "epinow2_lib.R"))

sim_data_dir <- rtestim_data_dir
# OUTPUT_TAG env var redirects writes to outputs/rtestim/<variant>/data/epinow2_<TAG>/.
.output_tag  <- Sys.getenv("OUTPUT_TAG", "")
.epinow2_dir <- if (nzchar(.output_tag)) sprintf("epinow2_%s", .output_tag) else "epinow2"
out_dir      <- file.path(rtestim_results_dir, variant,
                          "data", .epinow2_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load sim + delay params
# -----------------------------------------------------------------------------
sim     <- readRDS(file.path(sim_data_dir, "sim_combined.rds"))
sei_par <- readRDS(file.path(sim_data_dir, "seir_params.rds"))

combos <- unique(sim[, c("scenario", "si_type")])
combos <- combos[order(combos$si_type, combos$scenario), ]
if (scenario_idx > nrow(combos))
  stop(sprintf("scenario_idx=%d exceeds %d unique (scenario, si_type) combos",
               scenario_idx, nrow(combos)))
sc  <- combos$scenario[scenario_idx]
sit <- combos$si_type[scenario_idx]
tag <- sprintf("%s__%s", sit, gsub(":? ", "_", sc))

cat(sprintf("=== variant = %s  |  scenario_idx = %d\n", variant, scenario_idx))
cat(sprintf("    scenario = %s\n", sc))
cat(sprintf("    si_type  = %s\n", sit))
cat(sprintf("    tag      = %s\n", tag))

sub <- sim[sim$scenario == sc & sim$si_type == sit, , drop = FALSE]
sub <- sub[order(sub$time), ]
y   <- as.integer(pmax(0L, round(sub[[obs_channel]])))
date0 <- as.Date("2020-01-01")
dates <- date0 + seq(0, length(y) - 1L)
cat(sprintf("    obs_channel = %s   n_days = %d   max(y) = %d\n",
            obs_channel, length(y), max(y)))

p <- sei_par[sei_par$si_type == sit, ]
if (nrow(p) != 1L)
  stop(sprintf("Expected 1 row in seir_params for si_type='%s', got %d",
               sit, nrow(p)))

# TRUE SEIR-implied GI from seir_params.rds (matches MechRt + baselines on girt)
gi_mean <- p$theor_mean
gi_sd   <- p$theor_sd
delay_mean <- as.numeric(p$mean_EY)
delay_sd   <- as.numeric(p$sd_EY)

cat(sprintf("    GI: Gamma(mean=%.3f, sd=%.3f)\n",   gi_mean, gi_sd))
cat(sprintf("    Delay (E->Y): Gamma(mean=%.3f, sd=%.3f)\n",
            delay_mean, delay_sd))
cat(sprintf("    Env overrides: WARMUP=%s SAMPLES=%s CHAINS=%s ADAPT=%s TREE=%s\n",
            Sys.getenv("MCMC_WARMUP",    "(default)"),
            Sys.getenv("MCMC_SAMPLES",   "(default)"),
            Sys.getenv("MCMC_CHAINS",    "(default)"),
            Sys.getenv("MCMC_ADAPT_DELTA","(default)"),
            Sys.getenv("MCMC_TREEDEPTH", "(default)")))
cat(sprintf("    GP overrides: LS_MEAN=%s LS_SD=%s LS_MAX=%s BASIS=%s  rt_sdlog=%s\n",
            Sys.getenv("GP_LS_MEAN",    "(none)"),
            Sys.getenv("GP_LS_SD",      "(none)"),
            Sys.getenv("GP_LS_MAX",     "(default 60)"),
            Sys.getenv("GP_BASIS_PROP", "(default 0.2)"),
            Sys.getenv("RT_PRIOR_SDLOG","(default 1.0)")))
cat(sprintf("    GP kernel: KERNEL=%s MATERN_ORDER=%s\n",
            Sys.getenv("GP_KERNEL",       "(default matern)"),
            Sys.getenv("GP_MATERN_ORDER", "(default 3/2)")))
cat(sprintf("    GP prior: LS_DIST=%s ALPHA_SD=%s\n",
            Sys.getenv("GP_LS_DIST",  "(default normal)"),
            Sys.getenv("GP_ALPHA_SD", "(default 0.01)")))

# -----------------------------------------------------------------------------
# Build cfg object for epinow2_lib helpers
# -----------------------------------------------------------------------------
# --- Env-var overrides (parallel flu_based/epinow2_retrospective_fit_one.R) ---
.gp_basis_prop  <- suppressWarnings(as.numeric(Sys.getenv("GP_BASIS_PROP",  "")))
.gp_ls_mean     <- suppressWarnings(as.numeric(Sys.getenv("GP_LS_MEAN",     "")))
.gp_ls_sd       <- suppressWarnings(as.numeric(Sys.getenv("GP_LS_SD",       "")))
.gp_ls_max      <- suppressWarnings(as.numeric(Sys.getenv("GP_LS_MAX",      "")))
.mcmc_chains    <- suppressWarnings(as.integer(Sys.getenv("MCMC_CHAINS",    "")))
.mcmc_warmup    <- suppressWarnings(as.integer(Sys.getenv("MCMC_WARMUP",    "")))
.mcmc_samples   <- suppressWarnings(as.integer(Sys.getenv("MCMC_SAMPLES",   "")))
.mcmc_adapt     <- suppressWarnings(as.numeric(Sys.getenv("MCMC_ADAPT_DELTA","")))
.mcmc_treedepth <- suppressWarnings(as.integer(Sys.getenv("MCMC_TREEDEPTH", "")))
.rt_sdlog       <- suppressWarnings(as.numeric(Sys.getenv("RT_PRIOR_SDLOG", "")))
.gp_kernel       <- Sys.getenv("GP_KERNEL", "")
.gp_matern_order <- suppressWarnings(as.numeric(Sys.getenv("GP_MATERN_ORDER", "")))
.gp_ls_dist      <- Sys.getenv("GP_LS_DIST", "")
.gp_alpha_sd     <- suppressWarnings(as.numeric(Sys.getenv("GP_ALPHA_SD", "")))

# GRW (Gaussian random walk) overrides
.model_type      <- tolower(Sys.getenv("MODEL_TYPE", "gp"))
stopifnot(.model_type %in% c("gp", "grw"))
.rw_step         <- suppressWarnings(as.integer(Sys.getenv("RW_STEP", "")))
.bp_sd_prior_sd  <- suppressWarnings(as.numeric(Sys.getenv("BP_SD_PRIOR_SD", "")))
cat(sprintf("    MODEL_TYPE=%s RW_STEP=%s BP_SD_PRIOR_SD=%s\n",
            .model_type,
            Sys.getenv("RW_STEP",        "(default 1)"),
            Sys.getenv("BP_SD_PRIOR_SD", "(default 0.1, EpiNow2 hard-coded)")))

cfg <- list(
  gi_mean    = gi_mean,
  gi_sd      = gi_sd,
  delay_mean = delay_mean,
  delay_sd   = delay_sd,

  # Rt prior (slightly wider than flu/sim_v2 default; see header note).
  # Override via RT_PRIOR_SDLOG.
  rt_prior_meanlog = 0,
  rt_prior_sdlog   = if (is.finite(.rt_sdlog)) .rt_sdlog else 1.0,

  # Poisson-only sim, no DoW
  obs_family  = "poisson",
  week_effect = FALSE,

  # MCMC defaults (EpiNow2 defaults: 4/500/2000/0.95/12).  Override every
  # parameter via env (see flu_based for precedent).  Note EpiNow2's `samples`
  # is TOTAL post-warmup samples across chains, so e.g. samples=4000 chains=4
  # gives 1000 post-warmup draws/chain.
  mcmc_gp = list(
    chains        = if (!is.na(.mcmc_chains))    .mcmc_chains    else 4L,
    cores         = if (!is.na(.mcmc_chains))    .mcmc_chains    else 4L,
    warmup        = if (!is.na(.mcmc_warmup))    .mcmc_warmup    else 500L,
    samples       = if (!is.na(.mcmc_samples))   .mcmc_samples   else 2000L,
    adapt_delta   = if (is.finite(.mcmc_adapt))  .mcmc_adapt     else 0.95,
    max_treedepth = if (!is.na(.mcmc_treedepth)) .mcmc_treedepth else 12L,
    seed          = 42L + scenario_idx
  ),

  gp_basis_prop = if (is.finite(.gp_basis_prop)) .gp_basis_prop else 0.2,
  gp_ls_mean    = if (is.finite(.gp_ls_mean))    .gp_ls_mean    else NA_real_,
  gp_ls_sd      = if (is.finite(.gp_ls_sd))      .gp_ls_sd      else NA_real_,
  gp_ls_max     = if (is.finite(.gp_ls_max))     .gp_ls_max     else NA_real_,
  gp_kernel       = if (nzchar(.gp_kernel))            .gp_kernel       else NULL,
  gp_matern_order = if (is.finite(.gp_matern_order))   .gp_matern_order else NULL,
  gp_ls_dist      = if (nzchar(.gp_ls_dist))           .gp_ls_dist      else NULL,
  gp_alpha_sd     = if (is.finite(.gp_alpha_sd))       .gp_alpha_sd     else NULL,

  # GRW config (used when MODEL_TYPE=grw).  rw=1 means each day gets its own
  # breakpoint, so the random walk has 1 step per day -- enough resolution to
  # absorb large jumps (e.g. S1's 2.0->0.8 = -0.92 log-step) over a few days.
  # The step SD prior bp_sd ~ HalfNormal(0, 0.1) is hard-coded in EpiNow2's
  # Stan code; the hierarchical prior lets the data push it up if needed.
  grw_rw_step = if (!is.na(.rw_step)) .rw_step else 1L,
  mcmc_grw = list(
    chains        = if (!is.na(.mcmc_chains))    .mcmc_chains    else 4L,
    cores         = if (!is.na(.mcmc_chains))    .mcmc_chains    else 4L,
    warmup        = if (!is.na(.mcmc_warmup))    .mcmc_warmup    else 1000L,
    samples       = if (!is.na(.mcmc_samples))   .mcmc_samples   else 2000L,
    adapt_delta   = if (is.finite(.mcmc_adapt))  .mcmc_adapt     else 0.95,
    max_treedepth = if (!is.na(.mcmc_treedepth)) .mcmc_treedepth else 12L,
    seed          = 42L + scenario_idx
  ),

  model_types   = c("gp", "grw"),

  # Metadata fields used by epinow2_lib's fit wrapper (otherwise ignored)
  season        = sprintf("rtestim_%s_sc%d", variant, scenario_idx),
  mode          = "retrospective"
)

# -----------------------------------------------------------------------------
# Fit
# -----------------------------------------------------------------------------
reported_cases <- data.frame(date = dates, confirm = y,
                             stringsAsFactors = FALSE)

gen_time <- build_generation_time(cfg)
delays   <- build_delays(cfg)
st_opts  <- build_stan_opts(cfg, .model_type)
rt_o     <- build_rt_opts(cfg, .model_type)
gp_o     <- build_gp_opts_or_null(cfg, .model_type)
obs_o    <- build_obs_opts(cfg)

cat(sprintf("\nStarting EpiNow2 %s fit ...\n", toupper(.model_type)))
t_start <- Sys.time()
fit <- EpiNow2::estimate_infections(
  data            = reported_cases,
  generation_time = gen_time,
  delays          = delays,
  rt              = rt_o,
  stan            = st_opts,
  obs             = obs_o,
  gp              = gp_o,
  verbose         = TRUE
)
t_end <- Sys.time()
wall_sec <- as.numeric(difftime(t_end, t_start, units = "secs"))
cat(sprintf("Done in %.1f min\n", wall_sec / 60))

# -----------------------------------------------------------------------------
# Post-process + save (summary-only by default)
# -----------------------------------------------------------------------------
rt_df <- tryCatch(extract_rt_df(fit), error = function(e) {
  cat(sprintf("  Rt extraction failed: %s\n", conditionMessage(e)))
  NULL
})
diag <- tryCatch(extract_diagnostics(fit), error = function(e) {
  list(note = conditionMessage(e))
})

meta <- list(
  source          = "rtestim_paper",
  variant         = variant,
  obs_channel     = obs_channel,
  scenario        = sc,
  si_type         = sit,
  tag             = tag,
  scenario_idx    = scenario_idx,
  n_days          = length(y),
  first_date      = min(dates),
  last_date       = max(dates),
  gi_mean         = gi_mean,
  gi_sd           = gi_sd,
  delay_mean      = delay_mean,
  delay_sd        = delay_sd,
  rt_prior_meanlog = cfg$rt_prior_meanlog,
  rt_prior_sdlog   = cfg$rt_prior_sdlog,
  obs_family      = cfg$obs_family,
  week_effect     = cfg$week_effect,
  model_type      = .model_type,
  grw_rw_step     = if (.model_type == "grw") cfg$grw_rw_step else NA_integer_,
  mcmc            = if (.model_type == "gp")  cfg$mcmc_gp else cfg$mcmc_grw,
  gp_basis_prop   = cfg$gp_basis_prop,
  gp_ls_mean      = cfg$gp_ls_mean,
  gp_ls_sd        = cfg$gp_ls_sd,
  gp_ls_max       = cfg$gp_ls_max,
  gp_kernel       = cfg$gp_kernel,
  gp_matern_order = cfg$gp_matern_order,
  gp_ls_dist      = cfg$gp_ls_dist,
  gp_alpha_sd     = cfg$gp_alpha_sd,
  epinow2_version = as.character(utils::packageVersion("EpiNow2")),
  r_version       = R.version.string,
  fit_started_at  = format(t_start, "%Y-%m-%d %H:%M:%S"),
  fit_finished_at = format(t_end,   "%Y-%m-%d %H:%M:%S")
)
timings <- list(wall_total_sec = wall_sec, wall_total_min = wall_sec / 60)

result <- list(meta = meta, timings = timings, diagnostics = diag,
               rt_df = rt_df, fit = fit)

out_path <- file.path(out_dir, sprintf("fit_epinow2_%s.rds", tag))
save_fit_result(result, out_path)   # writes summary-only by default

cat(sprintf("\nSaved: %s\n", sub("\\.rds$", "__summary.rds", out_path)))
if (!is.null(diag)) {
  cat(sprintf("Diagnostics: max_rhat=%s  min_ess=%s  div=%s  tree=%s\n",
              format(diag$max_rhat,     digits = 3),
              format(diag$min_ess_bulk, digits = 3),
              format(diag$n_divergent),
              format(diag$n_max_treedepth)))
}
