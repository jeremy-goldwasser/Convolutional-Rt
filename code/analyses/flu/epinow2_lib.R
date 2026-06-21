# =============================================================================
# epinow2_lib.R
#
# Shared library for flu real-time EpiNow2 fitting experiments.
# Sourced by epinow2_poc.R, epinow2_fit_one.R.
#
# Exposes:
#   load_flu_daily(cfg)
#   build_generation_time(cfg)
#   build_delays(cfg)
#   build_stan_opts(cfg, model_type)
#   build_rt_opts(cfg, model_type)
#   build_gp_opts_or_null(cfg, model_type)
#   build_obs_opts(cfg)
#   fit_one_epinow2(cfg, hosps_all, week_end_date, model_type, verbose=TRUE)
#   extract_rt_df(fit)
#   extract_diagnostics(fit)
#   save_fit_result(result, path)
#
# Caller must have already library(EpiNow2)'d and source()'d epinow2_config.R.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# Data loader
# -----------------------------------------------------------------------------
load_flu_daily <- function(cfg, root_dir = NULL) {
  # Prefer cfg$daily_clean_path (absolute, set by realtime_config.R).  Fall
  # back to root_dir + cfg$daily_clean_relpath for legacy callers.
  if (!is.null(cfg$daily_clean_path)) {
    path <- cfg$daily_clean_path
  } else {
    path <- file.path(root_dir, cfg$daily_clean_relpath)
  }
  if (!file.exists(path))
    stop(sprintf("Cleaned daily flu data not found: %s", path))
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$date <- as.Date(df$date)
  df <- df[order(df$date), ]
  if (max(df$date) < cfg$max_data_date)
    stop(sprintf("Data file only covers up to %s, need >= %s",
                 max(df$date), cfg$max_data_date))
  if (!(cfg$data_start %in% df$date))
    stop(sprintf("data_start %s not in data file", cfg$data_start))
  df
}

# -----------------------------------------------------------------------------
# Distributions (generation interval and reporting delay)
# -----------------------------------------------------------------------------
build_generation_time <- function(cfg) {
  gi_shape <- (cfg$gi_mean / cfg$gi_sd)^2
  gi_rate  <- cfg$gi_mean / cfg$gi_sd^2
  gi_max   <- ceiling(qgamma(0.999, shape = gi_shape, rate = gi_rate))
  EpiNow2::generation_time_opts(
    EpiNow2::Gamma(mean = cfg$gi_mean, sd = cfg$gi_sd, max = gi_max)
  )
}

build_delays <- function(cfg) {
  d_shape <- (cfg$delay_mean / cfg$delay_sd)^2
  d_rate  <- cfg$delay_mean / cfg$delay_sd^2
  d_max   <- ceiling(qgamma(0.9995, shape = d_shape, rate = d_rate))
  EpiNow2::delay_opts(
    EpiNow2::Gamma(mean = cfg$delay_mean, sd = cfg$delay_sd, max = d_max)
  )
}

# -----------------------------------------------------------------------------
# MCMC / stan options (differ between GP and GRW)
# -----------------------------------------------------------------------------
build_stan_opts <- function(cfg, model_type) {
  stopifnot(model_type %in% cfg$model_types)
  m <- if (model_type == "gp") cfg$mcmc_gp else cfg$mcmc_grw
  # cores = chains -> run chains in parallel whenever the sbatch allocates
  # --cpus-per-task >= chains.  Falls back to serial if fewer cores available.
  EpiNow2::stan_opts(
    method  = "sampling",
    chains  = m$chains,
    cores   = if (!is.null(m$cores)) m$cores else m$chains,
    warmup  = m$warmup,
    samples = m$samples,
    seed    = m$seed,
    control = list(
      adapt_delta   = m$adapt_delta,
      max_treedepth = m$max_treedepth
    )
  )
}

build_rt_opts <- function(cfg, model_type) {
  stopifnot(model_type %in% cfg$model_types)
  prior <- EpiNow2::LogNormal(
    meanlog = cfg$rt_prior_meanlog,
    sdlog   = cfg$rt_prior_sdlog
  )
  if (model_type == "gp") {
    EpiNow2::rt_opts(prior = prior)
  } else {
    EpiNow2::rt_opts(prior = prior, rw = as.integer(cfg$grw_rw_step))
  }
}

build_gp_opts_or_null <- function(cfg, model_type) {
  stopifnot(model_type %in% cfg$model_types)
  if (model_type == "gp") {
    args <- list(basis_prop = cfg$gp_basis_prop)
    # Optional length-scale prior override.  Newer EpiNow2 deprecated
    # ls_mean/ls_sd; the prior is now passed as a distribution object via the
    # `ls` arg.  Default in EpiNow2 is Normal(mean = 21, sd = 7, max = 60).
    # Set cfg$gp_ls_mean + cfg$gp_ls_sd to enable a wigglier prior.
    if (!is.null(cfg$gp_ls_mean) && is.finite(cfg$gp_ls_mean) &&
        !is.null(cfg$gp_ls_sd)   && is.finite(cfg$gp_ls_sd)) {
      ls_max <- if (!is.null(cfg$gp_ls_max) && is.finite(cfg$gp_ls_max))
                  cfg$gp_ls_max else 60
      ls_dist <- if (!is.null(cfg$gp_ls_dist) && nzchar(cfg$gp_ls_dist))
                   tolower(cfg$gp_ls_dist) else "normal"
      args$ls <- if (ls_dist == "lognormal")
        EpiNow2::LogNormal(mean = cfg$gp_ls_mean,
                           sd   = cfg$gp_ls_sd, max = ls_max)
      else
        EpiNow2::Normal(mean = cfg$gp_ls_mean,
                        sd   = cfg$gp_ls_sd, max = ls_max)
      cat(sprintf("  GP ls override: %s(mean=%.1f, sd=%.1f, max=%.1f)\n",
                  if (ls_dist == "lognormal") "LogNormal" else "Normal",
                  cfg$gp_ls_mean, cfg$gp_ls_sd, ls_max))
    }
    if (!is.null(cfg$gp_alpha_sd) && is.finite(cfg$gp_alpha_sd)) {
      args$alpha <- EpiNow2::Normal(mean = 0, sd = cfg$gp_alpha_sd)
      cat(sprintf("  GP alpha override: Normal(mean=0, sd=%.3f)\n",
                  cfg$gp_alpha_sd))
    }
    if (!is.null(cfg$gp_kernel) && nzchar(cfg$gp_kernel)) {
      args$kernel <- cfg$gp_kernel
      cat(sprintf("  GP kernel override: %s\n", cfg$gp_kernel))
    }
    if (!is.null(cfg$gp_matern_order) && is.finite(cfg$gp_matern_order)) {
      args$matern_order <- cfg$gp_matern_order
      cat(sprintf("  GP matern_order override: %g\n", cfg$gp_matern_order))
    }
    do.call(EpiNow2::gp_opts, args)
  } else {
    NULL  # disable GP for GRW
  }
}

build_obs_opts <- function(cfg) {
  EpiNow2::obs_opts(
    family      = cfg$obs_family,
    week_effect = cfg$week_effect
  )
}

# -----------------------------------------------------------------------------
# Extract Rt summary data frame (the primary output)
# -----------------------------------------------------------------------------
extract_rt_df <- function(fit) {
  params <- tryCatch(
    as.data.frame(summary(fit, type = "parameters")),
    error = function(e) {
      # Fallback: use $estimates$summarised directly
      df <- fit$estimates$summarised
      if (is.null(df)) stop("No parameters summary found in fit")
      as.data.frame(df)
    }
  )

  rt_rows <- params[
    params$variable == "R" &
      params$type %in% c("estimate", "estimate based on partial data",
                         "forecast"),
    , drop = FALSE]

  if (nrow(rt_rows) == 0L) stop("No Rt rows found in EpiNow2 output.")

  # EpiNow2 exposes mean/median/sd + lower_N / upper_N quantiles.  Be robust:
  # pull whichever columns exist.
  cn <- colnames(rt_rows)
  pull_col <- function(nm) if (nm %in% cn) rt_rows[[nm]] else rep(NA_real_, nrow(rt_rows))

  data.frame(
    date      = as.Date(rt_rows$date),
    type      = as.character(rt_rows$type),
    Rt_mean   = pull_col("mean"),
    Rt_med    = pull_col("median"),
    Rt_sd     = pull_col("sd"),
    Rt_lo_90  = pull_col("lower_90"),
    Rt_hi_90  = pull_col("upper_90"),
    Rt_lo_50  = pull_col("lower_50"),
    Rt_hi_50  = pull_col("upper_50"),
    Rt_lo_20  = pull_col("lower_20"),
    Rt_hi_20  = pull_col("upper_20"),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Extract MCMC diagnostics from the stanfit object inside an EpiNow2 fit.
# Robust to both rstan::stanfit and CmdStanFit backends.
# -----------------------------------------------------------------------------
extract_diagnostics <- function(fit) {
  out <- list(
    backend         = NA_character_,
    n_chains        = NA_integer_,
    warmup          = NA_integer_,
    samples         = NA_integer_,
    max_rhat        = NA_real_,
    min_ess_bulk    = NA_real_,
    min_ess_tail    = NA_real_,
    n_divergent     = NA_integer_,
    n_max_treedepth = NA_integer_,
    elapsed_warmup  = NA_real_,
    elapsed_sample  = NA_real_,
    elapsed_total   = NA_real_,
    rt_max_rhat     = NA_real_,
    rt_min_ess_bulk = NA_real_,
    note            = ""
  )

  sf <- tryCatch(fit$fit, error = function(e) NULL)
  if (is.null(sf)) {
    out$note <- "fit$fit is NULL; no diagnostics available"
    return(out)
  }

  # rstan::stanfit path
  if (inherits(sf, "stanfit")) {
    out$backend <- "rstan"
    out$n_chains <- tryCatch(length(sf@stan_args), error = function(e) NA_integer_)
    args1 <- tryCatch(sf@stan_args[[1]], error = function(e) list())
    out$warmup  <- tryCatch(as.integer(args1$warmup),  error = function(e) NA_integer_)
    out$samples <- tryCatch(as.integer(args1$iter - args1$warmup),
                             error = function(e) NA_integer_)

    summ <- tryCatch(rstan::summary(sf)$summary, error = function(e) NULL)
    if (!is.null(summ)) {
      # Skip fixed/constant params (sd == 0): these are transformed-data
      # vectors like `gt_rev_pmf[]` that EpiNow2 exposes as Stan parameters.
      # rstan's summary formula gives them artifact Rhat/n_eff (n_eff ~ 0.5
      # on every fit), which contaminates the min/max headline numbers.
      live <- !is.na(summ[, "sd"]) & summ[, "sd"] > 0
      summ_live <- summ[live, , drop = FALSE]
      out$max_rhat     <- suppressWarnings(max(summ_live[, "Rhat"],  na.rm = TRUE))
      out$min_ess_bulk <- suppressWarnings(min(summ_live[, "n_eff"], na.rm = TRUE))
      out$min_ess_tail <- NA_real_  # rstan does not separately report ess_tail
      # Rt parameters: match names like "R[.*]"
      rt_rows <- grep("^R\\[", rownames(summ_live), value = TRUE)
      if (length(rt_rows) > 0) {
        out$rt_max_rhat     <- suppressWarnings(max(summ_live[rt_rows, "Rhat"],  na.rm = TRUE))
        out$rt_min_ess_bulk <- suppressWarnings(min(summ_live[rt_rows, "n_eff"], na.rm = TRUE))
      }
    }

    out$n_divergent <- tryCatch({
      div <- rstan::get_num_divergent(sf)
      as.integer(sum(div))
    }, error = function(e) NA_integer_)
    out$n_max_treedepth <- tryCatch({
      mt <- rstan::get_num_max_treedepth(sf)
      as.integer(sum(mt))
    }, error = function(e) NA_integer_)

    elapsed <- tryCatch(rstan::get_elapsed_time(sf), error = function(e) NULL)
    if (!is.null(elapsed)) {
      out$elapsed_warmup <- suppressWarnings(sum(elapsed[, "warmup"]))
      out$elapsed_sample <- suppressWarnings(sum(elapsed[, "sample"]))
      out$elapsed_total  <- out$elapsed_warmup + out$elapsed_sample
    }
    return(out)
  }

  # CmdStanFit path
  if (inherits(sf, c("CmdStanFit", "CmdStanMCMC"))) {
    out$backend <- "cmdstanr"
    out$n_chains <- tryCatch(sf$num_chains(), error = function(e) NA_integer_)
    meta <- tryCatch(sf$metadata(), error = function(e) list())
    out$warmup  <- tryCatch(as.integer(meta$iter_warmup),   error = function(e) NA_integer_)
    out$samples <- tryCatch(as.integer(meta$iter_sampling), error = function(e) NA_integer_)

    summ <- tryCatch(sf$summary(), error = function(e) NULL)
    if (!is.null(summ)) {
      # Drop fixed/constant params (sd == 0); see rstan branch comment.
      summ_live <- summ[!is.na(summ$sd) & summ$sd > 0, , drop = FALSE]
      out$max_rhat     <- suppressWarnings(max(summ_live$rhat,     na.rm = TRUE))
      out$min_ess_bulk <- suppressWarnings(min(summ_live$ess_bulk, na.rm = TRUE))
      out$min_ess_tail <- suppressWarnings(min(summ_live$ess_tail, na.rm = TRUE))
      rt_rows <- summ_live[grep("^R\\[", summ_live$variable), , drop = FALSE]
      if (nrow(rt_rows) > 0) {
        out$rt_max_rhat     <- suppressWarnings(max(rt_rows$rhat,     na.rm = TRUE))
        out$rt_min_ess_bulk <- suppressWarnings(min(rt_rows$ess_bulk, na.rm = TRUE))
      }
    }

    diag <- tryCatch(sf$diagnostic_summary(quiet = TRUE),
                     error = function(e) NULL)
    if (!is.null(diag)) {
      out$n_divergent     <- tryCatch(as.integer(sum(diag$num_divergent)),
                                      error = function(e) NA_integer_)
      out$n_max_treedepth <- tryCatch(as.integer(sum(diag$num_max_treedepth)),
                                      error = function(e) NA_integer_)
    }
    tt <- tryCatch(sf$time(), error = function(e) NULL)
    if (!is.null(tt) && !is.null(tt$total))
      out$elapsed_total <- as.numeric(tt$total)
    return(out)
  }

  out$note <- sprintf("unrecognised fit$fit class: %s",
                      paste(class(sf), collapse = ","))
  out
}

# -----------------------------------------------------------------------------
# Fit one EpiNow2 job
#
# Returns a list with the full $fit, extracted $rt_df, $diagnostics,
# $timings, $meta.  save_fit_result(...) writes it to disk.
# -----------------------------------------------------------------------------
fit_one_epinow2 <- function(cfg, hosps_all, week_end_date, model_type,
                            verbose = TRUE) {
  stopifnot(model_type %in% cfg$model_types)
  week_end_date <- as.Date(week_end_date)

  # --- Slice data: cfg$data_start ... week_end_date -------------------------
  reported_cases <- hosps_all %>%
    dplyr::filter(.data$date >= cfg$data_start,
                  .data$date <= week_end_date) %>%
    dplyr::transmute(date   = .data$date,
                     confirm = as.integer(pmax(0, round(.data$hosps_1d))))

  if (nrow(reported_cases) < 60L)
    stop(sprintf("Too few days (%d) for week_end=%s; need >= 60",
                 nrow(reported_cases), week_end_date))

  cat(sprintf("\n=== EpiNow2 fit: %s | %s | %s ===\n",
              cfg$season, week_end_date, model_type))
  cat(sprintf("  Data: %d days (%s -> %s)\n",
              nrow(reported_cases), min(reported_cases$date),
              max(reported_cases$date)))
  cat(sprintf("  Peak: %d\n", max(reported_cases$confirm, na.rm = TRUE)))

  gen_time <- build_generation_time(cfg)
  delays   <- build_delays(cfg)
  st_opts  <- build_stan_opts(cfg, model_type)
  rt_o     <- build_rt_opts(cfg, model_type)
  gp_o     <- build_gp_opts_or_null(cfg, model_type)
  obs_o    <- build_obs_opts(cfg)

  # --- Fit ------------------------------------------------------------------
  t_start <- Sys.time()
  fit <- EpiNow2::estimate_infections(
    data            = reported_cases,
    generation_time = gen_time,
    delays          = delays,
    rt              = rt_o,
    stan            = st_opts,
    obs             = obs_o,
    gp              = gp_o,
    verbose         = verbose
  )
  t_end <- Sys.time()

  wall_sec <- as.numeric(difftime(t_end, t_start, units = "secs"))
  cat(sprintf("  Wall time: %.1f sec (%.1f min)\n", wall_sec, wall_sec / 60))

  # --- Post-processing ------------------------------------------------------
  rt_df <- tryCatch(extract_rt_df(fit), error = function(e) {
    cat(sprintf("  Rt extraction failed: %s\n", conditionMessage(e)))
    NULL
  })

  diag <- tryCatch(extract_diagnostics(fit), error = function(e) {
    cat(sprintf("  Diagnostics extraction failed: %s\n", conditionMessage(e)))
    list(note = conditionMessage(e))
  })

  # Meta: compact record of everything needed to reproduce the fit
  meta <- list(
    season          = cfg$season,
    week_end_date   = week_end_date,
    model_type      = model_type,
    n_days          = nrow(reported_cases),
    first_date      = min(reported_cases$date),
    last_date       = max(reported_cases$date),
    gi_mean         = cfg$gi_mean,
    gi_sd           = cfg$gi_sd,
    delay_mean      = cfg$delay_mean,
    delay_sd        = cfg$delay_sd,
    rt_prior_meanlog = cfg$rt_prior_meanlog,
    rt_prior_sdlog   = cfg$rt_prior_sdlog,
    obs_family      = cfg$obs_family,
    week_effect     = cfg$week_effect,
    grw_rw_step     = if (model_type == "grw") cfg$grw_rw_step else NA_integer_,
    gp_basis_prop   = if (model_type == "gp")  cfg$gp_basis_prop else NA_real_,
    mcmc            = if (model_type == "gp")  cfg$mcmc_gp else cfg$mcmc_grw,
    epinow2_version = as.character(utils::packageVersion("EpiNow2")),
    r_version       = R.version.string,
    fit_started_at  = format(t_start, "%Y-%m-%d %H:%M:%S"),
    fit_finished_at = format(t_end,   "%Y-%m-%d %H:%M:%S")
  )

  timings <- list(
    wall_total_sec = wall_sec,
    wall_total_min = wall_sec / 60
  )

  list(
    meta        = meta,
    timings     = timings,
    diagnostics = diag,
    rt_df       = rt_df,
    fit         = fit
  )
}

# -----------------------------------------------------------------------------
# Save a fit result.
#
# Always writes:
#   <path with __summary.rds>  -> compact result list WITHOUT fit (~15 KB)
#                                 (meta + timings + diagnostics + rt_df)
# Writes only if save_full = TRUE:
#   <path>                     -> full result list WITH $fit (~5 MB)
#                                 Needed only by epinow2_inspect_fit.R for
#                                 digging into posterior draws / Stan diags.
#
# Default is summary-only: a full sweep (s1 + s2, every week, gp + grw) is
# ~100 fits, and keeping the stanfit on every one is ~500 MB which is painful
# to rsync/push.  Opt in via save_full = TRUE when you actively need to
# inspect a specific fit's posterior.
# -----------------------------------------------------------------------------
save_fit_result <- function(result, path, save_full = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  compact <- result
  compact$fit <- NULL
  summary_path <- sub("\\.rds$", "__summary.rds", path)
  saveRDS(compact, summary_path)

  if (isTRUE(save_full)) saveRDS(result, path)

  invisible(list(full = if (isTRUE(save_full)) path else NA_character_,
                 summary = summary_path))
}
