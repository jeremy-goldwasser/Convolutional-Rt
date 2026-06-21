## ================================================================================
## rtestim_paper_scenarios_seir.r
##
## SEIR simulator for the rtestim-paper scenarios.
##
##   SEIR simulation with explicit discrete-time gamma delays
##         - E -> I  (latent period)
##         - I -> R  (infectious duration)
##         - E -> Y  (exposure-to-report delay)
##       Latent mean/sd chosen so the INDUCED generation-interval mean/sd match
##       a fake / synthetic epidemic's SI (mean 8.4, sd 3.8), given a fixed
##       choice of infectious distribution (Gamma, shape 5, mean = 3).
##       Simulation uses the active-infectious kernel K[d] = P(L <= d, L + D > d):
##           Lambda_t = R_t * sum_{s<t} inc[s] * K[t-s] / mean_D
##           inc[t]  ~ Poisson(Lambda_t)
##       Observed case reports are obtained by dispatching each infection forward
##       by an individually sampled E->Y delay.
##
## Outputs in outputs/rtestim/:
##   - data/sim_combined.rds              one row per (scenario, si_type, time) with columns
##                                        seir_incidence, seir_reports, Rt
##   - data/seir_params.rds               chosen SEIR delay params + achieved induced GI
## ================================================================================

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep('^--file=', args, value = TRUE)
  if (length(m)) return(normalizePath(dirname(sub('^--file=', '', m[1]))))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(dirname(of)))
  }
  normalizePath(getwd())
})
source(file.path(script_dir, "..", "..", "_paths.R"))
rtestim_data_dir    <- file.path(data_dir,    "sim", "rtestim")
rtestim_results_dir <- file.path(results_dir, "sim", "rtestim")
rtestim_figures_dir <- file.path(figures_dir, "sim", "rtestim")
source(file.path(script_dir, "scenarios_Rt.R"))

# Simulated counts go to Experimental_data.
out_data_dir <- rtestim_data_dir
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------------
# 1. Discrete gamma PMF on {1, 2, ..., k} (midpoint convention -> unbiased mean)
# --------------------------------------------------------------------------------
# p_d = P(X in [d-0.5, d+0.5]) for d = 1,...,k. Any mass on [0, 0.5] is folded
# into d=1 so the support starts at 1 (appropriate for delays >= 1 day).
discrete_gamma <- function(mean_x, sd_x, k = 60) {
  shape <- (mean_x / sd_x)^2
  rate  <- mean_x / sd_x^2
  edges <- c(0, seq(1.5, k + 0.5, by = 1))  # cell right-edges: 0.5 rolled into d=1
  F_edges <- pgamma(edges, shape, rate)
  p <- diff(c(0, F_edges))                  # length k+1; first cell is [0, 0.5]? no
  # Build: d=1 cell = [0, 1.5); d=2 cell = [1.5, 2.5); ...; d=k cell = [k-0.5, k+0.5).
  p <- numeric(k)
  p[1] <- pgamma(1.5, shape, rate)
  if (k >= 2) p[2:k] <- pgamma((2:k) + 0.5, shape, rate) -
                       pgamma((2:k) - 0.5, shape, rate)
  p / sum(p)
}

# --------------------------------------------------------------------------------
# 2. Pick SEIR latent params to match a target GI mean/sd given (mean_D, shape_D)
# --------------------------------------------------------------------------------
# Discrete-time derivation for the convention "primary infectious on days
# {L, L+1, ..., L+D-1} and offspring-day uniform over that window":
#     GI = L + U,  U | D ~ Uniform{0, 1, ..., D-1},  L indep D
#     E[U]   = (mu_D - 1)/2 + sigma_D^2 / (2 mu_D)
#     Var[U] = sigma_D^2 / 3 + (mu_D^2 - 1) / 12
# With E[GI] = E[L] + E[U] and Var[GI] = Var[L] + Var[U], solve for (mu_L, sd_L).
seir_latent_params <- function(mean_si, sd_si, mean_D, shape_D) {
  var_D  <- mean_D^2 / shape_D
  mean_U <- (mean_D - 1) / 2 + var_D / (2 * mean_D)
  var_U  <- var_D / 3 + (mean_D^2 - 1) / 12
  mean_L <- mean_si - mean_U
  var_L  <- sd_si^2 - var_U
  if (var_L <= 0) stop("Infeasible: infectious gamma induces too much GI variance; raise shape_D.")
  list(mean_L = mean_L, sd_L = sqrt(var_L),
       mean_D = mean_D, sd_D = sqrt(var_D))
}

seir_params <- list(
  fake = c(seir_latent_params(8.4, 3.8, mean_D = 5, shape_D = 5),
           list(mean_EY = 7, sd_EY = 3))
)

# --------------------------------------------------------------------------------
# 3. Active-infectious kernel and induced GI (discrete)
# --------------------------------------------------------------------------------
# Convention: exposure on day 0 => infectious on days L, L+1, ..., L+D-1.
# K[d] = P(active at lag d >= 1) = sum_{l=1..d} pi_L[l] * P(D >= d - l + 1).
build_K <- function(pi_L, pi_D, max_d) {
  L <- length(pi_L)
  # survD[k+1] = P(D > k) = P(D >= k+1). survD[1] = 1.
  survD <- c(1, 1 - cumsum(pi_D))
  survD <- c(survD, rep(0, max_d + 1))
  K <- numeric(max_d)
  for (d in seq_len(max_d)) {
    lmax <- min(L, d)
    K[d] <- sum(pi_L[1:lmax] * survD[(d - 1:lmax) + 1])
  }
  K
}

induced_gi <- function(pi_L, pi_D, mean_D, max_d = 120) {
  K <- build_K(pi_L, pi_D, max_d)
  pmf <- K / mean_D
  pmf <- pmf / sum(pmf)        # renormalize over finite support
  d <- seq_along(pmf)
  m <- sum(d * pmf); s <- sqrt(sum(d^2 * pmf) - m^2)
  list(pmf = pmf, mean = m, sd = s)
}

# Monte-Carlo empirical GI: for each of n_primaries, sample L and D, then for each
# infectious day of that primary, draw Poisson(1) secondary infections at that
# absolute day; collect (secondary_exposure_day - primary_exposure_day) = GI sample.
empirical_gi <- function(pi_L, pi_D, n_primaries = 20000, rate_per_day = 1,
                         seed = 42) {
  set.seed(seed)
  L <- sample(seq_along(pi_L), n_primaries, replace = TRUE, prob = pi_L)
  D <- sample(seq_along(pi_D), n_primaries, replace = TRUE, prob = pi_D)
  # Build a flat vector of infectious-day offsets (one entry per primary-day)
  inf_days <- sequence(D, from = L)            # recycles: L_i, L_i+1, ..., L_i+D_i-1
  offspring <- rpois(length(inf_days), rate_per_day)
  gi_samples <- rep(inf_days, offspring)
  pmf <- tabulate(gi_samples, nbins = max(gi_samples))
  pmf <- pmf / sum(pmf)
  list(samples = gi_samples,
       mean = mean(gi_samples),
       sd   = sd(gi_samples),
       n    = length(gi_samples),
       pmf  = pmf)
}

# --------------------------------------------------------------------------------
# 4. Simulators
# --------------------------------------------------------------------------------
simulate_seir <- function(Rt, pi_L, pi_D, pi_EY, y1 = 2, warmup_days = 0L) {
  n <- length(Rt)
  mean_D <- sum(pi_D * seq_along(pi_D))
  N <- n + warmup_days
  K <- build_K(pi_L, pi_D, max_d = N)
  inc <- numeric(N)
  if (warmup_days > 0L) {
    # Pre-seed warmup_days of constant y1 exposures so that the active pool
    # is saturated at the scenario start. Without this, a single-day seed
    # dies out when mean_L is long (K[1] ~ 0 => no spread for ~mean_L days).
    inc[seq_len(warmup_days)] <- y1
    t_start <- warmup_days + 1L
  } else {
    inc[1] <- y1
    t_start <- 2L
  }
  for (t in t_start:N) {
    win    <- seq_len(t - 1L)
    active <- sum(inc[win] * K[t - win])
    Rt_t   <- Rt[t - warmup_days]
    inc[t] <- rpois(1, Rt_t * active / mean_D)
  }
  # Reports: per-infection report delay ~ pi_EY
  reports <- numeric(N)
  tot <- sum(inc)
  if (tot > 0) {
    expo_days <- rep(seq_len(N), inc)
    delays    <- sample(seq_along(pi_EY), tot, replace = TRUE, prob = pi_EY)
    rep_days  <- expo_days + delays
    keep      <- rep_days <= N
    reports   <- tabulate(rep_days[keep], nbins = N)
  }
  # Trim warmup so output aligns with Rt[1..n]
  keep_rng <- (warmup_days + 1L):N
  list(inc = inc[keep_rng], reports = reports[keep_rng])
}

# --------------------------------------------------------------------------------
# 5. Run ONE replicate per scenario x SI (seeds chosen so the SEIR model is
#    not stuck on a pure-extinction draw at the low-R start of Scenario 2)
# --------------------------------------------------------------------------------
Rt_scenarios <- true_Rt_list(300)

seeds <- list(
  fake = list(`Scenario 1: piecewise constant`    = 1,
              `Scenario 2: piecewise exponential` = 5,
              `Scenario 3: piecewise linear`      = 1,
              `Scenario 4: periodic`              = 1)
)

# Per-scenario initial seed. All four scenarios now start at y1=10 with a
# 30-day warmup. The warmup saturates the active-infectious pool so day 1 of
# the visible window sits at ~10 reports without the latent-period burn-in.
init_by_sc <- list(
  `Scenario 1: piecewise constant`    = list(y1 = 10L, warmup_days = 30L),
  `Scenario 2: piecewise exponential` = list(y1 = 10L, warmup_days = 30L),
  `Scenario 3: piecewise linear`      = list(y1 = 10L, warmup_days = 30L),
  `Scenario 4: periodic`              = list(y1 = 10L, warmup_days = 30L)
)

gi_rows  <- list()
out_rows <- list()

for (sit in names(seir_params)) {
  p     <- seir_params[[sit]]
  pi_L  <- discrete_gamma(p$mean_L,  p$sd_L,  k = 60)
  pi_D  <- discrete_gamma(p$mean_D,  p$sd_D,  k = 40)
  pi_EY <- discrete_gamma(p$mean_EY, p$sd_EY, k = 30)
  mean_D_discrete <- sum(pi_D * seq_along(pi_D))

  gi_th  <- induced_gi(pi_L, pi_D, mean_D_discrete)
  gi_emp <- empirical_gi(pi_L, pi_D, n_primaries = 20000)

  gi_rows[[sit]] <- data.frame(
    si_type = sit,
    target_mean  = 8.4,
    target_sd    = 3.8,
    theor_mean = gi_th$mean,  theor_sd = gi_th$sd,
    emp_mean   = gi_emp$mean, emp_sd   = gi_emp$sd, emp_n = gi_emp$n,
    mean_L = p$mean_L, sd_L = p$sd_L,
    mean_D = p$mean_D, sd_D = p$sd_D,
    mean_EY = p$mean_EY, sd_EY = p$sd_EY
  )

  for (sc in names(Rt_scenarios)) {
    Rt     <- Rt_scenarios[[sc]]
    y1     <- init_by_sc[[sc]]$y1
    w_days <- init_by_sc[[sc]]$warmup_days
    set.seed(seeds[[sit]][[sc]]); seir <- simulate_seir(Rt, pi_L, pi_D, pi_EY, y1 = y1,
                                                        warmup_days = w_days)
    out_rows[[paste(sit, sc, sep = "__")]] <- data.frame(
      time = seq_along(Rt), Rt = Rt,
      scenario = sc, si_type = sit,
      seir_incidence = seir$inc,
      seir_reports   = seir$reports
    )
  }
}

sim_combined <- do.call(rbind, out_rows); rownames(sim_combined) <- NULL
gi_table     <- do.call(rbind, gi_rows);  rownames(gi_table) <- NULL

saveRDS(sim_combined, file.path(out_data_dir, "sim_combined.rds"))
saveRDS(gi_table,     file.path(out_data_dir, "seir_params.rds"))

cat("\nSEIR delay params with theoretical vs empirical GI (midpoint discretization):\n")
print(gi_table[, c("si_type","target_mean","target_sd",
                   "theor_mean","theor_sd","emp_mean","emp_sd","emp_n",
                   "mean_L","sd_L","mean_D","sd_D","mean_EY","sd_EY")],
      row.names = FALSE, digits = 3)

cat("\nTotals and peaks per scenario x SI (single replicate):\n")
tot <- aggregate(cbind(seir_incidence, seir_reports) ~
                 scenario + si_type, data = sim_combined, FUN = sum)
pk  <- aggregate(seir_incidence ~ scenario + si_type, data = sim_combined, FUN = max)
names(pk)[3] <- "seir_peak"
print(merge(tot, pk, by = c("scenario","si_type")), row.names = FALSE)

cat("\nSaved:\n",
    "  ", file.path(out_data_dir, "sim_combined.rds"), "\n",
    "  ", file.path(out_data_dir, "seir_params.rds"), "\n", sep = "")

