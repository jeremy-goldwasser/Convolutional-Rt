# code/analyses/rtestim/regen.R
# Build girt MechRt artifacts for the rtestim figures, into results/sim/rtestim/seir/girt/data/:
#   retro_rt_estimates.rds  -- girt MechRt (smooth) + cached baselines (estimateR/EpiEstim/rtestim/EpiLPS)
#   fit_ours_jump_<tag>.rds -- girt MechRt(jumps) for S1-S3 (list with $rt)
# day is 1-based (= time) to match the truth/baseline convention in the plots.

sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/rtestim"
source(file.path(sd, "_common.R"))

sim <- readRDS(file.path(rt_data_dir, "sim_combined.rds")); sim$date <- rt_cfg$date0 + (sim$time - 1L)
dl  <- rt_delays("fake")
gtree <- file.path(rt_results_dir, "girt"); dir.create(gtree, recursive = TRUE, showWarnings = FALSE)
JUMPS <- list("Scenario 1: piecewise constant"=120L, "Scenario 2: piecewise exponential"=100L,
              "Scenario 3: piecewise linear"=c(75L,150L,225L))
tag_for <- function(sc) gsub(":? ", "_", sc)

mech_rows <- list()
for (sc in unique(sim$scenario)) {
  d <- sim[sim$scenario == sc & sim$si_type == "fake", ]; d <- d[order(d$time), ]
  y <- as.numeric(d$seir_reports); dates <- d$date
  gf <- fit_rtestim_girt(y, dates, dl)
  rt <- gf$rt; rt$method <- "MechRt"; rt$scenario <- sc; rt$si_type <- "fake"
  rt$day <- as.integer(rt$date - rt_cfg$date0) + 1L
  mech_rows[[length(mech_rows)+1L]] <- rt[, c("method","scenario","si_type","date","Rt_mean","Rt_lo","Rt_hi","day")]
  if (!is.null(JUMPS[[sc]])) {
    gj <- fit_rtestim_girt(y, dates, dl, jump_times = JUMPS[[sc]])
    rj <- gj$rt; rj$method <- "MechRt"; rj$scenario <- sc; rj$si_type <- "fake"
    saveRDS(list(rt = rj[, c("date","day","Rt_mean","Rt_lo","Rt_hi")]),
            file.path(gtree, sprintf("fit_ours_jump_fake__%s.rds", tag_for(sc))))
  }
}
mech_g <- do.call(rbind, mech_rows)

# bundle: cached baselines (from per-method fit files) + girt MechRt
data_dir_c <- rt_results_dir
base_keys <- c(estimateR="estimater", EpiEstim="epiestim", rtestim="rtestim", EpiLPS="epilps")
base <- do.call(rbind, lapply(names(base_keys), function(meth) do.call(rbind, lapply(unique(sim$scenario), function(sc) {
  x <- readRDS(file.path(data_dir_c, sprintf("fit_%s_fake__%s.rds", base_keys[[meth]], tag_for(sc))))
  rt <- if (is.data.frame(x)) x else x$rt
  data.frame(method = meth, scenario = sc, si_type = "fake", date = as.Date(rt$date),
             Rt_mean = rt$Rt_mean, Rt_lo = rt$Rt_lo, Rt_hi = rt$Rt_hi,
             day = as.integer(as.Date(rt$date) - rt_cfg$date0) + 1L)
}))))
bundle <- rbind(base[, names(mech_g)], mech_g)
bundle$method <- factor(bundle$method, levels = c("MechRt","estimateR","EpiEstim","rtestim","EpiLPS"))
saveRDS(bundle, file.path(gtree, "retro_rt_estimates.rds"))
cat("wrote girt retro_rt_estimates.rds + fit_ours_jump (S1-S3)\n")
cat("bundle methods:", paste(levels(droplevels(bundle$method)), collapse=", "), "\n")
