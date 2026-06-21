## ================================================================================
## rtestim_paper_Rt.R
##
## Ground-truth R_t scenarios from the rtestim paper
## (Liu, Zhang, McDonald, PLoS Comp Bio 2024, doi:10.1371/journal.pcbi.1012324).
## Values match the authors' dat/data_example.RDS to machine precision.
##
## Usage: source("Rt_scripts/sim_experiments/rtestim_paper_Rt.R")
##        Rt_list <- true_Rt_list(n = 300)
## ================================================================================

rtestim_scenario1 <- function(n = 300) {           # piecewise constant
  t <- seq_len(n)
  ifelse(t <= 120, 2, 0.8)
}

rtestim_scenario2 <- function(n = 300) {           # piecewise exponential
  # Note: intentional discontinuity at t=100, matches authors' data.
  full <- c(exp(0.01 * 1:100), exp(0.5 - 0.005 * 1:200))
  full[seq_len(n)]
}

rtestim_scenario3 <- function(n = 300) {           # piecewise linear (four segments)
  t <- seq_len(n)
  seg <- function(a, b, v) v * (t >= a & t < b)
  seg(1,   76,     2.5 - (0.5 / 74) * (t - 1)) +
    seg(76,  151,  0.8 - (0.2 / 74) * (t - 76)) +
    seg(151, 226,  1.7 + (0.3 / 74) * (t - 151)) +
    seg(226, n + 1, 0.9 - (0.4 / 74) * (t - 226))
}

rtestim_scenario4 <- function(n = 300) {           # periodic
  tt <- seq(0, 10, length.out = n)
  0.2 * ((sin(pi * tt / 12) + 1) +
         (2 * sin(5 * pi * tt / 12) + 2) +
         (3 * sin(5 * pi * tt / 6) + 3))
}

true_Rt_list <- function(n = 300) {
  list(
    `Scenario 1: piecewise constant`    = rtestim_scenario1(n),
    `Scenario 2: piecewise exponential` = rtestim_scenario2(n),
    `Scenario 3: piecewise linear`      = rtestim_scenario3(n),
    `Scenario 4: periodic`              = rtestim_scenario4(n)
  )
}
