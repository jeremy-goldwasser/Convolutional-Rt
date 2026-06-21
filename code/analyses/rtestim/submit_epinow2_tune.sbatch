#!/bin/bash
# Tuning sweep: rescue the S1 + S3 EpiNow2 fits that didn't converge at the true
# GI under the loose 1500/0.95/12 settings (R-hat in the hundreds-thousands).
# 3 settings x {S1,S3} x {gp,grw} = 12 jobs. Each writes to its own tagged dir
# results/sim/rtestim/seir/data/epinow2_<model>_<cfg>/ so configs don't collide.
# Pick the winning config from R-hat + MAE-vs-truth, then rerun S2/S4 with it.
#SBATCH --job-name=en2_tune
#SBATCH --time=02:30:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-12
#SBATCH --output=logs/en2_tune_%A_%a.out
#SBATCH --error=logs/en2_tune_%A_%a.err
set -euo pipefail
mkdir -p logs

module load R/4.6.0
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library-ubuntu-24.04/4.6"
FIT="$PWD/code/analyses/rtestim/epinow2_fit_one.R"

# fixed (GP hyperparams as before; chains=4)
export MCMC_CHAINS=4
export GP_KERNEL=matern GP_MATERN_ORDER=1.5 GP_ALPHA_SD=0.05 GP_LS_DIST=lognormal GP_LS_MEAN=10 GP_LS_SD=3 GP_LS_MAX=30 GP_BASIS_PROP=0.2

# 3 configs x (S1,gp)(S3,gp)(S1,grw)(S3,grw).  c1=warmup-heavy, c2=tightness-heavy, c3=more-of-everything.
W=(  3000 3000 3000 3000   2000 2000 2000 2000   2500 2500 2500 2500 )
S=(  2000 2000 2000 2000   2000 2000 2000 2000   2500 2500 2500 2500 )
A=(  0.97 0.97 0.97 0.97   0.99 0.99 0.99 0.99   0.97 0.97 0.97 0.97 )
T=(  13   13   13   13     14   14   14   14     13   13   13   13   )
SCN=( 1   3    1    3       1    3    1    3       1    3    1    3   )
MDL=( gp  gp   grw  grw     gp   gp   grw  grw     gp   gp   grw  grw )
CFG=( c1  c1   c1   c1      c2   c2   c2   c2      c3   c3   c3   c3  )

i=$((SLURM_ARRAY_TASK_ID - 1))
export MCMC_WARMUP=${W[$i]} MCMC_SAMPLES=${S[$i]} MCMC_ADAPT_DELTA=${A[$i]} MCMC_TREEDEPTH=${T[$i]}
export MODEL_TYPE=${MDL[$i]} OUTPUT_TAG=${MDL[$i]}_${CFG[$i]}
[ "$MODEL_TYPE" = grw ] && export RW_STEP=1

echo "=== task=$SLURM_ARRAY_TASK_ID cfg=${CFG[$i]} model=$MODEL_TYPE scen=${SCN[$i]} W=$MCMC_WARMUP S=$MCMC_SAMPLES A=$MCMC_ADAPT_DELTA T=$MCMC_TREEDEPTH $(date) ==="
Rscript "$FIT" ${SCN[$i]}
echo "Finished: $(date)"
