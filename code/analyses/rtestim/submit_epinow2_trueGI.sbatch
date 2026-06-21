#!/bin/bash
#SBATCH --job-name=en2_trueGI
#SBATCH --time=03:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-7
#SBATCH --output=logs/en2_%A_%a.out
#SBATCH --error=logs/en2_%A_%a.err
set -euo pipefail
mkdir -p logs

# Force the same R as the install (broken ~skel dotfiles make the bare batch R
# unpredictable). Module + explicit user library where EpiNow2 was installed.
module load R/4.6.0
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library-ubuntu-24.04/4.6"

FIT="$PWD/code/analyses/rtestim/epinow2_fit_one.R"

export MCMC_CHAINS=4 MCMC_WARMUP=1500 MCMC_SAMPLES=2000 MCMC_ADAPT_DELTA=0.95 MCMC_TREEDEPTH=12
export GP_KERNEL=matern GP_MATERN_ORDER=1.5 GP_ALPHA_SD=0.05 GP_LS_DIST=lognormal GP_LS_MEAN=10 GP_LS_SD=3 GP_LS_MAX=30 GP_BASIS_PROP=0.2

# array task -> (model, scenario):  1-4 = GP on S1-S4 ; 5-7 = GRW on S1-S3
MODELS=( gp gp gp gp grw grw grw )
SCEN=(   1  2  3  4  1   2   3   )
i=$((SLURM_ARRAY_TASK_ID - 1))
export MODEL_TYPE=${MODELS[$i]} OUTPUT_TAG=${MODELS[$i]}
[ "$MODEL_TYPE" = grw ] && export RW_STEP=1

echo "array=$SLURM_ARRAY_TASK_ID model=$MODEL_TYPE scen=${SCEN[$i]} $(date)"
Rscript "$FIT" ${SCEN[$i]}
echo "Finished: $(date)"
