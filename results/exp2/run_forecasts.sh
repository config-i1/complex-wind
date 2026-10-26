#!/bin/bash
# Experiment 2, stage 2: final forecasts with predictive distributions and the evaluation.
# Waits for stage 1 (R/10, R/11), then installs complex from the local clone (branch intervals),
# which is safe only when no R process has the package loaded.
cd /home/config/R/Projects/Experiments/cARIMA
unset EXP2_SMOKE
# One thread per process: must be set before R starts (OpenBLAS reads it when it is loaded)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
until [ -f results/exp2/stage1.done ]; do sleep 60; done
R CMD INSTALL /home/config/R/Projects/Packages/complex > results/exp2/log_install.txt 2>&1
Rscript R/12_forecasts_exp2.R > results/exp2/log_12.txt 2>&1
Rscript R/13_evaluate_exp2.R > results/exp2/log_13.txt 2>&1
touch results/exp2/stage2.done
