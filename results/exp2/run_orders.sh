#!/bin/bash
# Experiment 2, stage 1: all orders on the CV origins (8 variants) and on the test origins (2 variants)
cd /home/config/R/Projects/Experiments/cARIMA
unset EXP2_SMOKE
# One thread per process: must be set before R starts (OpenBLAS reads it when it is loaded)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript R/10_all_orders.R cv > results/exp2/log_10_cv.txt 2>&1
Rscript R/10_all_orders.R test carima_logZ_nd carima_logZ_cgn_nd > results/exp2/log_10_test.txt 2>&1
Rscript R/11_select_orders_cv.R > results/exp2/log_11.txt 2>&1
touch results/exp2/stage1.done
