#!/bin/bash
#SBATCH --job-name=tisr_templates
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/tisr_build_%j.log

# Builds the two reusable tisr templates (365-day, 366-day) at DART's native
# 6.6M-cell resolution -- see build_tisr_templates.py for why two templates
# are exact rather than an approximation.
#
# --exclusive on the `compute` partition reserves the ENTIRE node (256 cores,
# 257 GB RAM) for this job alone -- no other user's job can land on the same
# cores. Two earlier attempts hit real problems this avoids:
#   1. Ran directly on the login node: shared with ~90 other interactive
#      users, load ~27-29 already -- parallelism bought nothing, 24 workers
#      were as slow as 1.
#   2. Ran on SLURM's `shared` partition: CPUs were genuinely reserved, but
#      that specific node was 248/256 cores already allocated to OTHER
#      users' jobs, and the partition's fixed 940 MB/CPU ratio led to an
#      actual OUT_OF_MEMORY kill.
# --cpus-per-task=128 sizes the worker pool (see N_WORKERS in
# build_tisr_templates.py, capped at 24). --mem=0 is REQUIRED alongside
# --exclusive: without it, SLURM still applies its default per-CPU memory
# accounting underneath (940 MB/CPU on `compute`) even on an exclusive node,
# capping this job at ~128*940MB =~ 117 GB instead of the real 257 GB. That
# silent cap is what caused two runs to hang partway through: the first wave
# of ~24 chunks fit, but the second wave (each worker's 2nd task) pushed
# aggregate memory over ~117 GB, a worker got OOM-killed, and
# multiprocessing.Pool waited forever for a result that would never arrive.
# --mem=0 means "give this job all the memory on the node" in SLURM.
#
# Submit from a login node:  sbatch build_tisr_templates.sbatch
# Check status:              squeue --me
# Follow output:              tail -f tisr_build_<jobid>.log

set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

# Each process's numpy also tries to spin up its own OpenBLAS thread pool
# (128 threads by default on this node) purely on import, regardless of
# whether any BLAS routine is ever called. With N_WORKERS processes doing
# that simultaneously, total requested threads can exceed the account's
# RLIMIT_NPROC (2048), some pthread_create() calls fail, and OpenBLAS is
# left half-initialised -- confirmed directly: plain subprocesses (no
# Python multiprocessing involved at all) printed "pthread_create failed
# ... Resource temporarily unavailable" and several were later killed by
# the kernel. These three lines stop that thread pool from ever spinning up.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

python3 -u scripts/preprocessing/build_tisr_templates.py
