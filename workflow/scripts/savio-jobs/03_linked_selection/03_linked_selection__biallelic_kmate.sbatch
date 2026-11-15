#!/bin/bash
#SBATCH --job-name=03_biallelic_kmate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%j.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

# Streams the fixrm pool matrix in 20k-row batches, keeping only rows whose
# position carries exactly one record. 16 CPUs -> 64 GB under savio4_htc's
# 4 GB/core rule: the 2026-11-14 run was OOM-killed at 4 CPUs / 16 GB (MaxRSS
# 16.4 GB) after ~5 minutes, because arrow retains freed buffers in its memory
# pool across batches and 893 columns per batch add up. The job now also calls
# ReleaseUnused() after every batch.
module load r

Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/03_linked_selection__biallelic_kmate.R
