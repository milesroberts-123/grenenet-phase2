#!/bin/bash
#SBATCH --job-name=03_biallelic_kmate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%j.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

# Pure line filter: reads the fixrm pool matrix in 50k-line chunks and copies the
# kept lines verbatim, so memory stays under ~1 GB and no numeric parsing
# happens. Two earlier attempts failed here: readr reported EOF after ~6.9M of
# 7.56M rows and wrote a short table, and an arrow RecordBatchReader was
# OOM-killed at 16 GB while decoding and re-encoding 893 columns per batch.
# 4 CPUs / 16 GB is ample; the run is I/O bound (39 GB in, ~32 GB out).
module load r

Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/03_linked_selection__biallelic_kmate.R
