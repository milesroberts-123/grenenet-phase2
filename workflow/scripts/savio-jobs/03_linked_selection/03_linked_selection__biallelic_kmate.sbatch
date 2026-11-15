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

# Streams the fixrm pool matrix in 50k-row chunks, keeping only rows whose
# position carries exactly one record. Peak memory is one chunk (~360 MB), so
# the run is I/O bound: ~57 GB in, a little less out. 12 h covers it with margin.
module load r

Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/03_linked_selection__biallelic_kmate.R
