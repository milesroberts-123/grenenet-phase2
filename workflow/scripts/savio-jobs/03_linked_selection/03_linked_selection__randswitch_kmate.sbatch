#!/bin/bash
#SBATCH --job-name=03_randswitch_kmate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%j.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

# Streams the 60 GB trays CSV and rewrites it flipped. Peak memory is one 50k-row
# chunk (50000 x 893 doubles ~ 360 MB), so 4 CPUs / 16 GB is ample; the run is
# I/O bound. 12 h covers ~120 GB of traffic with margin.
module load r

Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/03_linked_selection__randswitch_kmate.R
