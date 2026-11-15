#!/bin/bash
#SBATCH --job-name=03_fix_site_labels
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%j.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

# Rewrites only the header line of the biallelic pool table: 1_4_* -> 5_4_* and
# 2_4_* -> 4_4_*. Single digit to single digit, so the line keeps its byte length
# and the 32 GB of data below it is never read. Seconds, not minutes; 1 CPU.
# Idempotent - a second run reports "nothing to rename" and exits.
module load r

Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/03_linked_selection__fix_site_labels.R
