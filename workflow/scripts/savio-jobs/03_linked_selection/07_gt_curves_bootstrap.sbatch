#!/bin/bash
#SBATCH --job-name=07_gt_curves_bootstrap
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# one task per site. the block bootstrap draws 500 times per replicate and there are 30
# sites, so a single all-sites task needs ~100 h even though the slowest site
# needs ~4 h. No throttle (%N): every task computes Sys.Date() at start, so a task
# that only begins after midnight writes into the next day's results folder and
# splits the output. Each writes gt/var_cov files for its own site, and
# 17_concat_gt_resampling joins them afterwards.
#SBATCH --array=0-29
# the biallelic table is 6.2M records and each task holds one site's columns
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A_%a.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

module load r

# the task index is 0-based and maps onto the sorted site list inside the script
Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/07_gt_curves_bootstrap.R "${SLURM_ARRAY_TASK_ID}"
