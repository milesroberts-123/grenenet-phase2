#!/bin/bash
#SBATCH --job-name=24_wright_s_manhattan_kmate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
# One task per site: reads that site's wright_s_estim table, keeps the longest
# segment of each variant, clips the |s| tail at the 99th percentile and renders.
# The wall clock is the slowest site. No throttle, so every task stamps the same
# Sys.Date() results folder.
#SBATCH --array=0-29
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A_%a.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

# requires rventures >= 3.5.0 (wright_s_estim(components =), plot_manhattan(plot_effects =))
module load r

# first argument: the results date whose wright_s_estim_site<site>.csv get plotted
Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/24_wright_s_manhattan_kmate.R "${1:-$(date +%F)}" "${SLURM_ARRAY_TASK_ID}"
