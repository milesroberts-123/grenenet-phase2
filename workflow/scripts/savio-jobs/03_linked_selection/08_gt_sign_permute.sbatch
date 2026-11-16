#!/bin/bash
#SBATCH --job-name=08_gt_sign_permute
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# one task per site. 500 sign permutations per replicate and there are 30
# sites, so a single all-sites task needs ~100 h even though the slowest site
# needs ~4 h. Tasks are throttled to 10 concurrent so the array does not take the
# whole partition; each writes gt/var_cov files for its own site, and
# 17_concat_gt_resampling joins them afterwards.
#SBATCH --array=0-29%10
# the biallelic table is 6.2M records and each task holds one site's columns
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A_%a.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

module load r

# the task index is 0-based and maps onto the sorted site list inside the script
Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/08_gt_sign_permute.R "${SLURM_ARRAY_TASK_ID}"
