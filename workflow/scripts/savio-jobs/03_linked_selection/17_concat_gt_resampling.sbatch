#!/bin/bash
#SBATCH --job-name=17_concat_gt_resampling
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
# Joins the 30 per-site outputs of 07 and 08 into the combined tables the
# notebook reads. Submit it only after both arrays have finished; it refuses to
# write a combined table if any site file is missing or a site appears twice.
#
# Args: <out_date> [<search_date> ...]. The arrays stamp every task with
# Sys.Date() at start, so tasks that began after midnight wrote into the next
# day's folder and the files are spread over two dates. Pass every date to read;
# the combined tables are written into the first. Leave the defaults to read only
# today.
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%j.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

module load r

# 2026-11-15 arrays: 23 sites in 2026-11-15, 7 tasks started on the 16th
Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/17_concat_gt_resampling.R \
  2026-11-15 2026-11-15 2026-11-16
