#!/bin/bash
#SBATCH --job-name=23_wright_s_estim_kmate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# one task per site. Each task collects that site's columns (~6.2M rows x
# generations x replicates), so it needs a wide memory margin and writes its own
# results/<today>/wright_s_estim_site<site>.csv; a failed site is resubmitted by
# index without touching the others.
#SBATCH --array=0-29
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A_%a.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

# R: 'module load r' uses Savio's default. If the code needs the same R as
# your OOD RStudio sessions, pin the version here (e.g. r/4.5.0).
module load r

# the task index is 0-based and maps onto the sorted site list inside the script.
# requires an rventures build that has wright_s_estim(components = )
Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/23_wright_s_estim_kmate.R "${SLURM_ARRAY_TASK_ID}"
