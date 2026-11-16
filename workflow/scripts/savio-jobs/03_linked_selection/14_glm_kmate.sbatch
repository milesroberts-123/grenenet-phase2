#!/bin/bash
#SBATCH --job-name=14_glm_kmate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# one task per site. Site 1 alone takes ~7 h and there are 30 sites, so a single
# task cannot finish inside the 168 h walltime. Each task writes its own
# results/<today>/glm_quasibinomial_coef_site<site>.txt; a failed site costs one
# site rather than the whole run, and a site can be resubmitted by index.
#SBATCH --array=0-29
# the biallelic table is 6.2M records and each task holds one site's columns
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A_%a.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

# R: 'module load r' uses Savio's default. If the code needs the same R as
# your OOD RStudio sessions, pin the version here (e.g. r/4.5.0).
module load r

# the task index is 0-based and maps onto the sorted site list inside the script
Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/14_glm_kmate.R "${SLURM_ARRAY_TASK_ID}"
