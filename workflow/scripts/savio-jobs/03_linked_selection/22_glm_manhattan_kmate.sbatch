#!/bin/bash
#SBATCH --job-name=22_glm_manhattan_kmate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
# One task per site: each reads four columns of its ~600 MB coefficient table,
# BH adjusts, and renders ~2.5 M points. A knit would take 25-35 minutes over all
# 30 sites; as an array the wall clock is the slowest site. No throttle, so every
# task stamps the same Sys.Date() results folder.
#SBATCH --array=0-29
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A_%a.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

module load r

# 2026-11-15 is the GLM run whose per-site tables get plotted
Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/22_glm_manhattan_kmate.R 2026-11-15 "${SLURM_ARRAY_TASK_ID}"
