#!/bin/bash
#SBATCH --job-name=03_linked_selection_feder_summary
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%j.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

# 8 CPUs -> 32 GB under savio4_htc's 4 GB/core rule. The 2026-11-06 feder run was
# OOM-killed at 4 CPUs / 16 GB while collecting site 5's 120M raw rows; this job
# aggregates inside arrow instead, and the extra headroom covers the per-site
# manhattan plot of ~2.5M variants.
module load r

Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/03_linked_selection__feder_kmate_summary.R
