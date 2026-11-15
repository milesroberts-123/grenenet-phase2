#!/bin/bash
#SBATCH --job-name=07_gt_curves_bootstrap
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# the biallelic table is 6.2M records, 2.5x the arch3 one these jobs were
# sized for; the default 4 GB/CPU (16 GB) OOM-killed 13 at site 5
#SBATCH --mem=64G
#SBATCH --time=336:00:00
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%j.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

# R: '# 500 bootstrap draws at 8.5M records project to ~258 h (the 2.5M-record runs took
# ~76 h). savio4_htc reports no partition or QoS wall limit, so this requests
# 336 h rather than the previous 168 h. If the site still enforces 7 days, the
# fallback is a per-site array with resume, not a longer wall clock.
module load r' uses Savio's default. If the code needs the same R as
# your OOD RStudio sessions, pin the version here (e.g. r/4.5.0).
module load r

Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/07_gt_curves_bootstrap.R
