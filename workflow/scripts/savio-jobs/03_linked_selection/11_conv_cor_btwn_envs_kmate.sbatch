#!/bin/bash
#SBATCH --job-name=11_conv_cor_btwn_envs_kmate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# the biallelic table is 6.2M records, 2.5x the arch3 one these jobs were
# sized for; the default 4 GB/CPU (16 GB) OOM-killed 13 at site 5
#SBATCH --mem=64G
# 24 h was chosen when this job only had to cache and pair 31 sites; it now builds
# 89 cache files and iterates 435 site pairs. scontrol update cannot extend a
# running job on Savio (access denied), so the limit is set here instead.
#SBATCH --time=168:00:00
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%j.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

# R: 'module load r' uses Savio's default. If the code needs the same R as
# your OOD RStudio sessions, pin the version here (e.g. r/4.5.0).
module load r

Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/11_conv_cor_btwn_envs_kmate.R
