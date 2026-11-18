#!/bin/bash
#SBATCH --job-name=18_feder_split
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
# Splits the 247 GB per-variant Feder table into one CSV per site, streaming, so
# the per-site summary array can read its own file with fread. Converting the
# table to parquet instead hit the 64 GB limit (arrow's pool cannot be drained in
# R 0.49/arrow 24). Single task: one sequential pass over the file.
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%j.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

module load r

# pass the results folder that holds the current feder table
Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/18_feder_split.R 2026-11-15
