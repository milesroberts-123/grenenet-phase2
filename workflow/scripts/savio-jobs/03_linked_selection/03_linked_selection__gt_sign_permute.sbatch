#!/bin/bash
#SBATCH --job-name=03_linked_selection_gt_sign_permute
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesdroberts@gmail.com
#SBATCH --output=logs/03_linked_selection_gt_sign_permute_%j.out
#SBATCH --chdir=/global/scratch/users/milesroberts/grenenet-phase2/resources

# R: 'module load r' uses Savio's default. If the code needs the same R as
# your OOD RStudio sessions, pin the version here (e.g. r/4.5.0).
module load r

Rscript --no-save /global/scratch/users/milesroberts/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/03_linked_selection__gt_sign_permute.R