#!/bin/bash
#SBATCH --job-name=01_merge_arch3v2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%j.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

# 16 CPUs -> 64 GB under savio4_htc's 4 GB/core rule. Peak is generation 1:
# 325 pools x 6,269,335 SNP records as doubles (16.3 GB) plus the non-SNP
# accumulator (5.9 GB), with a 200k-record block buffer during assembly.
# I/O: ~54 GB of float16 matrices read, ~60 GB of CSV written, ~30 GB of
# temporary float32 block files (removed on success).
module load r

Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/01_merge_arch3v2.R
