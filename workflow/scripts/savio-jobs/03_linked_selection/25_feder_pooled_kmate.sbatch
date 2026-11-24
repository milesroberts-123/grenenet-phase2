#!/bin/bash
#SBATCH --job-name=25_feder_pooled_kmate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
# One global aggregation, so this cannot be a per-site array: every variant needs
# rows from all sites at once. Arrow streams the 235 GB table and only the ~6M
# per-variant moments reach R, so the cost is CSV parsing, not statistics.
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

# requires rventures (plot_manhattan) and arrow >= 24
module load r

# first argument: the results date holding feder_t_test_kmate.tsv, or "latest"
# (the default) to use the folder whose table was written last
Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/25_feder_pooled_kmate.R "${1:-latest}"
