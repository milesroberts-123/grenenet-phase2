#!/bin/bash
#SBATCH --job-name=16_feder_kmate_summary
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# one task per site. The all-sites loop was OOM-killed at site 5 after 29.7 h
# (2026-11-16); per site the work is one arrow aggregation plus one Manhattan
# plot. Each task writes feder_t_test_kmate_summary_site<site>.csv and
# feder_t_test_kmate_site_<site>.png. No throttle: every task computes Sys.Date()
# at start, so a task that begins after midnight would write into the next day's
# folder. Requires 18_feder_split to have written data/feder_sites first.
#SBATCH --array=0-29
# a site's CSV is ~8 GB; fread holds it plus the aggregated table
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A_%a.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

module load r

# the task index is 0-based and maps onto the sorted site list inside the script
Rscript --no-save /global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/16_feder_kmate_summary.R "${SLURM_ARRAY_TASK_ID}"
