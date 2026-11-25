#!/bin/bash
#SBATCH --job-name=28_vw_site_va
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-29%1
# One task per site. Throttled to one task at a time: library(Vw) attaches asreml
# and the trial has a single seat, so parallel tasks would fail at load.
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A_%a.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

set -euo pipefail
export TMPDIR=/global/scratch/users/milesroberts/tmp
mkdir -p "$TMPDIR"
module load r

ROOT=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2
PANEL=$ROOT/data/vw_panel_100k_chr3.csv
STAB=$ROOT/data/vw_panel_100k_chr3_s.csv
OUT=$ROOT/results/$(date +%F)
mkdir -p "$OUT"

if [[ ! -s "$PANEL" ]]; then
  echo "ERROR: panel missing: $PANEL (run 27_vw_panel_extract.sbatch first)" >&2
  exit 1
fi

Rscript --no-save "$ROOT/workflow/scripts/savio-jobs/03_linked_selection/28_vw_site_va.R" \
  "$PANEL" "$STAB" "${SLURM_ARRAY_TASK_ID}" "$OUT"
