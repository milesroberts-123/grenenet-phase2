#!/bin/bash
#SBATCH --job-name=31_vw_window_sizes
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-89
# 90 tasks = 3 panel sizes x 30 sites. Task index: size_idx = ID / 30,
# site_idx = ID %% 30. One-generation windows only (base -> gen 1, and every
# adjacent measured pair). No licence seat needed: Q() and lm() only.
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=3:00:00
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
STAB=$ROOT/data/vw_panel_subsets_s.csv
SIZES=(10000 50000 100000)

ID=${SLURM_ARRAY_TASK_ID}
SIZE_IDX=$(( ID / 30 ))
SITE_IDX=$(( ID % 30 ))
SIZE=${SIZES[$SIZE_IDX]}

OUT=$ROOT/results/$(date +%F)/window_sizes
mkdir -p "$OUT"

if [[ ! -s "$STAB" ]]; then
  echo "ERROR: subset s table missing: $STAB (run 30_vw_subsets_s.sbatch first)" >&2
  exit 1
fi

Rscript --no-save "$ROOT/workflow/scripts/savio-jobs/03_linked_selection/28_vw_site_va.R" \
  "$PANEL" "$STAB" "$SITE_IDX" "$OUT" onen "$SIZE"
