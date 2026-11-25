#!/bin/bash
#SBATCH --job-name=30_vw_subsets_s
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
# s for the nested 10k / 50k / 100k sub-panels, one-generation windows only.
# Needs the asreml seat (library(Vw)).
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

set -euo pipefail
export TMPDIR=/global/scratch/users/milesroberts/tmp
mkdir -p "$TMPDIR"
module load r

ROOT=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2
CACHE=$TMPDIR/vw_founder_chr3_v2.tsv
PANEL=$ROOT/data/vw_panel_100k_chr3.csv
OUT=$ROOT/data/vw_panel_subsets_s.csv

if [[ ! -s "$CACHE" ]]; then
  echo "ERROR: founder cache missing: $CACHE" >&2
  exit 1
fi
if [[ ! -s "$PANEL" ]]; then
  echo "ERROR: panel missing: $PANEL" >&2
  exit 1
fi

Rscript --no-save "$ROOT/workflow/scripts/savio-jobs/03_linked_selection/30_vw_subsets_s.R" \
  "$CACHE" "$PANEL" "$OUT"
