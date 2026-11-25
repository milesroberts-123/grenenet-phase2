#!/bin/bash
#SBATCH --job-name=27_vw_panel_extract
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
# One task: freezes the 100k chr3 panel for the per-site Va runs. Streams the
# 42 GB trays table once with awk, and recomputes s for ngen 1..5 on the same
# window 26_vw_panel_compare used. Needs the asreml seat (library(Vw)).
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=6:00:00
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
if [[ ! -s "$CACHE" ]]; then
  echo "ERROR: founder cache missing: $CACHE (rerun 26_vw_panel_compare.sbatch)" >&2
  exit 1
fi

Rscript --no-save "$ROOT/workflow/scripts/savio-jobs/03_linked_selection/27_vw_panel_extract.R" \
  "$CACHE" 3 100000 "$ROOT"
