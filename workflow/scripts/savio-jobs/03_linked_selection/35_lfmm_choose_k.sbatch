#!/bin/bash
#SBATCH --job-name=35_lfmm_choose_k
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
# One task: reads the broken-stick table and every lambda sweep, and writes
# lfmm_k_chosen.csv, the single K per generation that 33 fits. Cheap (seconds).
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:20:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

set -euo pipefail
export TMPDIR=/global/scratch/users/milesroberts/tmp
module load r

ROOT=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2
JOBS=$ROOT/workflow/scripts/savio-jobs/03_linked_selection
OUT=$ROOT/results/$(date +%F)

# the newest folder that holds the broken-stick table, so a run started yesterday
# and read today still resolves
RES=$OUT
if [[ ! -s "$RES/lfmm_k_selection.csv" ]]; then
  RES=$(dirname "$(ls -t $ROOT/results/*/lfmm_k_selection.csv 2>/dev/null | head -1)")
fi
if [[ ! -s "$RES/lfmm_k_selection.csv" ]]; then
  echo "ERROR: no lfmm_k_selection.csv found (run 32_lfmm_k_selection.sbatch first)" >&2
  exit 1
fi

SWEEPS=$(ls -1 $RES/lfmm_k_sweep_gen*.csv 2>/dev/null || true)
mkdir -p "$OUT"

echo "broken-stick table: $RES/lfmm_k_selection.csv"
echo "sweeps: ${SWEEPS:-none}"

Rscript --no-save "$JOBS/35_lfmm_choose_k.R" \
  "$RES/lfmm_k_selection.csv" "$OUT/lfmm_k_chosen.csv" $SWEEPS
