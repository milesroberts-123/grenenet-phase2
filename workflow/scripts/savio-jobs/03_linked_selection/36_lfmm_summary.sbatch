#!/bin/bash
#SBATCH --job-name=36_lfmm_summary
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
# One task: concatenates the per-generation fits, computes BH within each
# generation, and draws QQ and Manhattan figures. No lfmm dependency.
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
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
JOBS=$ROOT/workflow/scripts/savio-jobs/03_linked_selection

# today's results if the fits are there, otherwise the newest folder that has them
OUT=$ROOT/results/$(date +%F)
if ! ls "$OUT"/lfmm_nonsnp_gen*.csv > /dev/null 2>&1; then
  OUT=$(dirname "$(ls -t $ROOT/results/*/lfmm_nonsnp_gen*.csv 2>/dev/null | head -1)")
fi
if [[ ! -d "${OUT:-}" ]]; then
  echo "ERROR: no lfmm_nonsnp_gen*.csv found (run 33_lfmm_kmate.sbatch first)" >&2
  exit 1
fi
echo "results directory: $OUT"

Rscript --no-save "$JOBS/36_lfmm_summary.R" "$OUT"
