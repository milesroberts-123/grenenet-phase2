#!/bin/bash
#SBATCH --job-name=32_lfmm_k_selection
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
# One task: extracts the non-SNP records once, then computes the covariance of the
# centred dp matrix for every generation and reports the broken-stick K for each.
# No lfmm dependency (eigen only), so it does not need the ASReml licence.
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
TRAYS=$ROOT/data/PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency_trays_randswitch_fixrm_biallelic.csv
SEEDMIX=$ROOT/data/PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency_seedmix_randswitch_fixrm_biallelic.csv
OUT=$ROOT/results/$(date +%F)/lfmm_k_selection.csv
WORK=$TMPDIR/lfmm_k_work
mkdir -p "$(dirname "$OUT")" "$WORK"

if [[ ! -s "$TRAYS" ]]; then
  echo "ERROR: trays table missing: $TRAYS" >&2
  exit 1
fi
if [[ ! -s "$SEEDMIX" ]]; then
  echo "ERROR: seedmix table missing: $SEEDMIX" >&2
  exit 1
fi

Rscript --no-save "$ROOT/workflow/scripts/savio-jobs/03_linked_selection/32_lfmm_k_selection.R" \
  "$TRAYS" "$SEEDMIX" "$OUT" "$WORK"
