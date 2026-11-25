#!/bin/bash
#SBATCH --job-name=33_lfmm_kmate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-5
# One task per generation. The task id IS the generation (1..5), so a failed
# generation is resubmitted by index without recomputing the others.
# No ASReml seat needed: lfmm uses plain R linear algebra.
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A_%a.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

set -euo pipefail
export TMPDIR=/global/scratch/users/milesroberts/tmp
mkdir -p "$TMPDIR"
module load r

ROOT=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2
JOBS=$ROOT/workflow/scripts/savio-jobs/03_linked_selection
TRAYS=$ROOT/data/PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency_trays_randswitch_fixrm_biallelic.csv
SEEDMIX=$ROOT/data/PHASE_1_2_COMBINED_kmate_arch3v2_allele_frequency_seedmix_randswitch_fixrm_biallelic.csv
# climate_sites_yearly.rda exists only in the home clone, not the scratch copy
CLIMATE=/global/home/users/milesroberts/moi_lab_projects/grenenet-phase2/data/climate_sites_yearly.rda
GEN=${SLURM_ARRAY_TASK_ID}

# K comes from 32_lfmm_k_selection; find the newest dated table if today's is absent
KCSV=$ROOT/results/$(date +%F)/lfmm_k_selection.csv
if [[ ! -s "$KCSV" ]]; then
  KCSV=$(ls -t $ROOT/results/*/lfmm_k_selection.csv 2>/dev/null | head -1 || true)
fi
if [[ -z "${KCSV:-}" || ! -s "$KCSV" ]]; then
  echo "ERROR: no lfmm_k_selection.csv found (run 32_lfmm_k_selection.sbatch first)" >&2
  exit 1
fi

OUT=$ROOT/results/$(date +%F)
WORK=$TMPDIR/lfmm_k_work
mkdir -p "$OUT" "$WORK"

if [[ ! -s "$CLIMATE" ]]; then
  echo "ERROR: climate file missing: $CLIMATE" >&2
  exit 1
fi

Rscript --no-save "$JOBS/33_lfmm_kmate.R" \
  "$TRAYS" "$SEEDMIX" "$CLIMATE" "$KCSV" "$GEN" \
  "$OUT/lfmm_nonsnp_gen${GEN}.csv" "$WORK" "$JOBS/run_lfmm_nogif.R"
