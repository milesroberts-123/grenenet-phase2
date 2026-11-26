#!/bin/bash
#SBATCH --job-name=34_lfmm_ksweep
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=2-5
# One task per generation, the task id being the generation. Generation 1 has
# already been swept (K=5..30 all gave lambda 1.02 and a flat null), so this
# array fills in the remaining generations with the same criterion and K list.
# No ASReml seat: lfmm uses plain linear algebra.
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=6:00:00
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
CLIMATE=/global/home/users/milesroberts/moi_lab_projects/grenenet-phase2/data/climate_sites_yearly.rda
GEN=${SLURM_ARRAY_TASK_ID}
OUT=$ROOT/results/$(date +%F)/lfmm_k_sweep_gen${GEN}.csv
mkdir -p "$(dirname "$OUT")"

if [[ ! -s "$CLIMATE" ]]; then
  echo "ERROR: climate file missing: $CLIMATE" >&2
  exit 1
fi

Rscript --no-save "$JOBS/34_lfmm_ksweep.R" \
  "$TRAYS" "$SEEDMIX" "$CLIMATE" "$GEN" "$OUT" "$TMPDIR/lfmm_k_work" \
  "$JOBS/run_lfmm_ksweep.R" 5,10,16,20,25,30
