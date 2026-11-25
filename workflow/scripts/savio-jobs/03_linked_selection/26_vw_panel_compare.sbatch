#!/bin/bash
#SBATCH --job-name=26_vw_panel_compare
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
# One task: builds the founder haplotype matrix for a single chromosome, then
# compares the s of a nested 10k/50k locus window. Needs the asreml seat because
# library(Vw) attaches asreml, so this must not run alongside another Vw task.
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=milesroberts@berkeley.edu
#SBATCH --output=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/workflow/scripts/savio-jobs/03_linked_selection/logs/%x_%A.out
#SBATCH --chdir=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/resources

set -euo pipefail

# R cannot create a temp dir on Savio's login/compute /tmp
export TMPDIR=/global/scratch/users/milesroberts/tmp
mkdir -p "$TMPDIR"
module load r

ROOT=/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2
VCF=$ROOT/resources/greneNet_final_v1.1.recode.vcf.gz
BCF=$ROOT/workflow/.snakemake/conda/7527142a2da7f7dea129024c4e749abc_/bin/bcftools
OUT=$ROOT/results/$(date +%F)
mkdir -p "$OUT"

# one chromosome, drawn reproducibly so the run can be repeated
CHROM=$(Rscript --no-save -e 'set.seed(20261124); cat(sample(c("1","2","3","4","5"), 1))')
echo "chromosome: $CHROM"

# haplotypes x loci, 0/1, MAF >= 0.05, complete calls only.
# The VCF is phased, so each sample contributes its two haplotypes directly.
CACHE=$TMPDIR/vw_founder_chr${CHROM}_v2.tsv
if [[ ! -s "$CACHE" ]]; then
  echo "building $CACHE"
  "$BCF" query -r "$CHROM" -f '%POS\t[%GT\t]\n' "$VCF" 2>/dev/null \
  | awk -F'\t' 'BEGIN{OFS="\t"}
      {
        # build the line from pos with OFS so no empty field can appear, and skip
        # the empty field left by the trailing tab of [%GT\t]
        pos=$1; s=0; nh=0; line=pos; bad=0
        for (i=2; i<=NF; i++) {
          gt=$i
          if (gt=="") continue
          if (gt=="." || gt=="./." || gt==".|.") { bad=1; break }
          a=substr(gt,1,1); b=substr(gt,3,1)
          if (a=="." || b==".") { bad=1; break }
          s += a + b; nh += 2
          line = line OFS a OFS b
        }
        if (bad || nh == 0) next
        af = s/nh
        if (af < 0.05 || af > 0.95) next
        print line
      }' > "$CACHE"
fi
NCOL=$(head -1 "$CACHE" | awk -F'\t' '{print NF-1}')
echo "cache: $(wc -l < "$CACHE") loci x ${NCOL} haplotypes"
if [[ "$NCOL" -ne 462 ]]; then
  echo "ERROR: expected 462 haplotypes (231 accessions x 2), got ${NCOL}" >&2
  exit 1
fi

Rscript --no-save "$ROOT/workflow/scripts/savio-jobs/03_linked_selection/26_vw_panel_compare.R" \
  "$CACHE" "$OUT" "$CHROM"
