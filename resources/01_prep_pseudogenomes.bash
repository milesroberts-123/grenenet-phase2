# input files
VCF="greneNet_final_v1.1.recode.vcf.gz"
REF="Athaliana_447_TAIR10.fa"

#HAPLOTYPES=("R" "A")

# get array of sample names
bcftools query -l $VCF > sample_names.txt


# read list into bash array
readarray -t sample_array < sample_names.txt

# loop over every haplotype and create pseudogenome sequences
for id in "${sample_array[@]}"; do
for hap in {1..2}; do
  echo "Current sample: $id"
  echo "Current haplotype: $hap"
  bcftools consensus --haplotype $hap -f $REF -s $id $VCF > $(echo $id)_$(echo $hap).fa
done
done
