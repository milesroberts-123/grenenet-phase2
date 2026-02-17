print("Loading packages...")
import fastdfe as fd

print("Annotating ancestral alleles...")
ann = fd.Annotator(
    vcf="greneNet_final_v1.1.recode.vcf.gz",
    fasta="Athaliana_447_TAIR10.fa",
    #gff="http://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/"
    #    "Homo_sapiens.GRCh38.109.chromosome.21.gff3.gz",
    output='greneNet.ancestral.vcf.gz',
    annotations=[fd.MaximumParsimonyAncestralAnnotation()]
)

ann.annotate()
