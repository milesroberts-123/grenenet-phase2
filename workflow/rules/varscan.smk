rule varscan:
    input:
        reffasta=config["reference"],
        trimbam="markdup_results/{ID}.bam",
    output:
        "varscan_results/{ID}.tsv",
    conda:
        "../envs/varscan.yaml"
    shell:
        """
        samtools mpileup -f {input.reffasta} {input.trimbam} | varscan pileup2snp --min-coverage 10 --min-avg-qual 40 --min-reads2 10 1> {output}
        """
