rule bwa_mem:
    input:
        ref=config["reference"],
        amb=config["reference"] + ".amb",
        ann=config["reference"] + ".ann",
        bwt=config["reference"] + ".bwt",
        pac=config["reference"] + ".pac",
        sa=config["reference"] + ".sa", 
        reads="filter_results/{sample}.fastq.gz",
    output:
        final="bwa_results/{sample}.bam",
        bai="bwa_results/{sample}.bam.bai"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        # align reads to reference
        echo Aligning reads...
        bwa mem -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}' -t {threads} {input.ref} {input.reads} | samtools sort -O bam > {output.final}
        
        # index final output
        echo Indexing final output...
        samtools index {output.final}
        """
