rule bwa_mem:
    input:
        ref=config["reference"],
        amb=config["reference"] + ".amb",
        ann=config["reference"] + ".ann",
        bwt=config["reference"] + ".bwt",
        pac=config["reference"] + ".pac",
        sa=config["reference"] + ".sa", 
        reads="no_contam_reads/{ID}.fastq",
    output:
        final="bwa_results/{ID}.bam",
        bai="bwa_results/{ID}.bam.bai"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        # align reads to reference
        echo Aligning reads...
        bwa mem -R '@RG\\tID:{wildcards.ID}\\tSM:{wildcards.ID}' -t {threads} {input.ref} {input.reads} | samtools sort -O bam > {output.final}
        
        # index final output
        echo Indexing final output...
        samtools index {output.final}
        """
