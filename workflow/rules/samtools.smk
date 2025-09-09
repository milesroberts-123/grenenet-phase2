rule samtools_flagstat:
    input:
        final="bwa_results/{ID}.bam",
        bai="bwa_results/{ID}.bam.bai"
    output:
        "samtools_results/{ID}.txt"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        samtools flagstat {input.final} > {output}
        """

rule samtools_markdup:
    input:
        "bwa_results/{ID}.bam"
    output:
        "markdup_results/{ID}.bam"
    conda:
        "../envs/bwa.yaml"
    shell:
        "samtools markdup {input} {output}"
