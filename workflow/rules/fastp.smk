rule fastp:
    input:
        read1=lookup(query="sample == '{sample}'", within=parameters, cols=["read1"]),
        read2=lookup(query="sample == '{sample}'", within=parameters, cols=["read2"]),
    output:
        pread1="fastp_results/trimmed_paired_R1_{sample}.fastq.gz",
        pread2="fastp_results/trimmed_paired_R2_{sample}.fastq.gz",
        uread1="fastp_results/trimmed_unpaired_R1_{sample}.fastq.gz",
        uread2="fastp_results/trimmed_unpaired_R2_{sample}.fastq.gz",
        jsonR1R2="fastp_results/{sample}_R1R2.json",
    conda:
        "../envs/fastp.yaml"
    params:
        unqualLimit=config["unqualLimit"],
        k=config["k"],
        qualThresh=config["qualThresh"],
        windowLength=config["windowLength"],
        nBaseLimit=config["nBaseLimit"]
    shell:
        """
        # remove duplicates, do read correction, drop low quality reads
        fastp --thread {threads} --n_base_limit {params.nBaseLimit} -u {params.unqualLimit} -q {params.qualThresh} -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --correction -i {input.read1} -I {input.read2} -o {output.pread1} -O {output.pread2} --unpaired1 {output.uread1} --unpaired2 {output.uread2} --json {output.jsonR1R2}
        """

rule qc_post_rm_contam:
    input:
        "no_contam_reads/{ID}.fastq"    
    output:
        temp("fastp_results/no_contam_{ID}.json")
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp -A -Q -L -G -i {input} --json {output}"
