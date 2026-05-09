rule kmc_filter:
    input:
        pre="after_{sample}.kmc_pre",
        suf="after_{sample}.kmc_suf",
        pread1="fastp_results/trimmed_paired_R1_{sample}.fastq.gz",
        pread2="fastp_results/trimmed_paired_R2_{sample}.fastq.gz",
        uread1="fastp_results/trimmed_unpaired_R1_{sample}.fastq.gz",
        uread2="fastp_results/trimmed_unpaired_R2_{sample}.fastq.gz",
    output:
        reads = "filter_results/{sample}.fastq.gz",
        filter_list = "filter_list_{sample}.txt"
    conda:
        "../envs/kmc.yaml"
    params:
        contam_match_limit_percent = config["contam_match_limit_percent"]
    shell:
        """
        # create file list
        echo {input.pread1} {input.pread2} {input.uread1} {input.uread2} | tr ' ' '\n' > {output.filter_list}

        kmc_tools filter after_{wildcards.sample} @{output.filter_list} -ci{params.contam_match_limit_percent} -cx1.0 {output.reads}
        """
