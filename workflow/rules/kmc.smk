rule kmc:
    input:
        pread1="fastp_results/trimmed_paired_R1_{ID}.fastq.gz",
        pread2="fastp_results/trimmed_paired_R2_{ID}.fastq.gz",
        uread1="fastp_results/trimmed_unpaired_R1_{ID}.fastq.gz",
        uread2="fastp_results/trimmed_unpaired_R2_{ID}.fastq.gz",
    output:
        list=temp("temp_input_{ID}.txt"),
        tmp_pre=temp("before_{ID}.kmc_pre"),
        tmp_suf=temp("before_{ID}.kmc_suf")
    conda:
        "../envs/kmc.yaml"
    params:
        mincount=config["mincount"],
        maxcount=config["maxcount"],
        k=config["k"],
    shell:
        """
        # create directory
        if [ -d "tmp_kmc_{wildcards.ID}" ]; then
            rm -r tmp_kmc_{wildcards.ID}
        fi

        mkdir tmp_kmc_{wildcards.ID}

        echo {input.pread1} {input.pread2} {input.uread1} {input.uread2} | tr ' ' '\n' > {output.list}

        # count k-mers
        kmc -m15 -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} @{output.list} before_{wildcards.ID} tmp_kmc_{wildcards.ID}

        # delete tmp directories
        rm -r tmp_kmc_{wildcards.ID}
        """
