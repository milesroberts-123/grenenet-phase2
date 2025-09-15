rule kmc_contam_db:
    input:
        read1 = expand("{read1}", read1 = parameters.loc[parameters["experiment"] == "seedmix", "read1"]),
        read2 = expand("{read2}", read2 = parameters.loc[parameters["experiment"] == "seedmix", "read2"]),
    output:
        list="contam.list",
        pre="contam.kmc_pre",
        suf="contam.kmc_suf"
    conda:
        "../envs/kmc.yaml"
    params:
        k = config["k"]
    shell:
        """
        if [ ! -d "kmc_tmp_dir" ]; then
            mkdir kmc_tmp_dir/
        fi

        # create file list
        echo {input.read1} {input.read2} | tr ' ' '\n' > {output.list}

        # build kmc db
        kmc -k{params.k} -m36 -t{threads} -ci1 -cs2 -fq @{output.list} contam kmc_tmp_dir/

        # clean up
        rm -r kmc_tmp_dir/
        """
