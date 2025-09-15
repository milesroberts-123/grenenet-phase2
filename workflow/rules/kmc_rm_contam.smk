rule kmc_rm_contam:
    input:
        sample_pre = "before_{sample}.kmc_pre",
        sample_suf = "before_{sample}.kmc_suf",
        contam_pre = "contam.kmc_pre",
        contam_suf = "contam.kmc_suf"
    output:
        before="kmc_results/before_{sample}.txt",
        after="kmc_results/after_{sample}.txt",
        pre=temp("after_{sample}.kmc_pre"),
        suf=temp("after_{sample}.kmc_suf")
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        # filter reads for contamination
        kmc_tools -t{threads} simple before_{wildcards.sample} contam kmers_subtract after_{wildcards.sample}
        
        # count number of k-mers before vs after
        kmc_tools transform before_{wildcards.sample} dump {output.before}
        kmc_tools transform after_{wildcards.sample} dump {output.after}
        """
