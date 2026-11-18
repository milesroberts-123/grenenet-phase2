rule plink:
    input: config["founder_vcf"]
    output:
    conda: "../envs/plink.yaml"
    shell:
        ""
