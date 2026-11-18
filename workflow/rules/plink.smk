rule plink_convert:
    input:
        config["founder_vcf"]
    output:
        fam = "results/plink/greneNet.fam",
        bed = "results/plink/greneNet.bed",
        bim = "results/plink/greneNet.bim"
    params:
        prefix = "results/plink/greneNet"
    conda: "../envs/plink.yaml"
    shell:
        "plink --vcf {input} --biallelic-only --make-bed --out {params.prefix}"


rule plink_ld:
    input:
        fam = "results/plink/greneNet.fam",
        bed = "results/plink/greneNet.bed",
        bim = "results/plink/greneNet.bim"
    output:
        temp("results/plink/greneNet.ld")
    params:
        plink_ld_window = config["plink_ld_window"],
        plink_ld_window_kb = config["plink_ld_window_kb"],
        plink_ld_window_r2 = config["plink_ld_window_r2"]
    threads: 32
    conda: "../envs/plink.yaml"
    shell:
        r"""
        plink --bed {input.bed} --bim {input.bim} --fam {input.fam} \
              --r2 d \
              --ld-window {params.plink_ld_window} \
              --ld-window-kb {params.plink_ld_window_kb} \
              --ld-window-r2 {params.plink_ld_window_r2} \
              --threads {threads} \
              --out results/plink/greneNet
        """


rule plink_tsv:
    input:
        "results/plink/greneNet.ld"
    output:
        "results/plink/greneNet.ld.tsv"
    # plink pads the columns with spaces; keep SNP_A, SNP_B and R2 (fields 3, 6, 7
    # of CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2)
    shell:
        r"tr -s ' ' '\t' < {input} | cut -f3,6,7 > {output}"


rule plink_all:
    input:
        "results/plink/greneNet.ld.tsv"
