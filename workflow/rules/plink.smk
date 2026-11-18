rule plink_convert:
    intput:
        config["founder_vcf"]
    output:
        fam = "results/plink/greneNet.fam",
        bed = "results/plink/greneNet.bed",
        bim = "results/plink/greneNet.bim"
    params:
        prefix="results/plink/greneNet"
    conda: "../envs/plink.yaml"
    shell:
        "plink --vcf {input} --make-bed --out {params.prefix}"

rule plink_ld:
    input:
        fam = "results/plink/greneNet.fam",
        bed = "results/plink/greneNet.bed",
        bim = "results/plink/greneNet.bim"
    output:
        temp("plink.ld")
    params:
        plink_ld_window=config["plink_ld_window"],
        plink_ld_window_kb=config["plink_ld_window_kb"],
        plink_ld_window_r2=config["plink_ld_window_r2"]
    conda: "../envs/plink.yaml"
    shell:
        "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --r2 --ld-window {params.plink_ld_window} --ld-window-kb {params.plink_ld_window_kb} --ld-window-r2 {params.plink_ld_window_r2}

rule plink_tsv:
    input: "plink.ld"
    output: "plink.ld.tsv"
    shell: "sed 's/ \+/\t/g' {input} | cut -f 1 --complement > {output}"


rule plink_all:
    input: "plink.ld.tsv"
