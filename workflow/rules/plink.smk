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
        "results/plink/greneNet.vcor"
    params:
        bfile = "results/plink/greneNet",
        # plink2 counts --ld-window as "max variant count + 1", so the configured
        # spacing of 500 gets 501 to stay equivalent to the plink 1.9 setting
        plink_ld_window = config["plink_ld_window"] + 1,
        plink_ld_window_kb = config["plink_ld_window_kb"],
        plink_ld_window_r2 = config["plink_ld_window_r2"],
        out = "results/plink/greneNet"
    threads: 32
    conda: "../envs/plink2.yaml"
    # --r-unphased rather than --r2: the signed correlation is what a conditional
    # (multi-locus) estimate needs. The sign is positive when the two MAJOR alleles
    # are positively correlated, and CHROM/POS/ID/MAJ columns are emitted for both
    # variants in each pair so the sign can be realigned to the allele the model
    # carries. --ld-window-r2 filters |r| >= sqrt(threshold), so 0 keeps every pair;
    # --ld-window-kb must be raised explicitly because plink2 defaults to 1000 kb.
    shell:
        r"""
        plink2 --bfile {params.bfile} \
               --r-unphased \
               --ld-window {params.plink_ld_window} \
               --ld-window-kb {params.plink_ld_window_kb} \
               --ld-window-r2 {params.plink_ld_window_r2} \
               --threads {threads} \
               --out {params.out}
        """


rule plink_all:
    input:
        "results/plink/greneNet.vcor"
