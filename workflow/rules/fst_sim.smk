rule grep:
    input:
        expand("msprime_results/{ID}.txt", ID = fst_params["ID"])
    output:
        "fst_results.txt"
    shell:
        """
        grep "" {input} | sed 's/:/\t/g' > {output}
        """

rule gcta:
    group: "fst"
    input:
        bed="plink_results/{ID}.bed",
        popfile="popfiles/{ID}.txt",
        fam="plink_results/{ID}.fam"
    output:
        "gcta_results/{ID}.fst",
        temp("gcta_results/{ID}.log")
    conda: 
        "../envs/gcta.yaml"
    shell:
        "gcta64 --bfile plink_results/{wildcards.ID} --fst --sub-popu {input.popfile} --out gcta_results/{wildcards.ID}"

rule popfile:
    group: "fst"
    input:
        "plink_results/{ID}.fam"
    output:
        temp("popfiles/{ID}.txt")
    shell:
        r"""
        sed 's:\(mod_.*\):\1\tmodern:g' {input} | sed 's:\(his_.*\):\1\thistorical:g' | cut -f 1,2,7 > {output}
        """

rule plink:
    group: "fst"
    input:
        "bcftools_merge_results/merged_{ID}.vcf"
    output:
        temp("plink_results/{ID}.bed"),
        temp("plink_results/{ID}.fam")
    conda:
        "../envs/plink2.yaml"
    shell:
        "plink2 --vcf {input} --make-bed --out plink_results/{wildcards.ID}"

rule bcftools:
    group: "fst"
    input:
        "msprime_results/modern_{ID}.vcf",
        "msprime_results/historical_{ID}.vcf"
    output:
        temp("bcftools_merge_results/merged_{ID}.vcf")
    conda: 
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools merge -0 -O z -o {output} --no-index {input}
        """

rule msprime:
    group: "fst"
    input:
        "slim_fst_results/{ID}.trees"
    output:
        temp("msprime_results/{ID}.txt"),
        temp("msprime_results/modern_{ID}.vcf"),
        temp("msprime_results/historical_{ID}.vcf")
    conda:
        "../envs/msprime.yaml"
    params:
        tau=lookup(query="ID == '{ID}'", within=fst_params, cols="tau"),
        mu=lookup(query="ID == '{ID}'", within=fst_params, cols="mu"),
        R=lookup(query="ID == '{ID}'", within=fst_params, cols="R"),
        N=lookup(query="ID == '{ID}'", within=fst_params, cols="msprimeN"),
        L=lookup(query="ID == '{ID}'", within=fst_params, cols="L"),
        sigma=lookup(query="ID == '{ID}'", within=fst_params, cols="sigma"),
    shell:
        """
        python scripts/burnin.py --mu {params.mu} --tau {params.tau} -N {params.N} -L {params.L} -R {params.R} --ID {wildcards.ID}
        """

rule fst_slim:
    group: "fst"
    output:
        temp("slim_fst_results/{ID}.trees")
    conda:
        "../envs/msprime.yaml"
    params:
        mu=lookup(query="ID == '{ID}'", within=fst_params, cols="mu"),
        R=lookup(query="ID == '{ID}'", within=fst_params, cols="R"),
        N=lookup(query="ID == '{ID}'", within=fst_params, cols="slimN"),
        L=lookup(query="ID == '{ID}'", within=fst_params, cols="L"),
        sigma=lookup(query="ID == '{ID}'", within=fst_params, cols="sigma"),
    shell:
        "slim -d N={params.N} -d L={params.L} -d mu={params.mu} -d R={params.R} -d sigma={params.sigma} -d ID={wildcards.ID} scripts/fst_expectations.slim"
