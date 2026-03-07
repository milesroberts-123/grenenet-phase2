rule msprime:
    group: "fst"
    output:
        temp("msprime_results/{ID}.trees")
    conda:
        "../envs/msprime.yaml"
    params:
        mu=lookup(query="ID == '{ID}'", within=fst_params, cols="mu"),
        R=lookup(query="ID == '{ID}'", within=fst_params, cols="R"),
        N=lookup(query="ID == '{ID}'", within=fst_params, cols="N"),
        L=lookup(query="ID == '{ID}'", within=fst_params, cols="L"),
        sigma=lookup(query="ID == '{ID}'", within=fst_params, cols="sigma"),
    shell:
        """
        python scripts/burnin.py -N {params.N} -L {params.L} -R {params.R} --ID {wildcards.ID}
        """

rule slim:
    group: "fst"
    input:
        "msprime_results/{ID}.trees"
    output:
        "slim_fst_results/{ID}_fst.tsv"
    conda:
        "../envs/msprime.yaml"
    params:
        mu=lookup(query="ID == '{ID}'", within=fst_params, cols="mu"),
        R=lookup(query="ID == '{ID}'", within=fst_params, cols="R"),
        N=lookup(query="ID == '{ID}'", within=fst_params, cols="N"),
        L=lookup(query="ID == '{ID}'", within=fst_params, cols="L"),
        sigma=lookup(query="ID == '{ID}'", within=fst_params, cols="sigma"),
    shell:
        "slim -d N={params.N} -d L={params.L} -d mu={params.mu} -d R={params.R} -d sigma={params.sigma} -d ID={wildcards.ID} scripts/fst_expectations.slim"

rule bgzip:
    group: "fst"
    input:
        "slim_results/{ID}_{population}_{time}.vcf"
    output:
        "slim_results/{ID}_{population}_{time}.vcf.gz",
        "slim_results/{ID}_{population}_{time}.vcf.gz.tbi"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bgzip {input}
        tabix {input}.gz
        """

rule bcftools_merge:
    group: "fst"
    input:
        expand("slim_results/{{ID}}_{population}_{time}.vcf.gz", time = fst_params["time"])
    output:
        "bcftools_results/{ID}_{population}.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools merge --force-samples --missing-to-ref -Oz -o {output} {input}"

rule bcftools_query:
    group: "fst"
    input: "bcftools_results/{ID}_{population}.vcf.gz"
    output: "allele_frequencies/{ID}_{population}.tsv"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # output allele depths
        bcftools view -m2 -M2 -v snps {input} | bcftools query -f '%CHROM %POS %REF %ALT %NS %AF %AC\n' > {output}
        """

#rule tskit:
#    group: "simulation"
#    input:
#        "slim_results/{ID}.trees"
#    output:
#        "tskit_results/{ID}.vcf"
#    conda:
#        "../envs/tskit.yaml"
#    shell:
#        "tskit vcf {input} > {output}"
