#rule msprime:
#    group: "simulation"
#    output:
#        "msprime_results/{ID}.trees"
#    conda:
#        "../envs/msprime.yaml"
#    params:
#        mu=
#        r=
#        N=
#        L=
#        n=
#    shell:
#        "msp ancestry {params.n} -N {params.N} -L {params.L} -r {params.r} | msp mutations {params.mu} > {output}"

rule bcftools_random:
    group: "simulation"
    input:
        fasta=config["sim_fasta"],
        vcf=config["sim_vcf"]
    output:
        bed=temp("random_{ID}.bed"),
        vcf=temp("bcftools_random_results/{ID}.vcf")
    conda:
        "../envs/bcftools.yaml"
    params:
        L=lookup(query="ID == '{ID}'", within=gt_params, cols="L")
    shell:
        """
        bedtools random -l {params.L} -n 1 -g {input.fasta}.fai > {output.bed}
        bcftools view -R {output.bed} -o {output.vcf} {input.vcf}
        """

rule bcftools_duplicate:
    group: "simulation"
    input:
        "bcftools_random_results/{ID}.vcf"
    output:
        "bcftools_duplicate_results/{ID}.vcf"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools merge --force-samples {input} {input} {input} {input} {input} {input} {input} -Ou -o {output}
        """

rule slim:
    group: "simulation"
    input:
        vcf = "bcftools_duplicate_results/{ID}.vcf",
        fasta = config["sim_fasta"]
    output:
        "slim_results/{ID}_neutfreqs.tsv"
    conda:
        "../envs/msprime.yaml"
    params:
        nmu=lookup(query="ID == '{ID}'", within=gt_params, cols="nmu"),
        tmu=lookup(query="ID == '{ID}'", within=gt_params, cols="tmu"),
        R=lookup(query="ID == '{ID}'", within=gt_params, cols="R"),
        N=lookup(query="ID == '{ID}'", within=gt_params, cols="N"),
        L=lookup(query="ID == '{ID}'", within=gt_params, cols="L"),
        sigma=lookup(query="ID == '{ID}'", within=gt_params, cols="sigma"),
        alpha=lookup(query="ID == '{ID}'", within=gt_params, cols="alpha"),
        gamma=lookup(query="ID == '{ID}'", within=gt_params, cols="gamma"),
        tau=lookup(query="ID == '{ID}'", within=gt_params, cols="tau"),
        type=lookup(query="ID == '{ID}'", within=gt_params, cols="type")
    shell:
        """
        if [[ "unstruct" == "{params.type}" ]]; then
            slim -d N={params.N} -d L={params.L} -d nmu={params.nmu} -d tmu={params.tmu} -d R={params.R} -d sigma={params.sigma} -d alpha={params.alpha} -d ID={wildcards.ID} -d gamma={params.gamma} -d tau={params.tau} scripts/gt_expectations.slim
        elif [[ "struct" == "{params.type}" ]]; then
            slim -d fastaFile='"{input.fasta}"' -d vcfFile='"{input.vcf}"' -d N={params.N} -d sigma={params.sigma} -d tmu={params.tmu} -d R={params.R} -d alpha={params.alpha} -d ID={wildcards.ID} scripts/gt_expectations_structured.slim
        elif [[ "bank" == "{params.type}" ]]; then
            slim -d vcfFile='"{input.vcf}"' -d SIGMA={params.sigma} -d MU={params.nmu} -d R={params.R} -d ID={wildcards.ID} scripts/gt_expectations_seed_bank.slim 
        else
            echo "Invalid simulation type!"
        fi
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
