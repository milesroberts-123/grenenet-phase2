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

rule bgzip:
    group: "simulation"
    input:
        "bcftools_random_results/{ID}.vcf"
    output:
        temp("bcftools_random_results/{ID}.vcf.gz"),
        temp("bcftools_random_results/{ID}.vcf.gz.tbi")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bgzip {input}
        tabix {input}.gz
        """

rule bcftools_duplicate:
    group: "simulation"
    input:
        gz="bcftools_random_results/{ID}.vcf.gz",
        tbi="bcftools_random_results/{ID}.vcf.gz.tbi"
    output:
        temp("bcftools_duplicate_results/{ID}.vcf")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools merge --force-samples {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} {input.gz} -Ou -o {output}
        """

rule unstruct_prep:
    group: "simulation"
    output:
        temp(touch("unstruct_prep/{ID}.done"))

def get_simtype(wildcards):
    simtype = gt_params.loc[gt_params["ID"] == wildcards.ID, "type"]
    simtype = str(simtype.iloc[0])
    return(simtype)

rule slim:
    group: "simulation"
    input:
        #vcf = "bcftools_duplicate_results/{ID}.vcf",
        vcf = branch(
            get_simtype,
            cases={
                "unstruct": ("unstruct_prep/{ID}.done"),
                "struct": ("bcftools_random_results/{ID}.vcf"),
                "bank": ("bcftools_duplicate_results/{ID}.vcf")
            }
            ),
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
        germ_rate=lookup(query="ID == '{ID}'", within=gt_params, cols="GERM_RATE"),
        bank_surv=lookup(query="ID == '{ID}'", within=gt_params, cols="BANK_SURV"),
        K=lookup(query="ID == '{ID}'", within=gt_params, cols="K"),
        n_offspring=lookup(query="ID == '{ID}'", within=gt_params, cols="N_OFFSPRING"),
        type=lookup(query="ID == '{ID}'", within=gt_params, cols="type"),
        surv=lookup(query="ID == '{ID}'", within=gt_params, cols="SURVIVAL_SELECTION"),
        min_age=lookup(query="ID == '{ID}'", within=gt_params, cols="MIN_AGE"),
        max_age=lookup(query="ID == '{ID}'", within=gt_params, cols="MAX_AGE")
    shell:
        """
        if [[ "unstruct" == "{params.type}" ]]; then
            slim -d N={params.N} -d L={params.L} -d nmu={params.nmu} -d tmu={params.tmu} -d R={params.R} -d sigma={params.sigma} -d alpha={params.alpha} -d ID={wildcards.ID} -d gamma={params.gamma} -d tau={params.tau} scripts/gt_expectations.slim
        elif [[ "struct" == "{params.type}" ]]; then
            slim -d fastaFile='"{input.fasta}"' -d vcfFile='"{input.vcf}"' -d N={params.N} -d sigma={params.sigma} -d tmu={params.tmu} -d R={params.R} -d alpha={params.alpha} -d ID={wildcards.ID} scripts/gt_expectations_structured.slim
        elif [[ "bank" == "{params.type}" ]]; then
            slim -d FASTA_FILE='"{input.fasta}"' -d VCF_FILE='"{input.vcf}"' -d N_VCF=5082 -d SIGMA={params.sigma} -d NMU={params.nmu} -d TMU={params.tmu} -d R={params.R} -d SIM_LENGTH={params.tau} -d GERM_RATE={params.germ_rate} -d BANK_SURV={params.bank_surv} -d K={params.K} -d MIN_AGE={params.min_age} -d MAX_BANK_AGE={params.max_age} -d N_OFFSPRING={params.n_offspring} -d ALPHA={params.alpha} -d ID={wildcards.ID} -d SURVIVAL_SELECTION={params.surv} scripts/gt_expectations_seed_bank.slim 
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
