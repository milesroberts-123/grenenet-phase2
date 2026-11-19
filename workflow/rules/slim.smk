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

rule wf_prep:
    group: "simulation"
    output:
        temp(touch("wf_prep/{ID}.done"))

def get_simcase(wildcards):
    row = gt_params.loc[gt_params["ID"] == wildcards.ID, ["type", "USE_VCF"]]
    family = str(row["type"].iloc[0])
    use_vcf = 0 if pd.isna(row["USE_VCF"].iloc[0]) else int(row["USE_VCF"].iloc[0])
    if family == "wf":
        return "wf_vcf" if use_vcf == 1 else "wf_novcf"
    return "bank_vcf" if use_vcf == 1 else "bank_novcf"


def bank_flag(wildcards, column, default=0):
    """Integer parameter for the bank script; WF rows have no value in these columns."""
    value = gt_params.loc[gt_params["ID"] == wildcards.ID, column].iloc[0]
    return default if pd.isna(value) else int(value)


def get_vcf_path(wildcards):
    if bank_flag(wildcards, "USE_VCF") == 0:
        return "none"   # never read, but SLiM needs the string constant defined
    family = str(gt_params.loc[gt_params["ID"] == wildcards.ID, "type"].iloc[0])
    prefix = "bcftools_random_results" if family == "wf" else "bcftools_duplicate_results"
    return "%s/%s.vcf" % (prefix, wildcards.ID)


rule slim:
    group: "simulation"
    input:
        #vcf = "bcftools_duplicate_results/{ID}.vcf",
        vcf = branch(
            get_simcase,
            cases={
                "wf_novcf": ("wf_prep/{ID}.done"),
                "wf_vcf": [],   # not wired up yet; see the shell branch
                "bank_vcf": ("bcftools_duplicate_results/{ID}.vcf"),
                # no VCF is built for a run that starts from the ancestral
                # sequence, which keeps 1000s of duplicate VCFs out of the DAG
                "bank_novcf": []
            }
            ),
        fasta = config["sim_fasta"]
    output:
        "slim_results/{ID}_neutfreqs.tsv"
    conda:
        "../envs/msprime.yaml"
    params:
        # core parameters for wf sims
        MU=lookup(query="ID == '{ID}'", within=gt_params, cols="MU"),
        R=lookup(query="ID == '{ID}'", within=gt_params, cols="R"),
        N=lookup(query="ID == '{ID}'", within=gt_params, cols="N"),
        L=lookup(query="ID == '{ID}'", within=gt_params, cols="L"),
        G=lookup(query="ID == '{ID}'", within=gt_params, cols="G"),
        VA=lookup(query="ID == '{ID}'", within=gt_params, cols="VA"),
        SIGMA=lookup(query="ID == '{ID}'", within=gt_params, cols="SIGMA"),
        TAU=lookup(query="ID == '{ID}'", within=gt_params, cols="TAU"),
        # extra parameters for non-wf sims. VA and L are shared with the wf grid:
        # L is the number of selected loci and VA the additive variance they carry
        va=lookup(query="ID == '{ID}'", within=gt_params, cols="VA"),
        #gamma=lookup(query="ID == '{ID}'", within=gt_params, cols="gamma"),
        germ_rate=lookup(query="ID == '{ID}'", within=gt_params, cols="GERM_RATE"),
        bank_surv=lookup(query="ID == '{ID}'", within=gt_params, cols="BANK_SURV"),
        K=lookup(query="ID == '{ID}'", within=gt_params, cols="K"),
        n_offspring=lookup(query="ID == '{ID}'", within=gt_params, cols="N_OFFSPRING"),
        type=lookup(query="ID == '{ID}'", within=gt_params, cols="type"),
        surv=lookup(query="ID == '{ID}'", within=gt_params, cols="SURVIVAL_SELECTION"),
        min_age=lookup(query="ID == '{ID}'", within=gt_params, cols="MIN_AGE"),
        max_age=lookup(query="ID == '{ID}'", within=gt_params, cols="MAX_AGE"),
        # USE_VCF = 1 starts from standing variation read out of the VCF and runs
        # no burn-in; USE_VCF = 0 starts from the ancestral FASTA sequence and
        # burns in for BURN generations with the selected mutation types neutral
        use_vcf=lambda wildcards: bank_flag(wildcards, "USE_VCF"),
        burn=lambda wildcards: bank_flag(wildcards, "BURN"),
        vcf_path=get_vcf_path
    shell:
        """
        if [[ "wf" == "{params.type}" && "0" == "{params.use_vcf}" ]]; then
            slim -d N={params.N} \\
                -d L={params.L} \\
                -d G={params.G} \\
                -d MU={params.MU} \\
                -d R={params.R} \\
                -d SIGMA={params.SIGMA} \\
                -d VA={params.VA} \\
                -d ID={wildcards.ID} \\
                -d TAU={params.TAU} scripts/gt_expectations.slim
        elif [[ "wf" == "{params.type}" && "1" == "{params.use_vcf}" ]]; then
            echo "wf with USE_VCF=1 is not wired up: gt_expectations_structured.slim" >&2
            echo "still takes -d tmu/-d alpha, and it hardcodes 231 founders and 10 generations." >&2
            echo "Port it to the VA/L architecture before enabling this branch." >&2
            exit 1
        elif [[ "bank" == "{params.type}" ]]; then
            slim -d FASTA_FILE='"{input.fasta}"' -d VCF_FILE='"{params.vcf_path}"' -d N_VCF=5082 -d SIGMA={params.SIGMA} -d NMU={params.MU} -d R={params.R} -d SIM_LENGTH={params.TAU} -d GERM_RATE={params.germ_rate} -d BANK_SURV={params.bank_surv} -d K={params.K} -d MIN_AGE={params.min_age} -d MAX_BANK_AGE={params.max_age} -d N_OFFSPRING={params.n_offspring} -d VA={params.va} -d L={params.L} -d ID={wildcards.ID} -d SURVIVAL_SELECTION={params.surv} -d USE_VCF={params.use_vcf} -d BURN={params.burn} scripts/gt_expectations_seed_bank.slim
        else
            echo "Invalid simulation type!"
        fi
        """
