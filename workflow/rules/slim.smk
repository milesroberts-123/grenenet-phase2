rule bcftools_random:
    group: "simulation"
    input:
        fasta=config["sim_fasta"],
        vcf=config["sim_vcf"]
    output:
        bed=temp("random_{ID}.bed"),
        vcf=temp("results/bcftools_random/{ID}.vcf")
    conda:
        "../envs/bcftools.yaml"
    params:
        # G is the genome size: the simulated genome of the wf runs, and the
        # length of the chromosome window the VCF-seeded runs read
        window=lambda wildcards: bank_flag(wildcards, "G")
    shell:
        """
        bedtools random -l {params.window} -n 1 -g {input.fasta}.fai > {output.bed}
        bcftools view -R {output.bed} -o {output.vcf} {input.vcf}
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
    # seed-bank runs always seed from the VCF, whatever the table says: the
    # fresh-genome burn-in path is no longer supported
    return "bank_vcf"


def bank_flag(wildcards, column, default=0):
    """Integer parameter for the bank script; WF rows have no value in these columns."""
    value = gt_params.loc[gt_params["ID"] == wildcards.ID, column].iloc[0]
    return default if pd.isna(value) else int(value)


def float_param(wildcards, column, default=0.0):
    """Float parameter. SLiM rejects an integer literal where a float-typed
    argument is expected (so VA must always render with a decimal point) and
    rejects "nan" outright (undefined identifier), so missing values must be
    replaced by a real number before they reach the command line."""
    value = gt_params.loc[gt_params["ID"] == wildcards.ID, column].iloc[0]
    return default if pd.isna(value) else float(value)


def get_fasta_path(wildcards):
    # a run that simulates its own genome of size G needs no reference FASTA
    if get_simcase(wildcards) in ("wf_vcf", "bank_vcf"):
        return config["sim_fasta"]
    return "none"


def get_vcf_path(wildcards):
    # seed-bank rows are forced to USE_VCF=1, so the VCF is always read there;
    # a run that simulates its own genome (wf_novcf) never reaches this function
    return "results/bcftools_random/%s.vcf" % wildcards.ID


def use_vcf_for(wildcards):
    """The USE_VCF value the SLiM script sees: seed-bank runs are forced to 1."""
    family = str(gt_params.loc[gt_params["ID"] == wildcards.ID, "type"].iloc[0])
    if family == "wf":
        return bank_flag(wildcards, "USE_VCF")
    return 1


rule slim:
    group: "simulation"
    input:
        vcf = branch(
            get_simcase,
            cases={
                "wf_novcf": ("wf_prep/{ID}.done"),
                "wf_vcf": [],   # not wired up yet; see the shell branch
                "bank_vcf": ("results/bcftools_random/{ID}.vcf"),
            }
            ),
        fasta = branch(
            get_simcase,
            cases={
                "wf_novcf": [],
                "wf_vcf": (config["sim_fasta"]),
                "bank_vcf": (config["sim_fasta"]),
            }
            ),
    output:
        "results/slim/{ID}_neutfreqs.tsv"
    conda:
        "../envs/msprime.yaml"
    params:
        # core parameters for wf sims
        MU=lookup(query="ID == '{ID}'", within=gt_params, cols="MU"),
        R=lookup(query="ID == '{ID}'", within=gt_params, cols="R"),
        N=lookup(query="ID == '{ID}'", within=gt_params, cols="N"),
        L=lookup(query="ID == '{ID}'", within=gt_params, cols="L"),
        G=lookup(query="ID == '{ID}'", within=gt_params, cols="G"),
        SIGMA=lookup(query="ID == '{ID}'", within=gt_params, cols="SIGMA"),
        TAU=lookup(query="ID == '{ID}'", within=gt_params, cols="TAU"),
        # extra parameters for non-wf sims. VA and L are shared with the wf grid:
        # L is the number of selected loci and VA the additive variance they carry
        va=lambda wildcards: float_param(wildcards, "VA"),
        germ_rate=lookup(query="ID == '{ID}'", within=gt_params, cols="GERM_RATE"),
        bank_surv=lookup(query="ID == '{ID}'", within=gt_params, cols="BANK_SURV"),
        K=lookup(query="ID == '{ID}'", within=gt_params, cols="K"),
        n_offspring=lookup(query="ID == '{ID}'", within=gt_params, cols="N_OFFSPRING"),
        type=lookup(query="ID == '{ID}'", within=gt_params, cols="type"),
        # SURVIVAL_SELECTION is NA on VA = 0 rows, where it has no effect
        surv=lambda wildcards: float_param(wildcards, "SURVIVAL_SELECTION", 0.5),
        min_age=lookup(query="ID == '{ID}'", within=gt_params, cols="MIN_AGE"),
        max_age=lookup(query="ID == '{ID}'", within=gt_params, cols="MAX_AGE"),
        # USE_VCF = 1 seeds the run from standing variation read out of the VCF and
        # runs no burn-in; USE_VCF = 0 starts from the ancestral FASTA sequence and
        # burns in for BURN generations with the selected mutation types neutral.
        # Seed-bank rows are always forced to 1; only the wf family can be 0.
        use_vcf=use_vcf_for,
        burn=lambda wildcards: bank_flag(wildcards, "BURN"),
        vcf_path=get_vcf_path,
        n_founders=config["sim_founders"],
        fasta_path=get_fasta_path,
        # L is NA for VA = 0 rows (no selected loci); SLiM needs a number
        l_loci=lambda wildcards: bank_flag(wildcards, "L")
    shell:
        """
        if [[ "wf" == "{params.type}" && "0" == "{params.use_vcf}" ]]; then
            slim -d N={params.N} \\
                -d L={params.l_loci} \\
                -d G={params.G} \\
                -d MU={params.MU} \\
                -d R={params.R} \\
                -d SIGMA={params.SIGMA} \\
                -d VA={params.va} \\
                -d ID={wildcards.ID} \\
                -d TAU={params.TAU} scripts/gt_expectations.slim
        elif [[ "wf" == "{params.type}" && "1" == "{params.use_vcf}" ]]; then
            echo "wf with USE_VCF=1 is not wired up: gt_expectations_structured.slim" >&2
            echo "still takes -d tmu/-d alpha, and it hardcodes 231 founders and 10 generations." >&2
            echo "Port it to the VA/L architecture before enabling this branch." >&2
            exit 1
        elif [[ "wf" != "{params.type}" ]]; then
            slim -d FASTA_FILE='"{params.fasta_path}"' -d VCF_FILE='"{params.vcf_path}"' -d N_VCF=5082 -d N_FOUNDERS={params.n_founders} -d SIGMA={params.SIGMA} -d NMU={params.MU} -d G={params.G} -d R={params.R} -d SIM_LENGTH={params.TAU} -d GERM_RATE={params.germ_rate} -d BANK_SURV={params.bank_surv} -d K={params.K} -d MIN_AGE={params.min_age} -d MAX_BANK_AGE={params.max_age} -d N_OFFSPRING={params.n_offspring} -d VA={params.va} -d L={params.l_loci} -d ID={wildcards.ID} -d SURVIVAL_SELECTION={params.surv} -d USE_VCF={params.use_vcf} -d BURN={params.burn} scripts/gt_expectations_seed_bank.slim
        else
            echo "Invalid simulation type!"
        fi
        """
