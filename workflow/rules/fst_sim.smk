rule grep:
    input:
        expand("msprime_results/{ID}.txt", ID = fst_params["ID"])
    output:
        "fst_results.txt"
    shell:
        """
        grep "" {input} | sed 's/:/\t/g' > {output}
        """

rule msprime:
    group: "fst"
    input:
        "slim_fst_results/{ID}.trees"
    output:
        temp("msprime_results/{ID}.txt")
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
