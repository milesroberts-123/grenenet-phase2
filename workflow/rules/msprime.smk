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

rule slim:
    group: "simulation"
    #input:
    #    "msprime_results/{ID}.trees"
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
    shell:
        "slim -d N={params.N} -d L={params.L} -d nmu={params.nmu} -d tmu={params.tmu} -d R={params.R} -d sigma={params.sigma} -d alpha={params.alpha} -d ID={wildcards.ID} -d gamma={params.gamma} -d tau={params.tau} scripts/gt_expectations.slim"

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
