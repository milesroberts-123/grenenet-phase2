rule msprime:
    group: "simulation"
    output:
        "msprime_results/{ID}.trees"
    conda:
        "../envs/msprime.yaml"
    params:
        mu=
        r=
        N=
        L=
        n=
    shell:
        "msp ancestry {params.n} -N {params.N} -L {params.L} -r {params.r} | msp mutations {params.mu} > {output}"

rule slim:
    group: "simulation"
    input:
        "msprime_results/{ID}.trees"
    output:
        "slim_results/{ID}.trees"
    conda:
        "../envs/slim.yaml"
    params:
        mu=
        r=
        N=
        L=
        n=
    shell:
        "slim"

rule tskit:
    group: "simulation"
    input:
        "slim_results/{ID}.trees"
    output:
        "tskit_results/{ID}.vcf"
    conda:
        "../envs/tskit.yaml"
    shell:
        "tskit vcf {input} > {output}"
