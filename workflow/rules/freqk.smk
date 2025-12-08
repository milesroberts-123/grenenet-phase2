rule cat:
    input:
        expand(config["seed_root"] + "{{sample}}{mate}", mate=config["seed_suffix"]),
    output:
        temp("cat_results/{sample}.fq"),
    shell:
        "zcat {input} > {output}"

rule freqk_count:
    input:
        index=config["freqk_index"],
        read="cat_results/{sample}.fq"
    output:
        freqs="freqk_count_results/counts_by_allele_{sample}.txt",
        kmers="freqk_count_results/counts_by_kmer_{sample}.txt"
    shell:
        "scripts/freqk count -n {threads} -i {input.index} -r {input.read} -f {output.freqs} -c {output.kmers}"

rule freqk_call:
    input:
        freqs="freqk_count_results/counts_by_allele_{sample}.txt",
        index=config["freqk_index"]
    output:
        "freqk_call_results/{sample}.txt"
    shell:
        "scripts/freqk call -i {input.index} -c {input.freqs} -o {output}"
