wildcard_constraints:
    sample="[^/]+"

rule genes_bed:
    group: "expression_bins"
    input:
        gff=config["gff"]
    output:
        "expression_bins/ref/genes_coords.bed"
    shell:
        """
        grep "gene" {input.gff} | cut -f 1,4,5,9 | sed 's:ID.*Name=::g' > {output}
        """

rule expression_bed:
    group: "expression_bins"
    input:
        rds=lambda wildcards: [
            p for p in config["seurat_rds"]
            if os.path.splitext(os.path.basename(p))[0] == wildcards.sample
        ][0],
        genes="expression_bins/ref/genes_coords.bed"
    output:
        "expression_bins/{sample}.bed"
    conda:
        "../envs/seurat.yaml"
    params:
        group_by=config["expression_bins_group_by"]
    shell:
        """
        Rscript scripts/expression_bins.R {input.rds} {input.genes} {output} {params.group_by}
        """

rule expression_bins:
    group: "expression_bins"
    input:
        windows=config["windows_bed"],
        bed="expression_bins/{sample}.bed"
    output:
        "expression_bins/{sample}.bins",
        temp("expression_bins/{sample}.sorted.bed")
    conda:
        "../envs/bcftools.yaml"
    params:
        map_cols=",".join(str(c) for c in config["expression_bins_map_cols"]),
        map_ops=",".join(["sum"] * len(config["expression_bins_map_cols"]))
    shell:
        """
        tail -n +2 {input.bed} | sort -k1,1 -k2,2n > {output[1]}
        bedtools map -a {input.windows} -b {output[1]} -c {params.map_cols} -o {params.map_ops} > {output[0]}
        """

rule plot_expression_bins:
    group: "expression_bins"
    input:
        bins="expression_bins/{sample}.bins",
        bed="expression_bins/{sample}.bed"
    output:
        "expression_bins/{sample}.png"
    conda:
        "../envs/seurat.yaml"
    params:
        width=config["expression_bins_plot_width"]
    shell:
        """
        Rscript scripts/plot_expression_bins.R {input.bins} {input.bed} {output} {params.width}
        """
