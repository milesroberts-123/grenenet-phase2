wildcard_constraints:
    sample="[^/]+"

rule ember_seurat_counts:
    input:
        rds=lambda wildcards: [
            p for p in config["seurat_rds"]
            if os.path.splitext(os.path.basename(p))[0] == wildcards.sample
        ][0]
    output:
        directory("results/ember/cells/{sample}")
    log:
        "logs/ember/counts/{sample}.log"
    conda:
        "../envs/seurat.yaml"
    shell:
        """
        Rscript scripts/seurat_to_counts.R {input.rds} {output} 2> {log}
        """

rule ember_entropy:
    input:
        cells=expand("results/ember/cells/{sample}", sample=seurat_samples)
    output:
        metrics="results/ember/entropy_metrics_stage_celltype.csv",
        psi_block="results/ember/psi_block_stage_celltype.csv"
    params:
        min_cells=config.get("ember_min_cells", 100)
    log:
        "logs/ember/entropy.log"
    conda:
        "../envs/seurat.yaml"
    shell:
        """
        Rscript scripts/ember_entropy.R {input.cells} {params.min_cells} {output.metrics} {output.psi_block} 2> {log}
        """