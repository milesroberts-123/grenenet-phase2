rule multiqc:
    input:
        expand("fastp_results/{ID}_R1R2.json", ID = parameters["sample"]),
    output:
        "multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc fastp_results/"
