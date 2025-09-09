rule bwa_index:
    input:
        config["reference"],
    output:
        amb=temp(config["reference"] + ".amb"),
        ann=temp(config["reference"] + ".ann"),
        bwt=temp(config["reference"] + ".bwt"),
        pac=temp(config["reference"] + ".pac"),
        sa=temp(config["reference"] + ".sa"),
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        # index reference
        bwa index {input}
        """
