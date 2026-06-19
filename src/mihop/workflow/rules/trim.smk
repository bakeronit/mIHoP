rule trim_illumina_pe:
    input:
        get_fastqs,
    output:
        trimmed = [
            "results/fastp/{sample}/{sample}_R1.fastq.gz",
            "results/fastp/{sample}/{sample}_R2.fastq.gz",
        ],
        html      = report("results/fastp/{sample}/{sample}.html"),
        json      = "results/fastp/{sample}/{sample}.json",
    params:
        extra = config["fastp"]["extra"]
    log:
        "logs/fastp/{sample}.log"
    threads: 4
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        fastp --thread {threads} \
            {params.extra} \
            --in1 {input[0]} \
            --in2 {input[1]} \
            --out1 {output.trimmed[0]} \
            --out2 {output.trimmed[1]} \
            --json {output.json} \
            --html {output.html} 2>&1 | tee > {log}
        """

rule trim_nanopore:
    input:
        get_fastqs,
    output:
        trimmed = "results/fastplong/{sample}/{sample}.fastq.gz",
        html    = report("results/fastplong/{sample}/{sample}.html"),
        json    = "results/fastplong/{sample}/{sample}.json",
    params:
        extra = config["fastplong"]["extra"]
    log:
        "logs/fastplong/{sample}.log"
    threads: 4
    conda:
        "../envs/fastplong.yaml"
    shell:
        """
        fastplong --thread {threads} \
            {params.extra} \
            --in {input[0]} \
            --out {output.trimmed} \
            --json {output.json} \
            --html {output.html} 2>&1 | tee > {log}
        """
