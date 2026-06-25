"""
Remove human reads from a sequencing run with nohuman which is k-mer based instead of alignment based,
so it should be faster and uses less memory than hostile.
"""
rule remove_host_nohuman_illumina_pe:
    input:
        get_trimmed,
    output:
        r1 = "results/nohuman/{sample}/{sample}_R1.clean_1.fastq.gz",
        r2 = "results/nohuman/{sample}/{sample}_R2.clean_2.fastq.gz",
    params:
        out_cmd = lambda w, output: f"--out1 {output.r1} --out2 {output.r2}",
        nohuman_db = config["nohuman_db"]
    log:
        "logs/nohuman/{sample}.log"
    threads: 12
    conda:
        "../envs/nohuman.yaml"
    shell:
        """
        nohuman --threads {threads} \
            {params.out_cmd} \
            --output-type g \
            --db {params.nohuman_db} \
            {input} 2>&1 | tee > {log}
        """

use rule remove_host_nohuman_illumina_pe as remove_host_nohuman_nanopore with:
    input:
        get_trimmed,
    output:
        "results/nohuman/{sample}/{sample}.clean.fastq.gz"
    params:
        out_cmd = lambda w, output: f"--out1 {output}",
        nohuman_db = config["nohuman_db"]

"""
Use hostile if nohuman db is not available, otherwise use nohuman which is faster and uses less memory. 
"""
rule remove_host_hostile_illumina_pe:
    input:
        get_trimmed,
    output:
        r1 = "results/hostile/{sample}/{sample}_R1.clean_1.fastq.gz",
        r2 = "results/hostile/{sample}/{sample}_R2.clean_2.fastq.gz",
    params:
        read_cmd = lambda w, input: f"--fastq1 {input[0]} --fastq2 {input[1]}",
        outdir = "results/hostile/{sample}",
        index = config["hostile_minimap2_index"],
        aligner_args = "-x sr"
    log:
        "logs/hostile/{sample}.log"
    threads: 12
    conda:
        "../envs/hostile.yaml"
    shell:
        """
        hostile clean \
            --threads {threads} \
            {params.read_cmd} \
            --aligner minimap2 \
            --aligner-args "{params.aligner_args}" \
            --index {params.index} \
            --output {params.outdir} \
            --reorder \
            --airplane \
            --force 2>&1 | tee > {log}
        """

use rule remove_host_hostile_illumina_pe as remove_host_hostile_nanopore with:
    input:
        get_trimmed,
    output:
        "results/hostile/{sample}/{sample}.clean.fastq.gz",
    params:
        read_cmd = lambda w, input: f"--fastq1 {input[0]}",
        outdir = "results/hostile/{sample}",
        index = config["hostile_minimap2_index"],
        aligner_args = "-x map-ont"
