configfile: workflow.source_path("config.yaml")
workdir: config["workdir"]
database_files: config["number_database_files"]

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")



rule all:
    input:
        expand("lbbc/{sample}.{ext}", sample=config["sample_name"], ext=["tblat.1", "grammy", "gi_tax_info.tab"])

rule trim_reads:
    params:
        read_dir = config["read_dir"]
    output:
        read1 = "qc_reads/{sample}_r1.fq.gz",
        read2 = "qc_reads/{sample}_r2.fq.gz"
    threads:
        config["threads"]
    run:
        import os, sys, subprocess
        read1, read2 = None, None
        for folder, subfolder, files in os.walk(params.read_dir):
            for f in files:
                print(f)
                if f.startswith(wildcards.sample) and f.endswith((".fq.gz", ".fq", ".fastq", ".fastq.gz")):
                    if "R1" in f.upper().lstrip(wildcards.sample) and "R2" in f.upper().lstrip(wildcards.sample):
                        sys.exit("Read file found that satisified both R1 and R2 regex")
                    elif "R1" in f.lstrip(wildcards.sample).upper():
                        read1 = os.path.join(folder, f)
                    elif "R2" in f.lstrip(wildcards.sample).upper():
                        read2 = os.path.join(folder, f)
                print(read1, read2)
        if read1 is None or read2 is None:
            sys.exit("Could't find read pairs in read directory.")
        subprocess.Popen("trimmomatic PE -threads {} {} {} {} /dev/null {} /dev/null "
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:100".format(threads, read1, read2, output.read1, output.read2), shell=True).wait()


rule remove_host:
    input:
        read1 = "qc_reads/{sample}_r1.fq.gz",
        read2 = "qc_reads/{sample}_r2.fq.gz"
    output:
        read1 = "qc_reads/{sample}_r1.clean_1.fastq.gz",
        read2 = "qc_reads/{sample}_r2.clean_2.fastq.gz"
    params:
        hostile_minimap2_index = config["hostile_minimap2_index"]
    threads:
        config["threads"]
    shell:
        "hostile clean --fastq1 {input.read1} --fastq2 {input.read2} --aligner minimap2 --index {params.hostile_minimap2_index}" \
        " --out-dir qc_reads --threads {threads}"




# takes a pair of read files and aligns them to a directory of
rule align:
    params:
        database_dir = config["database_dir"]
    input:
        read1 = "qc_reads/{sample}_r1.clean_1.fastq.gz",
        read2 = "qc_reads/{sample}_r2.clean_2.fastq.gz"
    output:
        paf1 = "alignments/{sample}.R1.{dbnum}.paf.gz",
        paf2 = "alignments/{sample}.R2.{dbnum}.paf.gz",
    threads:
        config["threads_align"]
    shell:
        "minimap2 -t {threads} -N 1000000 {params.database_dir}/cmdd.{dbnum}.bwa.fa.gz {input.read1} | gzip > {output.bam1} && "
        "minimap2 -t {threads} -N 1000000 {params.database_dir}/cmdd.{dbnum}.bwa.fa.gz {input.read2} | gzip > {output.bam2}"

rule merge_bams:
    input:
        bam1 = expand("alignments/{sample}.R1.{dbnum}.paf.gz", dbnum=[x for x in range(config["number_database_files"])]),
        bam2 = expand("alignments/{sample}.R2.{dbnum}.paf.gz", dbnum=[x for x in range(config["number_database_files"])])
    output:
        bam1 = "alignments/{sample}_R1.paf",
        bam2 = "alignments/{sample}_R2.paf"
    threads:
        config["threads"]
    run:
        shell("zcat {} | sort | gzip > {}".format(" ".join(input.bam1), output.bam1))
        shell("zcat {} | sort | gzip > {}".format(" ".join(input.bam2), output.bam2))


rule combine_bams:
    input:
        bam1 = "alignments/{sample}_R1.bam",
        bam2 = "alignments/{sample}_R2.bam"
    output:
        bam = "final_alignment/{sample}.bam"
    script:
        "scripts/combine.py"


