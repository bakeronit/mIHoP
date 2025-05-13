configfile: workflow.source_path("config.yaml")
workdir: config["workdir"]
database_files: config["number_database_files"]

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")



def get_samples(sample_name, sample_file):
    if sample_file != "none" and sample_name != "none":
        sys.exit("Please specify either sample_name OR sample_list")
    elif sample_file != "none":
        with open(sample_file) as f:
            sample_list = []
            for line in f:
                if line.rstrip() != "":
                    sample_list.append(line.rstrip())
    elif not sample_name != "none":
        sample_list = [sample_name]
    else:
        sys.exit("Please specify either sample_name or sample_list")
    return sample_list


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
        hostile_minimap2_index = config["hostile_minimap2_index"],
        nohuman_db = config["nohuman_database"]
    threads:
        config["threads"]
    run:
        if params.nohuman_db != "none":
            shell("nohuman -t {threads} -D {params.nohuman_db} --out1 {out.put.read1} --out2 {output.read2} {input.read1} {intput.read2}")
        else:
            shell("hostile clean --fastq1 {input.read1} --fastq2 {input.read2} --aligner minimap2 --index {params.hostile_minimap2_index}" \
        " --out-dir qc_reads --threads {threads}")





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
        "minimap2 -t {threads} -N 1000000 {params.database_dir}/cmdd.{wildcards.dbnum}.bwa.fa.gz {input.read1} | gzip > {output.paf1} && "
        "minimap2 -t {threads} -N 1000000 {params.database_dir}/cmdd.{wildcards.dbnum}.bwa.fa.gz {input.read2} | gzip > {output.paf2}"


def aggregate_input_R1(wildcards):
     return expand("alignments/{sample}.R1.{dbnum}.paf.gz", sample=wildcards.sample, dbnum=[x for x in range(config["number_database_files"])])

def aggregate_input_R2(wildcards):
     return expand("alignments/{sample}.R2.{dbnum}.paf.gz", sample=wildcards.sample, dbnum=[x for x in range(config["number_database_files"])])

rule merge_bams:
    input:
        paf1 = aggregate_input_R1,
        paf2 = aggregate_input_R2
    output:
        paf1 = "alignments/{sample}_R1.paf.gz",
        paf2 = "alignments/{sample}_R2.paf.gz"
    threads:
        config["threads"]
    run:
        shell("zcat {} | sort | gzip > {}".format(" ".join(input.paf1), output.paf1))
        shell("zcat {} | sort | gzip > {}".format(" ".join(input.paf2), output.paf2))


rule combine_pafs:
    input:
        paf1 = "alignments/{sample}_R1.paf.gz",
        paf2 = "alignments/{sample}_R2.paf.gz"
    params:
        taxonomy_file = config["taxonomy_file"]
    output:
        alignment = "final_alignment/{sample}.align.gz"
    script:
        "scripts/combine.py"

rule bin_alignments:
    input:
        alignment = "final_alignment/{sample}.align.gz" 
    output:
        bins = "bins/{sample}.bins"
    script:
        "scripts/bin.py"

rule process_bins:
    input:
        paf = expand("bins/{sample}.bins", sample=get_samples(config["sample_name"], config["sample_list"]))
