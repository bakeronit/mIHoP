rule minimap2_align:
    input:
        get_host_removed,
    output:
        temp("results/minimap2/{sample}/{sample}.{dbnum}.paf.gz"),
    params:
        extra = config.get("minimap", {}).get("extra", ""),
        preset = "-x sr" if config["platform"] == "illumina" else "-x map-ont",
        fastq_in = lambda w, input: f"-1 {input[0]} -2 {input[1]}" if len(input) == 2 else f"{input[0]}",
        db = lambda w: f"{config['database_dir']}/cmdd.{w.dbnum}.fa.gz"
    log:
        "logs/minimap2/{sample}.{dbnum}.log",
    threads: 2
    conda:
        "../envs/align.yaml"
    shell: # use samtools import | fastq to ensure the /1 and /2 for paired-end reads, which supports minimap2 <= 2.28
        """
        samtools import -T'*' {params.fastq_in} 2> {log} | \
        samtools fastq 2>> {log} | \
        minimap2 -t {threads} \
            {params.preset} \
            {params.extra} \
            {params.db} \
            - 2>> {log} | gzip > {output}
        """

rule merge_paf:
    input:
        (f"results/minimap2/{{sample}}/{{sample}}.{dbnum}.paf.gz" for dbnum in range(config['n_database_files'])),
    output:
        "results/minimap2/{sample}/{sample}.all.sorted.paf.gz",
    shell:
        """
        zcat {input} | sort | gzip > {output}
        """

rule filter_paf:
    input:
        "results/minimap2/{sample}/{sample}.all.sorted.paf.gz",
    output:
        "results/minimap2/{sample}/{sample}.all.sorted.filtered.paf.gz",
    params:
        min_length   = config.get("minimap", {}).get("min_length", 50),
        min_coverage = config.get("minimap", {}).get("min_coverage", 0.8),
        min_identity = config.get("minimap", {}).get("min_identity", 0.9)
    conda:
        "../envs/align.yaml"
    shell:
        """
        csvtk cut -t -f1-12 {input} | \
        csvtk add-header -H -t \
            --names "query_name,query_len,query_start,query_end,strand,target_name,target_len,target_start,target_end,matches,length,mapq" | \
        csvtk filter2 -t -f '$length > {params.min_length} && $matches / $length > {params.min_identity} && $length / $query_len > {params.min_coverage}' -o {output}
        """
