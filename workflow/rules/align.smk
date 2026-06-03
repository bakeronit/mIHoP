rule minimap2_align:
    input:
        get_host_removed,
    output:
        "results/minimap2/{sample}/{sample}.{dbnum}.paf.gz",
    params:
        extra = config.get("minimap", {}).get("extra", ""),
        preset = "-x sr" if config["platform"] == "illumina" else "-x map-ont",
        db = lambda w: f"{config['database_dir']}/cmdd.{w.dbnum}.fa.gz"
    log:
        "logs/minimap2/{sample}.{dbnum}.log",
    threads: 2
    shell:
        """
        minimap2 -t {threads} \
            {params.preset} \
            {params.extra} \
            {params.db} \
            {input} 2> {log} | gzip > {output}
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
    shell:
        """
        csvtk cut -t -f1-12 {input} | \
        csvtk add-header -H -t \
            --names "query_name,query_len,query_start,query_end,strand,target_name,target_len,target_start,target_end,matches,length,mapq" | \
        csvtk filter2 -t -f '$length > {params.min_length} && $matches / $length > {params.min_identity} && $length / $query_len > {params.min_coverage}' | \
        gzip > {output}
        """