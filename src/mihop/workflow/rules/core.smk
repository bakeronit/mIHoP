
rule classify_alignments:
    input:
        taxonomy = config["taxonomy_file"],
        paf = "results/minimap2/{sample}/{sample}.all.sorted.filtered.paf.gz"
    output:
        tsv = "results/mihop_core/{sample}/{sample}.alignments.tsv.gz"
    params:
        paired = config.get("platform", "") == "illumina" # Assume paired-end if platform is Illumina, single-end otherwise
    threads: 1
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/classify_alignments.py"


rule genome_bin:
    input:
        tsv = "results/mihop_core/{sample}/{sample}.alignments.tsv.gz"
    output:
        bins = "results/mihop_core/{sample}/{sample}.bins"
    params:
        bin_size = config.get("bin_size", 500)
    threads: 1
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/bin.py"