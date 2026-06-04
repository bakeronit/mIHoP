import pandas as pd
from pathlib import Path
import os

from snakemake.utils import validate

validate(config, schema=workflow.source_path("../schemas/config.schema.yaml"))

samplesheet = config['samples'] if Path(config['samples']).is_absolute() else f"{os.environ['PWD']}/{config['samples']}"
sep = "\t" if samplesheet.endswith(".tsv") else ","
samples = pd.read_csv(samplesheet, sep=sep, comment="#").set_index(
    "sample",
    drop=False,
    verify_integrity=True
)

samples.index.names = ["sample_id"]

def drop_contant_cols(df):
    counts = df.nunique(dropna=True)
    constant_cols = counts[counts <= 1].index
    return df.drop(columns=constant_cols)

samples = drop_contant_cols(samples)
validate(samples, schema=workflow.source_path("../schemas/samples.schema.yaml"))

wildcard_constraints:
    sample="|".join(samples.index),
    dbnum="\\d+"

def is_single_end(sample):
    return pd.isna(samples.loc[sample, "fastq_2"])

def get_fastqs(wildcards):
    if is_single_end(wildcards.sample):
        return [samples.loc[wildcards.sample, "fastq_1"]]
    else:
        return [samples.loc[wildcards.sample, "fastq_1"], samples.loc[wildcards.sample, "fastq_2"]]

def get_trimmed(wildcards):
    if not is_single_end(wildcards.sample):
        return [
            f"results/fastp/{wildcards.sample}/{wildcards.sample}_R1.fastq.gz",
            f"results/fastp/{wildcards.sample}/{wildcards.sample}_R2.fastq.gz",
        ]
    ## NOTE: we currently only support Illumina (paired-end) and Nanopore platforms. 
    #elif config['platform'] == "illumina":
    #    return [f"results/fastp/{wildcards.sample}/{wildcards.sample}.fastq.gz"]
    elif config['platform'] == "nanopore":
        return [f"results/fastplong/{wildcards.sample}/{wildcards.sample}.fastq.gz"]
    else:
        raise ValueError(f"Unsupported platform: {config['platform']}")

def get_host_removed(wildcards):
    tool = "nohuman" if config.get("nohuman_db", None) else "hostile"
    if not is_single_end(wildcards.sample):
        return [
            f"results/{tool}/{wildcards.sample}/{wildcards.sample}_R1.clean_1.fastq.gz",
            f"results/{tool}/{wildcards.sample}/{wildcards.sample}_R2.clean_2.fastq.gz",
        ]
    elif config['platform'] in ["illumina", "nanopore"]:
        return [f"results/{tool}/{wildcards.sample}/{wildcards.sample}.clean.fastq.gz"]
    else:
        raise ValueError(f"Unsupported platform: {config['platform']}")