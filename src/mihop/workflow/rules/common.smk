import os
from pathlib import Path

import pandas as pd
from snakemake.utils import validate

# ------------------------------------------------------------------------------- #
# 1. Config
# ------------------------------------------------------------------------------- #

validate(config, schema=workflow.source_path("../schemas/config.schema.yaml"))

# ------------------------------------------------------------------------------- #
# 2. Samplesheet
# ------------------------------------------------------------------------------- #

samplesheet = config['samples'] if Path(config['samples']).is_absolute() else f"{os.environ['PWD']}/{config['samples']}"
sep = "\t" if samplesheet.endswith(".tsv") else ","
samples = pd.read_csv(samplesheet, sep=sep, comment="#").set_index(
    "sample",
    drop=False,
    verify_integrity=True
)
samples.index.names = ["sample_id"]

def drop_constant_cols(df: pd.DataFrame) -> pd.DataFrame:
    counts = df.nunique(dropna=True)
    constant_cols = counts[counts <= 1].index
    return df.drop(columns=constant_cols)

samples = drop_constant_cols(samples)
validate(samples, schema=workflow.source_path("../schemas/samples.schema.yaml"))

if config["platform"] == "illumina" and samples["fastq_2"].isna().any():
    raise ValueError("Currently this pipeline doesn't support single-end reads for Illumina platform. Please check your samplesheet and config.")

wildcard_constraints:
    sample="|".join(samples.index),
    dbnum="\\d+"

# ------------------------------------------------------------------------------- #
# 3. Helper functions for per-sample input
# ------------------------------------------------------------------------------- #

def is_single_end(samples: pd.DataFrame, sample: str) -> bool:
    return pd.isna(samples.loc[sample, "fastq_2"])

def get_fastqs(wildcards):
    if is_single_end(samples, wildcards.sample):
        return [samples.loc[wildcards.sample, "fastq_1"]]
    else:
        return [samples.loc[wildcards.sample, "fastq_1"], samples.loc[wildcards.sample, "fastq_2"]]

def get_trimmed(wildcards):
    if not is_single_end(samples, wildcards.sample):
        return [
            f"results/fastp/{wildcards.sample}/{wildcards.sample}_R1.fastq.gz",
            f"results/fastp/{wildcards.sample}/{wildcards.sample}_R2.fastq.gz",
        ]
    elif config['platform'] == "nanopore":
        return [f"results/fastplong/{wildcards.sample}/{wildcards.sample}.fastq.gz"]
    else:
        raise ValueError(f"Unsupported platform: {config['platform']}")

def get_host_removed(wildcards):
    tool = "nohuman" if config.get("nohuman_db", None) else "hostile"
    if not is_single_end(samples, wildcards.sample):
        return [
            f"results/{tool}/{wildcards.sample}/{wildcards.sample}_R1.clean_1.fastq.gz",
            f"results/{tool}/{wildcards.sample}/{wildcards.sample}_R2.clean_2.fastq.gz",
        ]
    elif config['platform'] in ["illumina", "nanopore"]:
        return [f"results/{tool}/{wildcards.sample}/{wildcards.sample}.clean.fastq.gz"]
    else:
        raise ValueError(f"Unsupported platform: {config['platform']}")