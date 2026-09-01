# SPDX-FileCopyrightText: 2026-present QPHIRE-genomics <Mitchell.Sullivan@health.qld.gov.au>
#
# SPDX-License-Identifier: MIT

from pathlib import Path

from snk_cli import CLI

mihop = CLI(Path(__file__).parent.parent)

@mihop.app.command()
def metadb(directory: str|Path = "."):
    """
    Download clinical meta DB from s3
    """
    s3_url = "s3://xxxxx"
    print(f"Downloading clinical meta DB from {s3_url}..., to {directory}")

@mihop.app.command()
def prepare(fastq_dir: str|Path, output: str|Path = "samples.csv") -> Path:
    """
    Generate a samplesheet from FASTQ files found in `fastq_dir`.
    """
    fastq_dir, output = Path(fastq_dir), Path(output)
    if not fastq_dir.is_dir():
        raise FileNotFoundError(f"FASTQ directory does not exist: {fastq_dir}")

    suffixes = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
    sep = "\t" if output.suffix.lower() == ".tsv" else ","

    samples: dict[str, list[Path]] = {}
    for path in fastq_dir.rglob("*"):
        name = path.name.lower()
        if name.endswith(suffixes) and not name.startswith("undetermined") and path.is_file():
            samples.setdefault(path.name.split("_", 1)[0], []).append(path)

    with output.open("w", newline="") as out:
        out.write(sep.join(("sample", "fastq_1", "fastq_2")) + "\n")
        for sample, reads in sorted(samples.items()):
            reads = sorted(map(str, reads))
            if len(reads) > 2:
                raise ValueError(f"Expected at most two FASTQs for {sample!r}, found: {reads}")
            r1, r2, *_ = *reads, "", ""
            out.write(sep.join((sample, r1, r2)) + "\n")

    return output