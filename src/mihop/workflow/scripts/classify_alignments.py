import gzip
import csv
import sys
from collections import defaultdict
from dataclasses import dataclass
import pandas as pd


@dataclass
class ReadAlignment:
    read_name: str
    taxon: str
    best_hit_name: str
    best_hit_start: int
    matches: int
    is_core: bool
    is_repeat: bool
    is_shared: bool
    is_unmatched: bool

def load_taxonomy(taxonomy_file):
    taxonomy_of_accession = {}
    accession_counts_per_tax = defaultdict(int)
    with open(taxonomy_file) as f:
        for line in f:
            accession, taxonomy = line.rstrip().split("\t")
            taxonomy_of_accession[accession] = taxonomy
            accession_counts_per_tax[taxonomy] += 1
    return taxonomy_of_accession, accession_counts_per_tax

def group_alignments(alignments: list, taxonomy_of_accession: dict, accession_counts_per_tax: defaultdict[int], paired: bool=True) -> ReadAlignment | None:
    is_core = is_unmatched = is_shared = False
    has_left = any(row['query_name'].endswith("/1") for row in alignments)
    has_right = any(row['query_name'].endswith("/2") for row in alignments)
    if paired and not (has_left and has_right):
        return None  #NOTE: for paired-end reads, both ends have to be aligned.
    
    df = pd.DataFrame(alignments)
    df['accession'] = df['target_name'].apply(lambda x: x.split("|")[0])
    df['taxonomy'] = df['accession'].apply(lambda acc: taxonomy_of_accession.get(acc, None))
    df = df[df['taxonomy'].notnull()]  # filter out alignments to accessions that are not in the taxonomy file.
    if df.empty:
        return None
    accessions = set(df['accession'].tolist())
    taxonomies = set(df['taxonomy'].tolist())

    if len(taxonomies) == 1:
        is_core = len(accessions) == accession_counts_per_tax[list(taxonomies)[0]]
    elif len(taxonomies) > 1:
        is_shared = True
    else:
        is_unmatched = True  ## NOTE:only if it aligns to something that is not in the taxonomy file.

    if paired:
        df_r1 = df[df['query_name'].str.endswith("/1")]
        accessions_r1 = set(df_r1['accession'].tolist())
        df_r2 = df[df['query_name'].str.endswith("/2")]
        accessions_r2 = set(df_r2['accession'].tolist())
        is_repeat = len(accessions_r1) < len(df_r1) and len(accessions_r2) < len(df_r2)
        if accessions_r1.isdisjoint(accessions_r2) and len(taxonomies) > 1:
            print(f"Warning: read {df.iloc[0]['query_name'].replace('/1', '')} has no shared accessions between R1 and R2", file=sys.stderr)
    else:
        is_repeat = len(accessions) < len(df)

    # NOTE: get the best hit by the highest identity, and then the highest alignment block length, then mapq.
    # still might have multiple best hits, try to pick the first one in a deterministic way (by target_name).
    #df['identity'] = df['matches'].astype(int) / df['length'].astype(int)
    #df_sorted = df.sort_values(by=['identity', 'length', 'mapq', 'target_name'], ascending=False)  #NOTE: it might be better to sort based on num_of_matches, length, then mapq, without calculating identity, to avoid issues with very short alignments having high identity. 
    df = df.astype({'matches': int, 'length': int, 'mapq': int})
    df_sorted = df.sort_values(by=['matches', 'length', 'mapq', 'target_name'], ascending=[False, False, False, True])
    best_hit = df_sorted.iloc[0]

    read_alignment = ReadAlignment(
        read_name=best_hit['query_name'].replace("/1", "").replace("/2", ""),
        taxon=','.join(taxonomies),
        best_hit_name=best_hit['target_name'],
        best_hit_start=int(best_hit['target_start']),
        matches=len(accessions),
        is_core=is_core,
        is_repeat=is_repeat,
        is_shared=is_shared,
        is_unmatched=is_unmatched
    )

    return read_alignment


def main():

    # to make this code work for both snakemake and CLI.
    if "snakemake" in globals():
        print("Running via Snakemake...", file=sys.stderr)
        taxonomy_file = snakemake.input.taxonomy
        paf = snakemake.input.paf
        output = snakemake.output.tsv
        paired = snakemake.params.paired
    else:
        print("Running as a standalone CLI script...", file=sys.stderr)
        import argparse
        parser = argparse.ArgumentParser(description="Classify reads based on their alignments to reference genomes.")
        parser.add_argument("-t", "--taxonomy_file", required=True, help="Path to the taxonomy file mapping accessions to taxonomies.")
        parser.add_argument("-p", "--paf", required=True, help="Path to the input PAF file containing read alignments.")
        parser.add_argument("--paired", action='store_true', help="Whether the input reads are paired-end (default: False).")
        parser.add_argument("-o", "--output", required=True, help="Path to the output file where classified reads will be written.")
        args = parser.parse_args()

        taxonomy_file = args.taxonomy_file
        paf = args.paf
        output = args.output
        paired = args.paired

    taxonomy_of_accession, accession_counts_per_tax = load_taxonomy(taxonomy_file)
    read_alignments = []
    with gzip.open(paf, 'rt') as fh, gzip.open(output, 'wt') if output.endswith('.gz') else open(output, 'wt') as out:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["read_name", "taxonomies", "best_hit_name", "best_hit_start", "matches", "core", "repeat", "shared", "unmatched"])

        current_read_base = None
        alignments = []
        for row in reader:
            read = row['query_name']
            if paired:
                if not read.endswith("/1") and not read.endswith("/2"):
                    raise ValueError(f"Read {read} does not end with /1 or /2, but paired_mode is True")
            read_base = read.replace("/1", "").replace("/2", "")
            if current_read_base is None:
                alignments.append(row)
                current_read_base = read_base
            elif current_read_base == read_base:
                alignments.append(row)
            else:
                read_alignment = group_alignments(alignments, taxonomy_of_accession, accession_counts_per_tax, paired)
                if read_alignment:
                    writer.writerow([
                        read_alignment.read_name,
                        read_alignment.taxon,
                        read_alignment.best_hit_name,
                        read_alignment.best_hit_start,
                        read_alignment.matches,
                        read_alignment.is_core,
                        read_alignment.is_repeat,
                        read_alignment.is_shared,
                        read_alignment.is_unmatched
                    ])
                alignments = [row]
                current_read_base = read_base
        
        # process alignments of the last read
        read_alignment = group_alignments(alignments, taxonomy_of_accession, accession_counts_per_tax, paired)
        if read_alignment:
            writer.writerow([
                read_alignment.read_name,
                read_alignment.taxon,
                read_alignment.best_hit_name,
                read_alignment.best_hit_start,
                read_alignment.matches,
                read_alignment.is_core,
                read_alignment.is_repeat,
                read_alignment.is_shared,
                read_alignment.is_unmatched
            ])

if __name__ == "__main__":
    main()
