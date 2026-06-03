import gzip
import csv
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


def group_alignments(alignments: list, taxonomy_of_accession: dict, accession_counts_per_tax: defaultdict[int], paired: bool=True):
    is_core = is_unmatched = is_shared = False
    has_left = any(row['query_name'].endswith("/1") for row in alignments)
    has_right = any(row['query_name'].endswith("/2") for row in alignments)
    if paired and not (has_left and has_right):
        return False, None  # don't keep this read, if not both mates aligned.
    
    df = pd.DataFrame(alignments)
    accessions = [t.split("|")[0] for t in df['target_name'].to_list() ]
    taxonomies = [ taxonomy_of_accession.get(acc, None) for acc in accessions ]
    taxonomies = [ tax for tax in taxonomies if tax is not None ]

    if len(set(taxonomies)) == 1:
        is_core = len(set(accessions)) == accession_counts_per_tax[taxonomies[0]]
    elif len(set(taxonomies)) > 1:
        is_shared = True
    else:
        is_unmatched = True  ## dead code?

    if paired:
        df_r1 = df[df['query_name'].str.endswith("/1")]
        accessions_r1 = [ t.split("|")[0] for t in df_r1['target_name'].to_list() ]
        df_r2 = df[df['query_name'].str.endswith("/2")]
        accessions_r2 = [ t.split("|")[0] for t in df_r2['target_name'].to_list() ]
        is_repeat = len(set(accessions_r1)) < len(accessions_r1) and len(set(accessions_r2)) < len(accessions_r2)
    else:
        is_repeat = len(set(accessions)) < len(accessions)

    # get the best hit by the highest identity, and then the highest alignment block length, then mapq.
    # still might have multiple best hits, try to pick the first one in a deterministic way (by target_name).
    df['identity'] = df['matches'].astype(int) / df['length'].astype(int)
    df_sorted = df.sort_values(by=['identity', 'length', 'mapq', 'target_name'], ascending=False)
    best_hit = df_sorted.iloc[0]

    read_alignment = ReadAlignment(
        read_name=best_hit['query_name'].replace("/1", "").replace("/2", ""),
        taxon=','.join(set(taxonomies)),
        best_hit_name=best_hit['target_name'],
        best_hit_start=int(best_hit['target_start']),
        matches=len(accessions),
        is_core=is_core,
        is_repeat=is_repeat,
        is_shared=is_shared,
        is_unmatched=is_unmatched
    )

    return True, read_alignment


def main():

    # to make this code work for both snakemake and CLI.
    if "snakemake" in globals():
        print("Running via Snakemake...")
        taxonomy_file = snakemake.input.taxonomy
        paf = snakemake.input.paf
        output = snakemake.output.tsv
        paired = snakemake.params.paired
    else:
        print("Running as a standalone CLI script...")
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
    
    if not output.endswith(".gz"):
        output += ".gz"

    taxonomy_of_accession, accession_counts_per_tax = load_taxonomy(taxonomy_file)
    read_alignments = []
    with gzip.open(paf, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter="\t")
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
                keep, read_alignment = group_alignments(alignments, taxonomy_of_accession, accession_counts_per_tax, paired)
                if keep:
                    read_alignments.append(read_alignment)
                alignments = [row]
                current_read_base = read_base
        
        #alignments of the last read
        keep, read_alignment = group_alignments(alignments, taxonomy_of_accession, accession_counts_per_tax, paired)
        if keep:
            read_alignments.append(read_alignment)

    with gzip.open(output, 'wt' ) as o:
        writer = csv.writer(o, delimiter="\t")
        writer.writerow(["read_name", "taxonomies", "best_hit_name", "best_hit_start", "matches", "core", "repeat", "shared", "unmatched"])
        for rl in read_alignments:
            writer.writerow([
                rl.read_name, 
                rl.taxon, 
                rl.best_hit_name,
                rl.best_hit_start,
                rl.matches, 
                rl.is_core, 
                rl.is_repeat, 
                rl.is_shared, 
                rl.is_unmatched
            ])

if __name__ == "__main__":
    main()