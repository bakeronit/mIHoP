import gzip
import sys
min_length = 50
min_length_frac = 0.8
min_ident = 0.9



with open(snakemake.params.taxonomy_file) as f:
    for line in f:
        accession, tax = line.rstrip().split("\t")
        tax_dict[accession] = tax
        accession_dict[tax] = accession


with gzip.open(snakemake.output.paf, 'wt') as o, gzip.open(snakemake.input.paf1) as f1, gzip.open(snakemake.input.paf2) as f2:
    read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = f1.readline().split("\t")[:11]
    next_read_1 = [read, read_len, subject, sub_start, matches, length]
    read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = f2.readline().split(
        "\t")[:11]
    next_read_2 = [read, read_len, subject, sub_start, matches, length]
    while True:
        curr_read_1 = next_read_1
        curr_read_2 = next_read_2
        while curr_read_1[0] < curr_read_2[0]:
            read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = f1.readline().split(
                "\t")[:11]
            curr_read_1 = [read, read_len, subject, sub_start, matches, length]
        while curr_read_2[0] < curr_read_1[0]:
            read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = f2.readline().split(
                "\t")[:11]
            curr_read_2 = [read, read_len, subject, sub_start, matches, length]
        additional_hits_1 = []
        while next_read_1[0] == curr_read_1[0]:
            read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = f1.readline().split(
                "\t")[:11]
            next_read_1 = [read, read_len, subject, sub_start, matches, length]
            if next_read_1[0] == curr_read_1:
                additional_hits_1.append(next_read_1)
        additional_hits_2 = []
        while next_read_2[0] == curr_read_2[0]:
            read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = f2.readline().split(
                "\t")[:11]
            next_read_2 = [read, read_len, subject, sub_start, matches, length]
            if next_read_2[0] == curr_read_2:
                additional_hits_2.append(next_read_1)


