import sys
import gzip
from collections import defaultdict

if "snakemake" in globals():
    print("Running via Snakemake...", file=sys.stderr)
    alignment_file = snakemake.input.tsv
    mihop_out = snakemake.output.bins
    bin_size = snakemake.params.bin_size
else:
    print("Running as a standalone CLI script...", file=sys.stderr)
    alignment_file = sys.argv[1]
    mihop_out = sys.argv[2]
    bin_size = int(sys.argv[3])

if bin_size <= 0:
    raise ValueError("Bin size must be a positive integer.")


with gzip.open(alignment_file, 'rt') if alignment_file.endswith('.gz') else open(alignment_file, 'rt') as f:
    f.readline()
    bins = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0))) 
    for line in f:
        read, species, chrom, pos, genomes_hit, core, repeat, shared, unmatched = line.rstrip().split("\t")
        pos = int(pos)
        if repeat == "False" and shared == "False" and unmatched == "False":
            bins[species][chrom][pos//bin_size] += 1

outlist = []
for species in bins:
    total_bins = 0
    for chrom in bins[species]:
        for the_bin in bins[species][chrom]:
            total_bins += 1
    outlist.append([species.split(';')[-1], species, total_bins])

outlist.sort(key=lambda x: x[2], reverse=True)
with open(mihop_out, 'w') as o:
    for i in outlist:
        o.write("\t".join([i[0], i[1], str(i[2])]) + "\n")