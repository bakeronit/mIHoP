import sys
import gzip
from collections import defaultdict




alignment_file = snakemake.input.alignment
mihop_out = snakemake.output.bins

#alignment_file = sys.argv[2]
#mihop_out = sys.argv[3]

tax_dict = {}
full_tax_dict = {}




bin_size = 500


with gzip.open(alignment_file, 'rt') as f:
    bins = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0))) 
    for line in f:
        read, species, chrom, pos, genomes_hit, core, repeat, shared, unmatched = line.rstrip().split("\t")
        pos = int(pos)
        if core == "T" and repeat == "F" and shared == "F" and unmatched == "F":
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




