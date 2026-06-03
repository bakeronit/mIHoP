import gzip
import sys
from collections import defaultdict
min_length = 50
min_length_frac = 0.8
min_ident = 0.9

tax_dict = {}
genome_counts = defaultdict(lambda:0)

with open(snakemake.params.taxonomy_file) as f:
    for line in f:
        accession, tax = line.rstrip().split("\t")
        tax_dict[accession] = tax
        genome_counts[tax] += 1

with gzip.open(snakemake.output.alignment, 'wt') as o, gzip.open(snakemake.input.paf1, 'rt') as f1, gzip.open(snakemake.input.paf2, 'rt') as f2:
    # get first L read that passes min thresholds
    while True:
        line = f1.readline()
        if line == '':
            next_read_1 = None
            break
        read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = line.split("\t")[:11]
        if read.endswith('/1'):
            read = read[:-2]
        matches, read_len, length, sub_start = int(matches), int(read_len), int(length), int(sub_start)
        if matches / length >= min_ident and length / read_len >= min_length_frac and length > min_length:
            next_read_1 = [read, read_len, subject, sub_start, matches, length]
            break
    # get first R read that passes min thresholds
    while True:
        line = f2.readline()
        if line == '':
            next_read_2 = None
            break
        read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = line.split("\t")[:11]
        if read.endswith('/2'):
            read = read[:-2]
        matches, read_len, length, sub_start = int(matches), int(read_len), int(length), int(sub_start)
        if matches / length >= min_ident and length / read_len >= min_length_frac and length > min_length:
            next_read_2 = [read, read_len, subject, sub_start, matches, length]
            break
    # read the rest of the alignment file
    while not next_read_1 is None and not next_read_2 is None:
        curr_read_1 = next_read_1
        curr_read_2 = next_read_2
        # If we are behind in L reads
        while curr_read_1[0] < curr_read_2[0]:
            line = f1.readline()
            if line == '':
                next_read_1 = None
                break
            read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = line.split(
                "\t")[:11]
            if read.endswith('/1'):
                read = read[:-2]
            matches, read_len, length, sub_start = int(matches), int(read_len), int(length), int(sub_start)
            if matches / length >= min_ident and length / read_len >= min_length_frac and length > min_length:
                curr_read_1 = [read, read_len, subject, sub_start, matches, length]

        # if we are behind in R reads
        while curr_read_2[0] < curr_read_1[0]:
            line = f2.readline()
            if line == '':
                next_read_2 = None
                break
            read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = line.split(
                "\t")[:11]
            if read.endswith('/2'):
                read = read[:-2]
            matches, read_len, length, sub_start = int(matches), int(read_len), int(length), int(sub_start)
            if matches / length >= min_ident and length / read_len >= min_length_frac and length > min_length:
                curr_read_2 = [read, read_len, subject, sub_start, matches, length]
        #if we're at the end of one of the alignment files and didn't find a match
        if next_read_1 is None or next_read_2 is None:
            break
        next_read_1 = curr_read_1
        next_read_2 = curr_read_2

        if curr_read_2[0] > curr_read_1[0]:
            # if the R read got in front of the L read search forward in the L reads
            continue
        #print(curr_read_1[0], curr_read_2[0])
        # add the additional hits
        additional_hits_1 = []
        while next_read_1[0] == curr_read_1[0]:
            line = f1.readline()
            if line == '':
                next_read_1 = None
                break
            read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = line.split(
                "\t")[:11]
            if read.endswith('/1'):
                read = read[:-2]
            read_len, sub_start, matches, length = int(read_len), int(sub_start), int(matches), int(length)
            if matches / length >= min_ident and length / read_len >= min_length_frac and length > min_length:
                next_read_1 = [read, read_len, subject, sub_start, matches, length]
                if next_read_1[0] == curr_read_1: #REVIEW: comparing a string to list? addtional_hits_1 will always be empty
                    additional_hits_1.append(next_read_1)
        additional_hits_2 = []
        while next_read_2[0] == curr_read_2[0]:
            line = f2.readline()
            if line == '':
                next_read_2 = None
                break
            read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = line.split(
                "\t")[:11]
            if read.endswith('/2'):
                read = read[:-2]
            read_len, sub_start, matches, length = int(read_len), int(sub_start), int(matches), int(length)
            if matches / length >= min_ident and length / read_len >= min_length_frac and length > min_length:
                next_read_2 = [read, read_len, subject, sub_start, matches, length]
                if next_read_2[0] == curr_read_2:  #REVIEW: same as above, comparing a string to list? addtional_hits_2 will always be empty
                    additional_hits_2.append(next_read_2)
        # Process the alignments
        best_hit = [curr_read_1]
        accession_set_1 = set([curr_read_1[2].split("|")[0]])  
        for i in additional_hits_1:                                 
            read, read_len, subject, sub_start, matches, length = i
            accession = subject.split("|")[0]
            if matches/length > best_hit[0][4]/best_hit[0][5]:
                best_hit = [i]
            elif matches/length == best_hit[0][4]/best_hit[0][5] and length > best_hit[0][5]:
                best_hit = [i]
            elif matches/length == best_hit[0][4]/best_hit[0][5] and length == best_hit[0][5]:
                best_hit.append(i)
            accession_set_1.add(accession)
        best_hit_1 = best_hit
        best_hit = [curr_read_2]  #REVIEW: overwriting best_hit from the L read with the R read, but did nothing with the best_hit_1 variable that was just created
        accession_set_2 = set([curr_read_2[2].split("|")[0]])
        for i in additional_hits_2:
            read, read_len, subject, sub_start, matches, length = i
            accession = subject.split("|")[0]
            if matches / length > best_hit[0][4] / best_hit[0][5]:
                best_hit = [i]
            elif matches / length == best_hit[0][4] / best_hit[0][5] and length > best_hit[0][5]:
                best_hit = [i]
            elif matches / length == best_hit[0][4] / best_hit[0][5] and length == best_hit[0][5]:
                best_hit.append(i)
            accession_set_2.add(accession)
        bets_hit_2 = best_hit
        if len(accession_set_1) == len(additional_hits_1) + 1 or len(accession_set_2) == len(additional_hits_2) + 1:
            repeat = "F"     #REVIEW: here the repeat will always be "F" because of the bugs in the code above where additional_hits_1 and additional_hits_2 will always be empty, so len(accession_set_1) will always be 1 more than len(addtional_hits_1) and same for 2. So repeat will always be "F" and never "T"
        else:
            repeat = "T"
        taxonomies = []
        for i in accession_set_1.union(accession_set_2):
            taxonomies.append(tax_dict[i])
        core = "F"
        unmatched = "F"
        if len(set(taxonomies)) == 1:
            if len(taxonomies) == genome_counts[taxonomies[0]]:
                core = "T"
            shared = "F"
        elif len(set(taxonomies)) > 1:
            shared = "T"
        else:
            unmatched = "T"  #REVIEW: this shouldn't happend? what is unmatched? 
        tax_out = []
        for i in set(taxonomies):
            tax_out.append(i)
        tax_out = ','.join(tax_out)
        o.write("\t".join([curr_read_1[0], tax_out, best_hit[0][2], str(best_hit[0][3]), str(len(additional_hits_1)+1), core, repeat, shared, unmatched]) + "\n")  #REVIEW: why the length of hits of r1 is records as number of hits?
        










            
            


