import gzip
import sys
min_length = 50
min_length_frac = 0.8
min_ident = 0.9




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
        while next_read_1[0] == curr_read_1[0]:
            read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = f1.readline().split(
                "\t")[:11]
            if next_read_1[0] == curr_read_1:
                
            next_read_1 = [read, read_len, subject, sub_start, matches, length]

        read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = f2.readline().split(
            "\t")[:11]
        next_read_2 = [read, read_len, subject, sub_start, matches, length]

        curr_read_1 = []

        while curr_read_1 == [] or curr_read_1[0][0] == next_read_1[0]:
            curr_read_1.append(next_read_1)
            read, read_len, read_start, read_end, strand, subject, sub_len, sub_start, sub_end, matches, length = f1.readline().split(
                "\t")[:]
            next_read_1 = [read, read_len, subject, sub_start, matches, length]

        curr_read_2 = next_read_2

        if curr_read1 = curr_read2:

    while True:
        while True:
            try:
                alignment_left = b1iter.__next__()
                if alignment_left.query_name == left_alignments[0].query_name:
                    left_alignments.append(alignment_left)
                else:
                    break
            except StopIteration:
                alignment_left = None
                break
        while True:
            try:
                alignment_right = b2iter.__next__()
                if alignment_right.query_name == right_alignments[0].query_name:
                    right_alignments.append(alignment_right)
                else:
                    break
            except StopIteration:
                alignment_right = None
                break
        #filter alignments
        ref_set_left = set()
        best_as = 0
        for i in left_alignments:
            if i.is_mapped:
                al = i.query_alignment_length
                mm =  i.get_tag("NM")
                ascore = i.get_tag("AS")
                if ascore > best_as:
                    hq_left_alignments = []
                    ref_set_left = set()
                    best_as = ascore
                if ascore == best_as and mm/al < 0.05 and al >= read_length * 0.9:
                    if not i.is_secondary or hq_left_alignments != []:
                        hq_left_alignments.append(i)
                        ref_set_left.add(i.reference_name)
        ref_set_right = set()
        best_as = 0
        for i in right_alignments:
            if i.is_mapped:
                al = i.query_alignment_length
                mm =  i.get_tag("NM")
                ascore = i.get_tag("AS")
                if ascore > best_as:
                    hq_right_alignments = []
                    ref_set_right = set()
                    best_as = ascore
                if ascore == best_as and mm/al < 0.05 and al >= read_length * 0.9:
                    if not i.is_secondary or hq_right_alignments != []:
                        hq_right_alignments.append(i)
                        ref_set_right.add(i.reference_name)
        valid_reference_names = ref_set_left.intersection(ref_set_right)
        if valid_reference_names == set():
            if alignment_left is None or alignment_right is None:
                sys.exit()
                break
            else:
                left_alignments = [alignment_left]
                right_alignments = [alignment_right]
                continue
        matched_left_alignments = []
        matched_right_alignments = []
        got_primary = True
        primary_seq, primary_qual = None, None
        for i in hq_left_alignments:
            if i.reference_name in valid_reference_names:
                if i.is_secondary and matched_left_alignments == []:
                    i.is_secondary = False
                    if i.cigartuples[0][0] == 5:
                        primary_seq = primary_seq[i.cigartuples[0][1]:]
                        primary_qual = list(primary_qual)[i.cigartuples[0][1]:]
                    if i.cigartuples[-1][0] == 5:
                        primary_seq = primary_seq[:-i.cigartuples[-1][1]]
                        primary_qual = list(primary_qual)[:-i.cigartuples[-1][1]]
                    i.query_sequence = primary_seq
                    i.query_qualities = primary_qual
                matched_left_alignments.append(i)
            if not i.is_secondary:
                primary_seq = i.query_sequence
                primary_qual = i.query_qualities
                if i.cigartuples[0][0] == 5:
                    primary_seq = "N" * i.cigartuples[0][1] + primary_seq
                    primary_qual = [0] * i.cigartuples[0][1] + list(primary_qual)
                if i.cigartuples[-1][0] == 5:
                    primary_seq += "N" * i.cigartuples[-1][1]
                    primary_qual = list(primary_qual) + [0] * i.cigartuples[-1][1]
        primary_seq, primary_qual = None, None
        for i in hq_right_alignments:
            if i.reference_name in valid_reference_names:
                if i.is_secondary and matched_left_alignments == []:
                    i.is_secondary = False
                    if i.cigartuples[0][0] == "H":
                        primary_seq = primary_seq[i.cigartuples[0][1]:]
                        primary_qual = list(primary_qual)[i.cigartuples[0][1]:]
                    if matched_left_alignments[0].cigartuples[-1][0] == "H":
                        primary_seq = primary_seq[:-i.cigartuples[-1][1]]
                        primary_qual = list(primary_qual)[:-i.cigartuples[-1][1]]
                    i.query_sequence = primary_seq
                    i.query_qualities = primary_qual
                matched_right_alignments.append(i)
            if not i.is_secondary:
                primary_seq = i.query_sequence
                primary_qual = i.query_qualities
                if i.cigartuples[0][0] == "H":
                    primary_seq = "N" * i.cigartuples[0][1] + primary_seq
                    primary_qual = [0] * i.cigartuples[0][1] + list(primary_qual)
                if i.cigartuples[-1][0] == "H":
                    primary_seq += "N" * i.cigartuples[-1][1]
                    primary_qual = list(primary_qual) + [0] * i.cigartuples[-1][1]
        for i in matched_left_alignments:
            i.is_paired = True
            i.is_read1 = True
            if i.query_name == "NB501781:467:HTGY5AFX5:1:11105:19538:10001":
                print(i)
            out_bam.write(i)
        for i in matched_right_alignments:
            i.is_paired = True
            i.is_read2 = True
            if i.query_name == "NB501781:467:HTGY5AFX5:1:11105:19538:10001":
                print(i)
            out_bam.write(i)
        if alignment_left is None or alignment_right is None:
            break
        else:
            left_alignments = [alignment_left]
            right_alignments = [alignment_right]
