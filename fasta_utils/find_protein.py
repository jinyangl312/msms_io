import re


def find_seq_in_fasta_db(pattern, offset, fm_decoder):
    res = []
    for (fasta_name, seq_offset) in fm_decoder.get_all_pattern(pattern):
        res.append(f"{fasta_name} ({offset+seq_offset})")
    return res


def find_protein_in_fasta_db(pattern, fm_decoder):
    res = []
    for (fasta_name, seq_offset) in fm_decoder.get_all_pattern(pattern):
        res.append(f"{fasta_name}")
    return res


def find_xl_seq_in_fasta_db(mixed_peptide, fm_decoder):
    line = re.split("\-|\(|\)", mixed_peptide)
    seq1, site1, seq2, site2 = line[0], int(line[1]), line[3], int(line[4])

    res1_list = find_seq_in_fasta_db(seq1, site1, fm_decoder)
    res2_list = find_seq_in_fasta_db(seq2, site2, fm_decoder)

    if len(res1_list) == 0 or len(res2_list) == 0:
        return ""

    res = []
    has_target = False
    for res1 in res1_list:
        for res2 in res2_list:
            if not "REV" in res1 and not "REV" in res2:
                # Only keep target if target exists
                has_target = True
                res.append(f"{res1}-{res2}/")
            elif not has_target:
                # Keep decoy only if target not exist
                res.append(f"{res1}-{res2}/")
    return "".join(res)


def contains_subseq_in_fasta_db(pattern, fm_decoder):
    for position in range(len(pattern)):
        new_pattern = pattern[0:position]+pattern[position+1:len(pattern)]
        if len(fm_decoder.get_all_pattern(new_pattern)) > 0:
            return True
    return False


def is_xl_seq_in_same_syn_group(mixed_peptide, fm_decoder):
    line = re.split("\-|\(|\)", mixed_peptide)
    seq1, site1, seq2, site2 = line[0], int(line[1]), line[3], int(line[4])

    res1_list = find_protein_in_fasta_db(seq1, fm_decoder)
    res2_list = find_protein_in_fasta_db(seq2, fm_decoder)

    dest_set = set(res2_list)
    for src_item in res1_list:
        if src_item in dest_set:
            return 2 # same group
    if len(res1_list) and len(res2_list):
        return 1 # different group
    else:
        if len(res1_list) and contains_subseq_in_fasta_db(seq2, fm_decoder):
            return 0 # similar peptides
        elif len(res2_list) and contains_subseq_in_fasta_db(seq1, fm_decoder):
            return 0 # similar peptides
        return -1 # Not in syn group
