import re
import pathlib
from msms_io.fasta.fasta_toolbox import encode_esa, DecoderEsa
from pyteomics import fasta

def find_seq_in_fasta_db(pattern, offset, decoder):
    res = []
    for (fasta_name, seq_offset) in decoder.get_all_pattern(pattern):
        res.append(f"{fasta_name} ({offset+seq_offset})")
    return res


def find_protein_in_fasta_db(pattern, decoder):
    res = []
    for (fasta_name, seq_offset) in decoder.get_all_pattern(pattern):
        res.append(f"{fasta_name}")
    return res


def find_xl_seq_in_fasta_db(mixed_peptide, decoder):
    line = re.split("\-|\(|\)", mixed_peptide)
    seq1, site1, seq2, site2 = line[0], int(line[1]), line[3], int(line[4])

    res1_list = find_seq_in_fasta_db(seq1, site1, decoder)
    res2_list = find_seq_in_fasta_db(seq2, site2, decoder)

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


def contains_subseq_in_fasta_db(pattern, decoder):
    for position in range(len(pattern)):
        new_pattern = pattern[0:position]+pattern[position+1:len(pattern)]
        if len(decoder.get_all_pattern(new_pattern)) > 0:
            return True
    return False


def is_xl_seq_in_syn_group(mixed_peptide, decoder):
    line = re.split("\-|\(|\)", mixed_peptide)
    seq1, site1, seq2, site2 = line[0], int(line[1]), line[3], int(line[4])

    res1_list = find_protein_in_fasta_db(seq1, decoder)
    res2_list = find_protein_in_fasta_db(seq2, decoder)

    return 2 if len(res1_list) and len(res2_list) else -1


def is_xl_seq_in_same_syn_group(mixed_peptide, decoder):
    line = re.split("\-|\(|\)", mixed_peptide)
    seq1, site1, seq2, site2 = line[0], int(line[1]), line[3], int(line[4])

    res1_list = find_protein_in_fasta_db(seq1, decoder)
    res2_list = find_protein_in_fasta_db(seq2, decoder)

    dest_set = set(res2_list)
    for src_item in res1_list:
        if src_item in dest_set:
            return 2  # same group
    if len(res1_list) and len(res2_list):
        return 1  # different group
    else:
        if len(res1_list) and contains_subseq_in_fasta_db(seq2, decoder):
            return 0  # similar peptides
        elif len(res2_list) and contains_subseq_in_fasta_db(seq1, decoder):
            return 0  # similar peptides
        return -1  # Not in syn group


def is_xl_protein_in_same_syn_group(mixed_protein, syn_protein_dict):
    line = re.sub("\(\d+\)", "", mixed_protein)
    line = re.split("\-", line)
    protein_1, protein_2 = line[0].strip(), line[1].strip()
    protein_1 = re.search('(?<=\|).*(?=\|)', protein_1).group()
    protein_2 = re.search('(?<=\|).*(?=\|)', protein_2).group()

    try:
        res1_list = syn_protein_dict[protein_1]
        res2_list = syn_protein_dict[protein_2]
    except:
        return -1  # Not in syn group

    dest_set = set(res2_list)
    for src_item in res1_list:
        if src_item in dest_set:
            return 2  # same group
    return 1  # different group


def get_decoder_I2L_rev(fasta_path):
    index_prefix = fasta_path.parent/'.index.tmp'/(fasta_path.stem+'.I2L.rev.esa')/'index'
    if not pathlib.Path(index_prefix).parent.exists():
        index_prefix.parent.mkdir(exist_ok=True, parents=True)
        with fasta.read(str(fasta_path)) as db:
            transferred_data = [
                (item.description, item.sequence.replace("I", "L")) for item in db]
        fasta.write(transferred_data, str(fasta_path.with_suffix(
            ".I2L.fasta")), file_mode="w")
        fasta.write_decoy_db(str(fasta_path.with_suffix(".I2L.fasta")),
                                str(fasta_path.with_suffix(".I2L.rev.fasta")),
            prefix="REV_",
            mode="reverse",
            file_mode="w",
        )

        # encode_esa(fasta_path.replace(
        #     ".fasta", ".I2L.rev.fasta"), index_prefix)
        encode_esa(str(fasta_path.with_suffix('.I2L.rev.fasta')), str(index_prefix))

    decoder = DecoderEsa(str(index_prefix))
    return decoder

    
def get_decoder_I2L(fasta_path):
    index_prefix = fasta_path.parent/'.index.tmp'/(fasta_path.stem+'.I2L.esa')/'index'
    if not pathlib.Path(index_prefix).parent.exists():
        index_prefix.parent.mkdir(exist_ok=True, parents=True)
        with fasta.read(str(fasta_path)) as db:
            transferred_data = [
                (item.description, item.sequence.replace("I", "L")) for item in db]
        fasta.write(transferred_data, str(fasta_path.with_suffix(
            ".I2L.fasta")), file_mode="w")

        encode_esa(str(fasta_path.with_suffix('.I2L.fasta')), str(index_prefix))

    decoder = DecoderEsa(str(index_prefix))
    return decoder


def is_in_sync_db(data, syn_db_path):
    syn_db_path = pathlib.Path(syn_db_path)
    syn_decoder = get_decoder_I2L(syn_db_path)

    data["in_syn_db"] = data["Peptide"].apply(
        lambda x: is_xl_seq_in_same_syn_group(x, syn_decoder))
    
    return data


def is_xl_seq_in_real_db(mixed_peptide, decoder):
    line = re.split("\-|\(|\)", mixed_peptide)
    seq1, seq2 = line[0], line[3]

    res1_list = find_protein_in_fasta_db(seq1, decoder)
    res2_list = find_protein_in_fasta_db(seq2, decoder)

    return len(res1_list) > 0 and len(res2_list) > 0


def is_in_real_db(data, real_db_path):
    real_db_path = pathlib.Path(real_db_path)
    real_db_decoder = get_decoder_I2L(real_db_path)

    data["in_real_db"] = data["Peptide"].apply(
        lambda x: is_xl_seq_in_real_db(x, real_db_decoder))
    
    return data
