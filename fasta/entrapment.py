from pyteomics.pyteomics import fasta, parser    
import pathlib
import ujson
import tqdm
import pandas as pd
import itertools
import re
from msms_io.tools.df_split_merge import split_column_into_columns, group_by_col


def append_entrapment_to_database(
    original_fasta_path,
    entrapment_fasta_path,
    res_fasta_path,
    log_path,
):
    log_file = open(log_path, 'w')
    with fasta.read(original_fasta_path) as db:
        transferred_data = [
            (item.description, item.sequence.replace("I", "L")) for item in db]
    fasta.write(transferred_data, original_fasta_path.replace(
        ".fasta", ".I2L.fasta"), file_mode="w")
    fasta.write_decoy_db(original_fasta_path.replace(".fasta", ".I2L.fasta"),
                            original_fasta_path.replace(
        ".fasta", ".I2L.rev.fasta"),
        prefix="REV_",
        mode="reverse",
        file_mode="w",
    )
    log_file.write('Length of original fasta: ')
    log_file.write(str(len(transferred_data)))
    log_file.write('\n')
    
    with fasta.read(entrapment_fasta_path) as db:
        transferred_data = [
            (item.description, item.sequence.replace("I", "L")) for item in db]
    fasta.write(transferred_data, entrapment_fasta_path.replace(
        ".fasta", ".I2L.fasta"), file_mode="w")
    fasta.write_decoy_db(entrapment_fasta_path.replace(".fasta", ".I2L.fasta"),
                            entrapment_fasta_path.replace(
        ".fasta", ".I2L.rev.fasta"),
        prefix="REV_",
        mode="reverse",
        file_mode="w",
    )
    log_file.write('Length of entrapment fasta: ')
    log_file.write(str(len(transferred_data)))
    log_file.write('\n')

    with open(res_fasta_path, 'w') as output:
        fasta.write(fasta.read(original_fasta_path), output)        
        fasta.write(fasta.read(entrapment_fasta_path), output)

    with fasta.read(res_fasta_path) as db:
        transferred_data = [
            (item.description, item.sequence.replace("I", "L")) for item in db]
    fasta.write(transferred_data, res_fasta_path.replace(
        ".fasta", ".I2L.fasta"), file_mode="w")
    fasta.write_decoy_db(res_fasta_path.replace(".fasta", ".I2L.fasta"),
                            res_fasta_path.replace(
        ".fasta", ".I2L.rev.fasta"),
        prefix="REV_",
        mode="reverse",
        file_mode="w",
    )
    log_file.write('Length of concated fasta: ')
    log_file.write(str(len(transferred_data)))
    log_file.write('\n')

    log_file.close()


def in_silico_digestion_common(
    fasta_path,
    output_path,
):
    if pathlib.Path(output_path).exists():       
        peptide_df = pd.read_csv(output_path)
        return peptide_df

    res_peptide_list = list()
    for record in tqdm.tqdm(fasta.read(fasta_path)):
        peptides = [[*x, record.description.split(' ')[0]] for x in parser.icleave_trypsin_P_common(record.sequence.replace('I', 'L'))]
        res_peptide_list.extend(peptides)
    
    peptide_df = pd.DataFrame(res_peptide_list, columns=[
        'start', 'sequence', 'is_N_term', 'is_C_term', 'missed_cleavages', 'protein'])
    peptide_df.drop(columns=['start'], inplace=True)

    peptide_df = peptide_df.groupby('sequence').agg({
        'sequence': 'first',
        'missed_cleavages': 'first',
        'is_N_term': 'max',
        'is_C_term': 'max',
        'protein': lambda x: ';'.join(x),
        })

    peptide_df.to_csv(output_path, index=False)
    
    return peptide_df



def add_modification(sequence, src_aa, dest_aa, max=3):
  result = []

  M_indices = [i for i in range(len(sequence)) if sequence[i] == src_aa]
  for r in range(len(M_indices) + 1):
    comb = itertools.combinations(M_indices, r)
    for indices in comb:
      if len(indices) > max:
        continue
      new_s = "".join([sequence[i] if i not in indices else dest_aa for i in range(len(sequence))])
      result.append(new_s)

  return result


def in_silico_digestion_add_modification(
    common_path,
    output_path,
    src_aa, dest_aa,
):
    if pathlib.Path(output_path).exists():       
        peptide_df = pd.read_csv(output_path)
        return peptide_df

    peptide_df = pd.read_csv(common_path)
    peptide_df['sequence'] = peptide_df['sequence'].apply(lambda x: ';'.join(add_modification(x, src_aa, dest_aa)))
    # peptide_df = split_column_into_columns(peptide_df, 'sequence', ';')
    
    peptide_df.to_csv(output_path, index=False)
    
    return peptide_df


def add_mono_linker(sequence, is_C_term, src_aa, dest_aa):
  result = []

  M_indices = [i for i in range(len(sequence)) if sequence[i] == src_aa]
  comb = itertools.combinations(M_indices, 1)
  for indices in comb:
      # Pep C term but not Pro C term
      if indices[0] == len(sequence) - 1 and not is_C_term:
        continue

      new_s = "".join([sequence[i] if i not in indices else dest_aa for i in range(len(sequence))])

      result.append(new_s)

  return result


def add_loop_linker(sequence, is_C_term, src_aa, dest_aa):
  result = []

  M_indices = [i for i in range(len(sequence)) if sequence[i] == src_aa]
  comb = itertools.combinations(M_indices, 2)
  for indices in comb:
      # Pep C term but not Pro C term
      if (indices[0] == len(sequence) - 1 or indices[1] == len(sequence) - 1) and not is_C_term:
        continue

      new_s = "".join([sequence[i] if i not in indices else dest_aa for i in range(len(sequence))])

      result.append(new_s)

  return result


def in_silico_digestion_add_mono(
    common_path,
    output_path,
    src_aa, dest_aa,
):
    if pathlib.Path(output_path).exists():       
        peptide_df = pd.read_csv(output_path)
        return peptide_df

    peptide_df = pd.read_csv(common_path)
    peptide_df = peptide_df[peptide_df['sequence'].str.contains(src_aa)]
    peptide_df['sequence'] = peptide_df[['sequence', 'is_C_term']].apply(lambda x: ';'.join(add_mono_linker(x[0], x[1], src_aa, dest_aa)), axis=1)
    peptide_df = peptide_df[peptide_df['sequence'] != '']
    peptide_df = split_column_into_columns(peptide_df, 'sequence', ';')
        
    peptide_df.to_csv(output_path, index=False)
    
    return peptide_df


def in_silico_digestion_add_loop(
    common_path,
    output_path,
    src_aa, dest_aa,
):
    if pathlib.Path(output_path).exists():       
        peptide_df = pd.read_csv(output_path)
        return peptide_df

    peptide_df = pd.read_csv(common_path)
    peptide_df = peptide_df[peptide_df['sequence'].apply(lambda x: len(re.findall(src_aa, x)) >= 2)]
    peptide_df['sequence'] = peptide_df[['sequence', 'is_C_term']].apply(lambda x: ';'.join(add_loop_linker(x[0], x[1], src_aa, dest_aa)), axis=1)
    peptide_df = peptide_df[peptide_df['sequence'] != '']
    peptide_df = split_column_into_columns(peptide_df, 'sequence', ';')
        
    peptide_df.to_csv(output_path, index=False)
    
    return peptide_df


def in_silico_digestion_asem_xl_homo(
    mono_path,
    output_path,      
):
    if pathlib.Path(output_path).exists():       
        peptide_df = pd.read_csv(output_path)
        return peptide_df

    peptide_df = pd.read_csv(mono_path)
    peptide_df.drop(columns=['is_N_term', 'is_C_term', 'missed_cleavages'], inplace=True)
    
    peptide_df['sequence_a'] = peptide_df['sequence'].copy()
    peptide_df['sequence_b'] = peptide_df['sequence'].copy()
    peptide_df.drop('sequence', axis=1, inplace=True)
    peptide_df['Peptide_Type'] = 'Homo'

    peptide_df.to_csv(output_path, index=False)


def in_silico_digestion_asem_xl(
    mono_path,
    output_path,      
):
    if pathlib.Path(output_path).exists():       
        peptide_df = pd.read_csv(output_path)
        return peptide_df

    peptide_df = pd.read_csv(mono_path)
    peptide_df.drop(columns=['is_N_term', 'is_C_term', 'missed_cleavages'], inplace=True)
    
    peptide_df = peptide_df.merge(peptide_df, how='cross', suffixes=('_a', '_b'))
    peptide_df
    # TODO: how to deal with xl search space?
   