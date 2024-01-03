from pyteomics.fasta import FASTA, write
import pandas as pd


def build_fasta_from_protein_table(protein_table_path, column, fasta_original_path, fasta_dest_path):
    protein_trace_table = pd.read_csv(protein_table_path)
    protein_set = set(
        protein_trace_table[pd.notna(protein_trace_table[column])]["_protein"].to_list())

    original_fasta = FASTA(fasta_original_path)
    original_fasta_len = 0
    dest_fasta = list()
    for descr, seq in original_fasta:
        original_fasta_len += 1
        # Only keep target database
        if descr.split(" ")[0] in protein_set:
            dest_fasta.append((descr, seq))
    print("Proteins num in original fasta:", original_fasta_len)
    print("Proteins num in res:", len(protein_set))
    print("Proteins num in dest fasta:", len(dest_fasta))

    write(dest_fasta, fasta_dest_path, file_mode="w")


def build_fasta_from_protein_set(protein_set, fasta_original_path, fasta_dest_path):

    original_fasta = FASTA(fasta_original_path)
    original_fasta_len = 0
    dest_fasta = list()
    for descr, seq in original_fasta:
        original_fasta_len += 1
        # Only keep target database
        if descr.split(" ")[0] in protein_set:
            dest_fasta.append((descr, seq))
    print("Proteins num in original fasta:", original_fasta_len)
    print("Proteins num in res:", len(protein_set))
    print("Proteins num in dest fasta:", len(dest_fasta))

    write(dest_fasta, fasta_dest_path, file_mode="w")
