import pandas as pd
from ..load_results import get_PSM_from_pF
from .xl_utils import rmv_xl_site, split_xl_protein


def get_protein_set_from_list(proteins_list):
    '''Turn pandas serier to python set'''
    proteins_list = list(map(lambda x: x.strip("/").split("/"), proteins_list))
    protein_set = set()
    for proteins_line in proteins_list:
        for protein in proteins_line:
            protein = protein.strip()
            if protein not in protein_set:
                protein_set.add(protein)
    return protein_set


def get_protein_set_from_pFind(pFind_res_path):
    '''Load and extract protein set from pFind'''
    PSM_res = get_PSM_from_pF(pFind_res_path)
    protein_set = get_protein_set_from_list(list(set(
        PSM_res["proteins"].tolist())))
    return protein_set


def get_xl_protein_set_from_pLink(pLink_res_prefix):
    '''Load and extract xl protein set from pLink'''
    xl_results = pd.read_csv(f"{pLink_res_prefix}filtered_cross-linked_spectra.csv")

    xl_protein_set = get_protein_set_from_list(list(set(
            xl_results["Proteins"].apply(rmv_xl_site).apply(split_xl_protein).tolist())))
    return xl_protein_set


def get_protein_set_from_pLink(pLink_res_prefix):
    '''Load and extract protein set from pLink'''
    xl_results = pd.read_csv(f"{pLink_res_prefix}filtered_cross-linked_spectra.csv")
    mono_results = pd.read_csv(f"{pLink_res_prefix}filtered_mono-linked_spectra.csv")
    loop_results = pd.read_csv(f"{pLink_res_prefix}filtered_loop-linked_spectra.csv")
    regular_results = pd.read_csv(f"{pLink_res_prefix}filtered_regular_spectra.csv")

    xl_protein_set = get_protein_set_from_list(list(set(
            xl_results["Proteins"].apply(rmv_xl_site).apply(split_xl_protein).tolist())))
    mono_protein_set = get_protein_set_from_list(list(set(
            mono_results["Proteins"].apply(rmv_xl_site).tolist())))
    loop_protein_set = get_protein_set_from_list(list(set(
            loop_results["Proteins"].apply(rmv_xl_site).tolist())))
    regular_protein_set = get_protein_set_from_list(list(set(
            regular_results["Proteins"].tolist())))
    
    protein_set = set()
    protein_set.update(xl_protein_set)
    protein_set.update(mono_protein_set)
    protein_set.update(loop_protein_set)
    protein_set.update(regular_protein_set)
    return protein_set