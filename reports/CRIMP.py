from functools import partial
import tqdm
import pandas as pd
import numpy as np
import re
from pyteomics import fasta
from msms_io.reports.functions_xl import *
import pathlib
from msms_io.fasta.find_protein import find_xl_seq_in_fasta_db
from msms_io.fasta.find_protein import get_decoder_I2L_rev, get_decoder_I2L


def transvert_peptide_modifications(
        peptide1, peptide2, site1, site2):
    mod_dict = {
        'M(35)': 'Oxidation[M]',
        'C(4)': 'Carbamidomethyl[C]',
    }

    # Formalize
    peptide1 = peptide1.replace('I', 'L')
    peptide2 = peptide2.replace('I', 'L')

    mods_list_1 = []
    for idx, aa in enumerate(re.findall('[A-Z][^A-Z]*', peptide1)):
        if aa in mod_dict:
            mods_list_1.append((idx+1, mod_dict[aa]))
    sequence1 = re.sub('\(\d+\)', '', peptide1)

    mods_list_2 = []
    for idx, aa in enumerate(re.findall('[A-Z][^A-Z]*', peptide2)):
        if aa in mod_dict:
            mods_list_2.append((idx+1, mod_dict[aa]))
    sequence2 = re.sub('\(\d+\)', '', peptide2)
    mixed_peptide, mixed_modification = f'{sequence1}({site1})-{sequence2}({site2})', ";".join(
        [f"{x[1]}({x[0]})" for x in mods_list_1] +
        [f"{x[1]}({x[0]+len(mods_list_1)+3})" for x in mods_list_2])

    # Sort order
    return sort_alpha_beta_order(mixed_peptide, mixed_modification)


def transvert_peptide_pp(
        peptide1, peptide2, site1, site2):

    # Formalize
    peptide1 = peptide1.replace('I', 'L')
    peptide2 = peptide2.replace('I', 'L')

    mixed_peptide = f'{peptide1}({site1})-{peptide2}({site2})'

    # Sort order
    return mixed_peptide


def transvert_site(peptide_start_1, coupling_site_1, peptide_start_2, coupling_site_2):

    coupling_site_1 = re.sub('[A-Z]', '', coupling_site_1)
    site_1 = int(coupling_site_1) - int(peptide_start_1) + 1
    coupling_site_2 = re.sub('[A-Z]', '', coupling_site_2)
    site_2 = int(coupling_site_2) - int(peptide_start_2) + 1
    return site_1, site_2


def extract_best_match_site(mixed_line):
    mixed_line = mixed_line.split(';')[0]
    return mixed_line.split('-')


def CRIMP_to_pLink(plink_template_path, inpath, fasta_path):   
    if pathlib.Path(inpath).with_suffix('.pLink.csv').exists():
        return

    df_plink_template = pd.read_csv(plink_template_path, sep=',')

    # CSM, prec
    data = pd.read_csv(inpath)
    data = data[data['Decision'] == 'Accepted']
    data = data[data['Decoy'] == False]
    data = data[(data['Category'] == 'Heterotypic Inter-Protein') | \
                (data['Category'] == 'Homotypic Inter-Protein') | \
                (data['Category'] == 'Intra-Protein')]
    data = data.reset_index(drop=True)
    print('======Total Spec Number:', data.shape[0])

    df_plink = pd.DataFrame(columns=df_plink_template.columns)
    df_plink['Order'] = list(range(1, data.shape[0]+1))
    df_plink['Charge'] = data['z']
    data = data.reset_index().rename(columns={'index': f"Scan"})
    df_plink['Title'] = data[['Run', 'Scan', 'z']].apply(
        lambda x: f'{pathlib.Path(x[0]).stem}.{x[1]}.{x[1]}.{x[2]}.0.dta', axis=1
    )

    data[['_site_1', '_site_2']] = list(
        data[['Peptide 1 Start', 'Best Coupling Site 1', 'Peptide 2 Start', 'Best Coupling Site 2']].swifter.apply(
            lambda x: transvert_site(*x), axis=1))
    
    df_plink[['Peptide', 'Modifications']] = list(
        data[['Peptide 1 Modified', 'Peptide 2 Modified', '_site_1', '_site_2']].swifter.apply(
            lambda x: transvert_peptide_modifications(*x), axis=1))

    df_plink['Peptide_Type'] = 'Cross-Linked'
    df_plink['Protein_Type'] = data['Category'].map({
        'Heterotypic Inter-Protein': 'Inter-Protein',
        'Homotypic Inter-Protein': 'Intra-Protein',
        'Intra-Protein': 'Intra-Protein',
    })

    decoder = get_decoder_I2L_rev(pathlib.Path(fasta_path))
    df_plink['Proteins'] = df_plink['Peptide'].apply(
        lambda x: find_xl_seq_in_fasta_db(x, decoder))
        
    df_plink = df_plink.fillna(0)
    df_plink.to_csv(pathlib.Path(inpath).with_suffix('.pLink.csv'), index=False)
    
    return data


def get_pep_table(path):
    pep_res = pd.read_csv(path)
    pep_res = pep_res[pep_res['Decision'] == 'Accepted']
    pep_res = pep_res[pep_res['Decoy'] == False]

    pep_res['Peptide_Type'] = pep_res['Category'].map({
        'Heterotypic Inter-Protein': 'Cross-Linked',
        'Homotypic Inter-Protein': 'Cross-Linked',
        'Intra-Protein': 'Cross-Linked',
        'Mono-Link': 'Mono-Linked',
        'Loop-Link': 'Loop-Linked',
        'Unbound/Free': 'Regular',
    })
    pep_res = pep_res[pep_res['Peptide_Type'] == 'Cross-Linked']

    pep_res[['Best Coupling Site 1', 'Best Coupling Site 2']] = list(pep_res['Selected Sites'].apply(
        extract_best_match_site))
    pep_res[['_site_1', '_site_2']] = list(
        pep_res[['Peptide 1 Start', 'Best Coupling Site 1', 'Peptide 2 Start', 'Best Coupling Site 2']].swifter.apply(
            lambda x: transvert_site(*x), axis=1))
    pep_res['Peptide'] = list(
        pep_res[['Peptide 1', 'Peptide 2', '_site_1', '_site_2']].apply(
            lambda x: transvert_peptide_pp(*x), axis=1))
    # pep_res['Modifications'] = list(
    #     pep_res[['Peptide 1 Modification', 'Peptide 2 Modification']].fillna('').apply(
    #         lambda x: ';'.join(x), axis=1))
    pep_res['Q-value'] = pep_res['q']

    pep_res = pep_res[[
        'Peptide', 
        # 'Modifications', 
        "Peptide_Type",
        "Score", "Q-value",
    ]]

    return pep_res


def get_rp_table(path):
    rp_res = pd.read_csv(path)
    rp_res = rp_res[rp_res['Decision'] == 'Accepted']
    rp_res = rp_res[rp_res['Decoy'] == False]

    rp_res['Protein_Type'] = rp_res['Category'].map({
        'Heterotypic Inter Protein': 'Inter-Protein',
        'Intra Protein': 'Intra-Protein',
    })

    rp_res['_protein_1'] = rp_res['Protein 1'].apply(
        lambda x: x.split(' ')[0])
    rp_res['_protein_2'] = rp_res['Protein 2'].apply(
        lambda x: x.split(' ')[0])
    rp_res['_site_1'] = rp_res['Site 1'].apply(
        lambda x: int(re.sub('[A-Z]', '', x)))
    rp_res['_site_2'] = rp_res['Site 2'].apply(
        lambda x: int(re.sub('[A-Z]', '', x)))
    rp_res['_rp'] = list(
        rp_res[['_protein_1', '_site_1', '_protein_2', '_site_2']].apply(
            lambda x: f'{x[0]} ({x[1]})-{x[2]} ({x[3]})/', axis=1))
    rp_res['_rp'] = rp_res['_rp'].apply(sort_site_order)
    rp_res['_rp'] = rp_res['_rp'].apply(lambda x: x.strip('/'))
    rp_res['Q-value'] = rp_res['q']

    rp_res = rp_res[[
        '_rp', 'Protein_Type', "Score", "Q-value",
    ]]

    return rp_res


def get_ppi_table(path):
    ppi_res = pd.read_csv(path)
    ppi_res = ppi_res[ppi_res['Decision'] == 'Accepted']
    ppi_res = ppi_res[ppi_res['Decoy'] == False]

    ppi_res['Protein_Type'] = 'Inter-Protein'

    ppi_res['_protein_1'] = ppi_res['Protein 1'].apply(
        lambda x: x.split(' ')[0])
    ppi_res['_protein_2'] = ppi_res['Protein 2'].apply(
        lambda x: x.split(' ')[0])
    ppi_res['_ppi'] = list(
        ppi_res[['_protein_1', '_protein_2']].apply(
            lambda x: f'{x[0]} -{x[1]} ', axis=1))
    ppi_res['_ppi'] = ppi_res['_ppi'].apply(sort_protein_order)
    ppi_res['_ppi'] = ppi_res['_ppi'].apply(lambda x: x.strip('/'))
    ppi_res['Q-value'] = ppi_res['q']

    ppi_res = ppi_res[[
        '_ppi', 'Protein_Type', "Score", "Q-value",
    ]]

    return ppi_res