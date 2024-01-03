from functools import partial
import tqdm
import pandas as pd
import numpy as np
import re
from pyteomics import fasta
from msms_io.reports.functions_xl import *
import pathlib
from msms_io.fasta.find_protein import find_xl_seq_in_fasta_db, is_xl_seq_in_same_syn_group
from msms_io.fasta.find_protein import get_decoder_I2L_rev, get_decoder_I2L

transfer_dict = {
    'Carbamidomethyl': 'Carbamidomethyl',
    'Oxidation of Methionine': 'Oxidation',
}


def format_Scout_modification(modification, transfer_dict, offset=0):
    if modification == '':
        return ''

    modification_list = modification.split(';')
    res_list = []
    for mod in modification_list:
        aa = mod[0]
        pos = int(re.search('\d+', mod).group()) + offset
        descr = re.search('(?<=\().*(?=\))', mod).group()
        descr = transfer_dict[descr]
        res_list.append(f'{descr}[{aa}]({pos})')
    return ';'.join(res_list)


def Scout_to_pLink(plink_template_path, inpath, fasta_path):   
    if pathlib.Path(inpath).with_suffix('.pLink.csv').exists():
        return

    df_plink_template = pd.read_csv(plink_template_path, sep=',')

    # CSM, prec
    data = pd.read_csv(inpath)
    data = data[data['Link-Type'] != 'Looplink']
    data = data.reset_index(drop=True)
    print('======Total Spec Number:', data.shape[0])

    df_plink = pd.DataFrame(columns=df_plink_template.columns)
    df_plink['Order'] = list(range(1, data.shape[0]+1))
    df_plink['Title'] = data[['File', 'Scan', 'Precursor charge']].apply(
        lambda x: f'{pathlib.Path(x[0]).stem}.{x[1]}.{x[1]}.{x[2]}.0.dta', axis=1
    )
    df_plink['Charge'] = data['Precursor charge']

    # Assert uCSM, which means no mixed spectra
    assert sum(data[['File', 'Scan']].duplicated()) == 0

    df_plink['Peptide'] = data[['Alpha peptide', 'Alpha peptide position', 'Beta peptide', 'Beta peptide position']].apply(
        lambda x: f'{x[0]}({x[1]})-{x[2]}({x[3]})', axis=1
    )
    df_plink['Peptide'] = df_plink['Peptide'].map(lambda x: x.replace("I", "L"))

    # C8 (Carbamidomethyl); M5 (Oxidation of Methionine) to: Carbamidomethyl[C](8);Oxidation[M](5)
    data['Modifications_1'] = data['Alpha modification(s)'].fillna('').map(
        lambda x: format_Scout_modification(x, transfer_dict)
    )
    data['Modifications_2'] = data[['Beta modification(s)', 'Alpha peptide']].fillna('').apply(
        lambda x: format_Scout_modification(x[0], transfer_dict, len(x[1]) + 3), axis=1
    )
    df_plink['Modifications'] = data[['Modifications_1', 'Modifications_2']].apply(
        lambda x: (x[0] + ';' + x[1]).strip(';'), axis=1
    )

    df_plink['Peptide_Type'] = 'Cross-Linked'
    df_plink['Protein_Type'] = data['Link-Type'].map({
        'Heteromeric-inter': 'Inter-Protein',
        'Homomeric-inter': 'Intra-Protein',
        'Intralink': 'Intra-Protein',
    })

    decoder = get_decoder_I2L_rev(pathlib.Path(fasta_path))
    df_plink['Proteins'] = df_plink['Peptide'].apply(
        lambda x: find_xl_seq_in_fasta_db(x, decoder))
        
    df_plink = df_plink.fillna(0)
    df_plink.to_csv(pathlib.Path(inpath).with_suffix('.pLink.csv'), index=False)
    
    return data


def load_seq_Scout(project_dir, fasta_path):    
    # Fields  need for evaluation:
    # Peptide, 
    # Proteins, Protein_Type

    project_dir = pathlib.Path(project_dir)
    fasta_path = pathlib.Path(fasta_path)

    if (project_dir/'rp_preprocessed.csv').exists():
        return pd.read_csv(project_dir/'rp_preprocessed.csv')

    data = pd.read_csv(project_dir / 'rp.csv').fillna('')

    # FAAYA [K] AYPQEAAEFTR -> FAAYAKAYPQEAAEFTR
    data['Peptide_1'] = data['Alpha peptide'].apply(
        lambda x: re.sub(' |\[|\]', '', x)
    )
    data['Peptide_2'] = data['Beta peptide'].apply(
        lambda x: re.sub(' |\[|\]', '', x)
    )

    data['Peptide'] = data[['Peptide_1', 'Alpha peptide position', 'Peptide_2', 'Beta peptide position']].apply(
        lambda x: f'{x[0]}({x[1]})-{x[2]}({x[3]})', axis=1
    )
    data['Peptide'] = data['Peptide'].map(lambda x: x.replace("I", "L"))
    data['Modifications'] = ''

    data[['Peptide', 'Modifications']] = list(data[['Peptide', 'Modifications']].apply(
        lambda x: sort_alpha_beta_order(*x), axis=1))

    data = data.drop(columns='Modifications')

    data['Peptide_Type'] = 'Cross-Linked'

    sequence_decoder = get_decoder_I2L_rev(fasta_path)
    data["Proteins"] = data["Peptide"].apply(
        lambda x: find_xl_seq_in_fasta_db(x, sequence_decoder))
    data["Proteins"] = data["Proteins"].map(
        sort_site_order)
        
    data['Protein_Type'] = data['Link-Type'].map({
        'Heteromeric-inter': 'Inter-Protein',
        'Homomeric-inter': 'Intra-Protein',
        'Intralink': 'Intra-Protein',
    })

    data.to_csv(project_dir/'rp_preprocessed.csv', index=False)
    return data


def load_ppi_Scout(project_dir, fasta_path):    
    # Fields  need for evaluation:
    # _ppi, Protein_Type

    project_dir = pathlib.Path(project_dir)
    fasta_path = pathlib.Path(fasta_path)

    if (project_dir/'ppi_preprocessed.csv').exists():
        return pd.read_csv(project_dir/'ppi_preprocessed.csv')

    data = pd.read_csv(project_dir / 'ppi.csv').fillna('')
    data = data[data['PPI'].apply(lambda x: len(re.findall(';', x)) == 0)] # rmv non unique

    # P0A7W1 <=> P27302 -> sp|P0A7W1|RS5_ECOLI - sp|P27302|TKT1_ECOLI /
    with fasta.read(str(fasta_path.with_suffix('.I2L.fasta'))) as db:
        descriptions = [
            (item.description.split(' ')[0]) for item in db]
    search_dict = {descrption.split('|')[1]: descrption for descrption in descriptions}

    def process(Scout_ppi):
        ppi_list = Scout_ppi.split(' <=> ')
        ppi_1 = search_dict[ppi_list[0]]
        ppi_2 = search_dict[ppi_list[1]]
        return f'{ppi_1} -{ppi_2} '

    data['_ppi'] = data['PPI'].apply(process)

    data['Protein_Type'] = data['Link-Type'].map({
        'Inter': 'Inter-Protein',
    })

    data.to_csv(project_dir/'ppi_preprocessed.csv', index=False)
    return data
