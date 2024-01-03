# -*- coding: UTF-8 -*-
# @Date     : 1st Arp, 2020
# @Author   : Northblue

"""
XlinkX交联鉴定结果文件转换为pLink格式
"""

"""20210426标明intra/inter"""

import sys
import pandas as pd
from msms_io.fasta.find_protein import get_decoder_I2L_rev
from msms_io.fasta.find_protein import find_xl_seq_in_fasta_db
from msms_io.reports.sort_order import sort_modification
import pathlib
import re


MAP_TRANS_MODIFICAIONT = {'Carbamidomethyl': 'Carbamidomethyl[C]', 'Oxidation': 'Oxidation[M]'}
DELIMITER = '\t'


def xlinkx2plink_title(ser):
    raw_name = ser['Spectrum File'][:-4]
    scan_id = str(ser['All Scans'])
    charge = str(ser['Charge'])
    pparse_id = '0'
    extended_name = 'dta'
    return '.'.join([raw_name, scan_id, scan_id, charge, pparse_id, extended_name])

def xlinkx2plink_xl_pep_mod(ser):
    seq1 = ser['Sequence A']
    seq2 = ser['Sequence B']

    link1 = int(ser['Crosslinker Position A'])
    link2 = int(ser['Crosslinker Position B'])
    if link1 == 0:
        link1 = 1
    if link2 == 0:
        link2 = 1

    linker = ser['Crosslinker']
    site2mod = {}

    mods1 = ser['Modifications A']
    mods2 = ser['Modifications B']

    for mod in mods1.split(';'):
        mod_name = mod.split('(')[1][:-1]
        if mod_name != linker:
            mod_site = int(mod.split('(')[0][1:])
            site2mod[mod_site] = MAP_TRANS_MODIFICAIONT[mod_name]

    for mod in mods2.split(';'):
        mod_name = mod.split('(')[1][:-1]
        if mod_name != linker:
            mod_site = int(mod.split('(')[0][1:])
            site2mod[mod_site+3+len(seq1)] = MAP_TRANS_MODIFICAIONT[mod_name]

    mods_ls = []
    if site2mod:
        for key, value in site2mod.items():
            mods_ls.append(value + '('+ str(key) + ')')
        mod_str = ';'.join(mods_ls)
    else:
        mod_str = 'null'

    seq_str = seq1 + '(' + str(link1) + ')-' + seq2 + '(' + str(link2) + ')'

    return [seq_str, mod_str]


def XlinkX_to_pLink(plink_template_path, inpath, fasta_path):
    if pathlib.Path(inpath).with_suffix('.pLink.csv').exists():
        return

    df_plink_template = pd.read_csv(plink_template_path, sep=',')
    print('===@()', sys._getframe().f_code.co_name)

    df_xlinkx = pd.read_csv(inpath, sep=DELIMITER)
    # print(df_xlinkx.iloc[0])
    print('======Total Spec Number:', df_xlinkx.shape[0])

    psm_num = df_xlinkx.shape[0]
    df_plink = pd.DataFrame(columns=df_plink_template.columns)

    df_plink['Order'] = list(range(1, psm_num+1))
    df_plink['Title'] = list(df_xlinkx.apply(xlinkx2plink_title, axis=1))
    # Should be uCSM
    df_plink["_scan_id"] = df_plink["Title"].swifter.apply(
        lambda x: re.search("(?<=^).*?\.\d+\.\d+(?=\.\d+\.\d+\.dta\;?)", x).group())
    assert df_plink['_scan_id'].duplicated().sum() == 0, 'rCSM rather than uCSM!'
    df_plink['Charge'] = df_xlinkx['Charge']
    
    df_plink[['Peptide', 'Modifications']] = list(df_xlinkx.apply(lambda x:xlinkx2plink_xl_pep_mod(x), axis=1))
    df_plink['Modifications'] = df_plink['Modifications'].map(sort_modification)
    df_plink['Peptide'] = df_plink['Peptide'].map(lambda x: x.replace("I", "L"))
    
    df_plink['Peptide_Type'] = 'Cross-Linked'
    df_plink['Linker'] = df_xlinkx['Crosslinker']
    df_plink['Protein_Type'] = df_xlinkx['Crosslink Type'].map({'Intra':'Intra-Protein', 'Inter':'Inter-Protein'})

    decoder = get_decoder_I2L_rev(pathlib.Path(fasta_path))
    df_plink['Proteins'] = df_plink['Peptide'].apply(
        lambda x: find_xl_seq_in_fasta_db(x, decoder))

    df_plink = df_plink.fillna(0)
    df_plink.to_csv(pathlib.Path(inpath).with_suffix('.pLink.csv'), index=False)
