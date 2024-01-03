import pathlib
import pandas as pd
from .functions_xl import *
from .functions import *
from pyteomics import fasta
from msms_io.fasta.fasta_toolbox import encode_esa, DecoderEsa
from msms_io.fasta.find_protein import find_xl_seq_in_fasta_db, is_xl_seq_in_same_syn_group
import os
import sys
from msms_io.fasta.find_protein import get_decoder_I2L_rev
import re


MAP_TRANS_MODIFICAIONT = {'B': 'Carbamidomethyl[C]', 'm': 'Oxidation[M]'}
MAP_TRANS_AA = {'B': 'C', 'm': 'M'}
AA_LIST = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
DELIMITER = ';'

def keep_xl_df(df):
    df = df[df['Peptide2'] != '0']
    df = df[df['Peptide2'] != '1']
    return df


def get_tt_df(df):
    def is_tt(ser):
        if ser['Protein 1'][:4] != 'DEC_' and ser['Protein 2'][:4] != 'DEC_':
            return 1
        else:
            return 0

    df['is_tt'] = df.apply(is_tt, axis=1)
    df = df[df['is_tt'] == 1]
    return df


def is_pparse_title(s):
    segs = s.split('.')

    if len(segs) >= 6:
        scan1 = segs[-5]
        scan2 = segs[-4]
        charge = segs[-3]
        pparseid = segs[-2]
        suffix = segs[-1]

        try:
            int(scan1)
            int(scan2)
            int(charge)
            int(pparseid)
            if suffix == 'dta':
                return True
        except ValueError:
            pass
        return False
    return False


def is_pxtract_title(s):
    segs = s.split('.')

    if len(segs) >= 5:
        scan1 = segs[-4]
        scan2 = segs[-3]
        charge = segs[-2]
        suffix = segs[-1]

        try:
            int(scan1)
            int(scan2)
            int(charge)
            if suffix == 'dta':
                return True
        except ValueError:
            pass
        return False
    return False


def merox2plink_title(ser):

    # pParse导出
    if is_pparse_title(ser['Scan number']):
        dta_name = ser['Scan number'][:-4]
        dta_name = '.'.join(dta_name.split('.')[:-1])

    # pXtract导出
    elif is_pxtract_title(ser['Scan number']):
        dta_name = ser['Scan number'][:-4]
    
    # MSConvert导出
    elif 'File:"' in ser['Scan number'][:-4]:
        dta_name = ser['Scan number'].split()[0]

    pparse_id = '0'
    extended_name = 'dta'
    return '.'.join([dta_name, pparse_id, extended_name])



def merox2plink_proteintype(ser):
    """获得intra/inter"""

    prt1_str = ser['Protein 1']
    prt2_str = ser['Protein 2']

    def get_prt_ls(prt_str):
        prt_ls = []
        if '(>' not in prt_str:
            prt_ls.append(prt_str.strip()[1:])
        elif '/>' not in prt_str:
            segs = prt_str.strip().split('(>')
            prt_ls.append(segs[0][1:])
            prt_ls.append(segs[1][:-1])
            # print(prt_str, '\n', prt_ls)
        else:
            segs = prt_str.strip().split('(>')
            prt_ls.append(segs[0][1:])
            prt_ls.extend(segs[1][:-1].split('/>'))
            # print(prt_str, '\n', prt_ls)
        return prt_ls

    prt1_ls = get_prt_ls(prt1_str)
    prt2_ls = get_prt_ls(prt2_str)

    if set(prt1_ls).intersection(set(prt2_ls)):
        return 'Intra-Protein'
    else:
        return 'Inter-Protein'
        

def merox2plink_xl_pep_mod(ser):

    pep1 = ser['Peptide 1'][1:-1]
    pep2 = ser['Peptide2'][1:-1]
    
    link1 = int(ser['best linkage position peptide 1'][1:]) # 从1开始，N端为{0
    link2 = int(ser['best linkage position peptide 2'][1:])
    if link1 == 0:
        link1 = 1
    if link2 == 0:
        link2 = 1
     
    site2mod = {}

    seq = []
    for i, aa in enumerate(pep1):
        if aa not in AA_LIST:
            seq.append(MAP_TRANS_AA[aa])
            site2mod[i+1] = MAP_TRANS_MODIFICAIONT[aa]
        else:
            seq.append(aa)
    seq1 = ''.join(seq)

    seq = []
    for i, aa in enumerate(pep2):
        if aa not in AA_LIST:
            seq.append(MAP_TRANS_AA[aa])
            site2mod[i+1+len(seq1)+3] = MAP_TRANS_MODIFICAIONT[aa]
        else:
            seq.append(aa)
    seq2 = ''.join(seq)

    mods_ls = []
    if site2mod:
        for key, value in site2mod.items():
            mods_ls.append(value + '('+ str(key) + ')')
        mod_str = ';'.join(mods_ls)
    else:
        mod_str = 'null'

    seq_str = seq1 + '(' + str(link1) + ')-' + seq2 + '(' + str(link2) + ')'

    return [seq_str, mod_str]


def MeroX_to_pLink(plink_template_path, inpath, fasta_path, linker_name='null'):
    if pathlib.Path(inpath).with_suffix('.pLink.csv').exists():
        return

    df_plink_template = pd.read_csv(plink_template_path, sep=',')
    print('===@()', sys._getframe().f_code.co_name)

    df_merox = pd.read_csv(inpath, sep=DELIMITER)
    # print(df_merox.iloc[0])
    print('======Total Spec Number:', df_merox.shape[0])
    df_merox = keep_xl_df(df_merox)
    print('======XL Spec Number:', df_merox.shape[0])
    df_merox = get_tt_df(df_merox)
    print('======TT Spec Number:', df_merox.shape[0])
    df_merox = df_merox.reset_index(drop=True)

    psm_num = df_merox.shape[0]
    df_plink = pd.DataFrame(columns=df_plink_template.columns)

    df_plink['Order'] = list(range(1, psm_num+1))
    df_plink['Title'] = df_merox.apply(merox2plink_title, axis=1)
    # Should be uCSM
    df_plink["_scan_id"] = df_plink["Title"].swifter.apply(
        lambda x: re.search("(?<=^).*?\.\d+\.\d+(?=\.\d+\.\d+\.dta\;?)", x).group())
    assert df_plink['_scan_id'].duplicated().sum() == 0, 'rCSM rather than uCSM!'
    df_plink.drop(columns=['_scan_id'], inplace=True)
    df_plink['Charge'] = df_merox['Charge']
    df_plink[['Peptide', 'Modifications']] = list(df_merox.apply(lambda x: merox2plink_xl_pep_mod(x), axis=1))
    df_plink['Peptide'] = df_plink['Peptide'].map(lambda x: x.replace("I", "L"))
    df_plink['Peptide_Type'] = 'Cross-Linked'
    df_plink['Linker'] = linker_name
    decoder = get_decoder_I2L_rev(pathlib.Path(fasta_path))
    df_plink['Proteins'] = df_plink['Peptide'].apply(
        lambda x: find_xl_seq_in_fasta_db(x, decoder))
    df_plink['Protein_Type'] = list(df_merox.apply(merox2plink_proteintype, axis=1))

    df_plink['Precursor_Mass'] = df_merox['M+H+'] # pLink中为不带电质量+(H+)
    df_plink['Peptide_Mass'] = df_merox['Calculated Mass'] # pLink中为不带电质量+(H+)

    df_plink = df_plink.fillna(0)

    df_plink.to_csv(pathlib.Path(inpath).with_suffix('.pLink.csv'), index=False)
