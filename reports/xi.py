# -*- coding: UTF-8 -*-
# @Date     : 1st Arp, 2020
# @Author   : Northblue

"""将xi的结果转换为pLink2结果格式
    1. 谱图层次用于定量和对比
    2. 肽段层次用于对比
"""
"""20210426标明intra/inter"""

import sys
import pandas as pd
from msms_io.fasta.find_protein import get_decoder_I2L_rev
from msms_io.fasta.find_protein import find_xl_seq_in_fasta_db
import pathlib

MAP_TRANS_MODIFICAIONT = {
    'cm': 'Carbamidomethyl[C]',
    'ox': 'Oxidation[M]',
    'leikermono': 'Leiker_Mono[K]',
    'dssmono': 'DSS_Mono[K]',
    'dssooh': 'DSSO_H2O[K]',
    'dssoh2': 'DSSO_NH2[K]',
    'dsbuoh': 'DSBU_H2O[K]',
    'dsbuh2': 'DSBU_NH2[K]',
    'dsbunh2': 'DSBU_NH2[K]',
    'dsbuloop': 'DSBU_LOOP[K]',
    'bs3oh':'BS3_OH[K]',
    'bs3nh2':'BS3_NH2[K]'
    }
# MAP_TRANS_MODIFICAIONT = {'cm': 'Carbamidomethyl[C]', 'ox': 'Oxidation[M]', 'bs3oh':'BS3_OH', 'bs3nh2':'BS3_NH2', 'leikeroh':'Leiker_OH', 'leikernh2':'Leiker_NH2'}
# MAP_TRANS_MODIFICAIONT = {'cm': 'Carbamidomethyl[C]', 'ox': 'Oxidation[M]', 'dssooh': 'DSSO_H2O[K]'}


MASS_H = 1.00782503214
MASS_PROTON = 1.007276466621

# 带mono修饰的直接当做单肽处理

# delimiter = ',' # 
delimiter = '\t'
print('Xi Result分隔符', '"%s"'%(delimiter))

def xi2plink_title(ser):
    """xi中提取plink的谱图Title"""

    mgf_name = ser['PeakListFileName'][:-4]
    mgf_name = mgf_name.split('.')[0] # 去除.HCD.FTMS.mgf

    if pd.isnull(ser['scan']):
        scan_id = str(0)
    else:
        scan_id = str(ser['scan'])
    # print(scan_id)

    exp_charge = str(ser['exp charge'])
    pparse_id = '0'
    extended_name = 'dta'

    return '.'.join([mgf_name, scan_id, scan_id, exp_charge, pparse_id, extended_name])

def xi2plink_proteintype(ser):
    """获得intra/inter"""

    fdr_group = ser['fdrGroup'].strip()
    fdr_group2protein_type = {'Self':'Intra-Protein', 'Between':'Inter-Protein'}

    return fdr_group2protein_type[fdr_group]


def xi2plink_linear_seq_mod(ser):
    """xi中提取单肽的序列和修饰"""

    pep = ser['PepSeq1']

    seq_ls = []
    mod_ls = []
    map_mod = {}
    site = 0
    for aa in pep:
        if aa.isupper():
            if mod_ls: # 修饰列表不为空
                map_mod[site] = MAP_TRANS_MODIFICAIONT[''.join(mod_ls)]
                mod_ls.clear()
            seq_ls.append(aa)
            site += 1
        else:
            mod_ls.append(aa)

    seqs = ''.join(seq_ls)

    mods_ls = []
    if map_mod:
        for key, value in map_mod.items():
            mods_ls.append(value + '('+ str(key) + ')')
        mod_str = ';'.join(mods_ls)
    else:
        mod_str = 'null'

    return [seqs, mod_str]


def xi2plink_crosslink_seq_mod(ser):
    """xi中提取交联的序列和修饰"""

    pep1 = ser['PepSeq1']
    pep2 = ser['PepSeq2']
    link_pos1 = ser['LinkPos1'] #从1开始, N端交联位点在1
    link_pos2 = ser['LinkPos2']

    map_mod = {}

    seq_ls = []
    mod_ls = []
    site = 0
    for aa in pep1:
        if aa.isupper():
            if mod_ls: # 修饰列表不为空
                map_mod[site] = MAP_TRANS_MODIFICAIONT[''.join(mod_ls)]
                mod_ls.clear()
            seq_ls.append(aa)
            site += 1
        else:
            mod_ls.append(aa)
    
    # 修饰在最后一位
    if mod_ls: # 修饰列表不为空
        map_mod[site] = MAP_TRANS_MODIFICAIONT[''.join(mod_ls)]

    seqs1 = ''.join(seq_ls)
    seqs1_len = len(seqs1)

    seq_ls = []
    mod_ls = []
    site = 0
    for aa in pep2:
        if aa.isupper():
            if mod_ls: # 修饰列表不为空
                map_mod[site + seqs1_len + 3] = MAP_TRANS_MODIFICAIONT[''.join(mod_ls)]
                mod_ls.clear()
            seq_ls.append(aa)
            site += 1
        else:
            mod_ls.append(aa)
    seqs2 = ''.join(seq_ls)

    # 修饰在最后一位
    if mod_ls: # 修饰列表不为空
        map_mod[site + seqs1_len + 3] = MAP_TRANS_MODIFICAIONT[''.join(mod_ls)]

    mods_ls = []
    if map_mod:
        for key, value in map_mod.items():
            mods_ls.append(value + '('+ str(key) + ')')
        mod_str = ';'.join(mods_ls)
    else:
        mod_str = 'null'

    # if map_mod:
    #     if ser['scan'] == 12066:
    #         print(pep1, pep2)
    #         print(mod_str)

    seqs_str = seqs1 + '(' + str(link_pos1) + ')-' + seqs2 + '(' + str(link_pos2) + ')'

    return [seqs_str, mod_str]

# def xi2plink_linear_seq_mod(ser):
    # """xi中提取单肽的蛋白AC"""


def set_scan_order(df):
    """xi如果没有scan号, 按照排序赋予scan"""


    def set_title(x):
        raw_name = '.'.join(x['Title'].split('.')[:-5])
        scan_num = str(x['scan'])
        exp_charge = x['Title'].split('.')[-3]
        pparse_id = x['Title'].split('.')[-2]
        extended_name = x['Title'].split('.')[-1]

        return '.'.join([raw_name, scan_num, scan_num, exp_charge, pparse_id, extended_name])

    df['scan'] = df['Title'].apply(lambda x:x.split('.')[-4])
    ls_scan = list(df['scan'])
    if len(set(ls_scan)) == 1 and int(ls_scan[0]) == 0:
        df['scan'] = list(range(1, df.shape[0]+1))
        df['Title'] = df.apply(lambda x:set_title(x), axis=1)
        df.drop(columns=['scan'])
    
    return df


def result_xi2plink_linear_psm(plink_template_path, inpath, outpath):
    """将xi的单肽PSM结果转换为pLink格式"""

    df_plink_template = pd.read_csv(plink_template_path, sep=',')
    df_xi = pd.read_csv(inpath, sep=',')

    print('===@()', sys._getframe().f_code.co_name)
    # print(df_plink_template.columns)
    print('======PSM Number:', df_xi.shape[0])
    # print(df_xi.head()['isTT'])
    # print(df_xi.head().apply(xi2plink_title, axis=1))
    # print(df_xi.head().apply(lambda x:xi2plink_linear_seq_mod(x)[1], axis=1))
    # print(df_xi.head().apply(xi2plink_linear_seq_mod, axis=1)[1][1])

    # 去掉反库
    df_xi = df_xi[df_xi['isTT']]
    print('======PSM Number(only Target):', df_xi.shape[0])

    psm_num = df_xi.shape[0]
    df_plink = pd.DataFrame(columns=df_plink_template.columns)

    df_plink['Order'] = list(range(1, psm_num+1))
    df_plink['Title'] = list(df_xi.apply(xi2plink_title, axis=1))
    df_plink['Charge'] = list(df_xi['match charge'])
    df_plink['Peptide'] = list(df_xi.apply(lambda x:xi2plink_linear_seq_mod(x)[0], axis=1))
    df_plink['Peptide'] = df_plink['Peptide'].map(lambda x: x.replace("I", "L"))
    df_plink['Peptide_Type'] = ['Common'] * psm_num
    df_plink['Linker'] = ['null'] * psm_num
    df_plink['Modifications'] = list(df_xi.apply(lambda x:xi2plink_linear_seq_mod(x)[1], axis=1))

    df_plink = df_plink.fillna(0)

    df_plink.to_csv(outpath, index=False)


def xi_to_pLink(plink_template_path, inpath, fasta_path):
    """将xi的交联PSM结果转换为pLink格式,去掉非TT"""
    if pathlib.Path(inpath).with_suffix('.pLink.csv').exists():
        return

    df_plink_template = pd.read_csv(plink_template_path, sep=',')
    df_xi = pd.read_csv(inpath)

    print('===@()', sys._getframe().f_code.co_name)
    # print(df_plink_template.columns)
    print('======PSM Number:', df_xi.shape[0])
    # print(df_xi.head()['isTT'])
    # print(df_xi.head().apply(xi2plink_title, axis=1))
    # print(df_xi.head().apply(lambda x:xi2plink_linear_seq_mod(x)[1], axis=1))
    # print(df_xi.head().apply(xi2plink_linear_seq_mod, axis=1)[1][1])

    # 去掉反库
    df_xi = df_xi[df_xi['isTT']]
    df_xi = df_xi.reset_index(drop=True)
    print('======PSM Number(only Target):', df_xi.shape[0])

    psm_num = df_xi.shape[0]
    df_plink = pd.DataFrame(columns=df_plink_template.columns)

    df_plink['Order'] = list(range(1, psm_num+1))
    df_plink['Title'] = list(df_xi.apply(xi2plink_title, axis=1))
    df_plink = set_scan_order(df_plink)
    df_plink['Charge'] = df_xi['match charge']
    df_plink['Peptide'] = list(df_xi.apply(lambda x:xi2plink_crosslink_seq_mod(x)[0], axis=1))
    df_plink['Peptide'] = df_plink['Peptide'].map(lambda x: x.replace("I", "L"))
    df_plink['Peptide_Type'] = 'Cross-Linked'
    df_plink['Linker'] = df_xi['Crosslinker']
    df_plink['Modifications'] = list(df_xi.apply(lambda x:xi2plink_crosslink_seq_mod(x)[1], axis=1))

    df_plink['Protein_Type'] = list(df_xi.apply(xi2plink_proteintype, axis=1))
    df_plink['Precursor_Mass'] = df_xi['exp mass']+MASS_PROTON # pLink中为不带电质量+(H+)
    df_plink['Peptide_Mass'] = df_xi['match mass']+MASS_PROTON # pLink中为不带电质量+(H+)

    decoder = get_decoder_I2L_rev(pathlib.Path(fasta_path))
    df_plink['Proteins'] = df_plink['Peptide'].apply(
        lambda x: find_xl_seq_in_fasta_db(x, decoder))
    
    df_plink = df_plink.fillna(0)
    df_plink.to_csv(pathlib.Path(inpath).with_suffix('.pLink.csv'), index=False)
