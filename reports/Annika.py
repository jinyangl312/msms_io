# -*- coding: UTF-8 -*-
# @Date     : 28rd Feb, 2022
# @Author   : Northblue

"""将MSAnnika的交联CSM结果转换为pLink格式
    经过FDR过滤
    去掉非TT
"""

import sys
import os
import pandas as pd
from msms_io.fasta.find_protein import get_decoder_I2L_rev
from msms_io.fasta.find_protein import find_xl_seq_in_fasta_db
from msms_io.reports.sort_order import *
from msms_io.reports.functions_xl import *
import pathlib
import re
from msms_io.fasta.find_protein import find_xl_seq_in_fasta_db, is_xl_seq_in_same_syn_group
from msms_io.fasta.find_protein import get_decoder_I2L_rev, get_decoder_I2L


MAP_TRANS_MODIFICAIONT = {'Carbamidomethyl':'Carbamidomethyl[C]', 'Oxidation':'Oxidation[M]'}
# , 'leikermono': 'Leiker_Mono[K]', 'dssmono': 'DSS_Mono[K]', 'dssooh': 'DSSO_H2O[K]', 'dssoh2': 'DSSO_NH2[K]', 'dsbuoh': 'DSBU_H2O[K]', 'dsbuh2': 'DSBU_NH2[K]', 'dsbunh2': 'DSBU_NH2[K]', 'dsbuloop': 'DSBU_LOOP[K]'}

H_MASS = 1.00782503214

delimiter = '\t'

def msannika2plink_title(ser):
    """MSAnnika中提出plink的谱图title"""

    mgf_name = ser['Spectrum File'][:-4]
    scan_id = str(ser['First Scan'])
    exp_charge = str(ser['Charge'])
    pparse_id = '0'
    extended_name = 'dta'

    return '.'.join([mgf_name, scan_id, scan_id, exp_charge, pparse_id, extended_name])

def msannika2plink_csm_pep(ser):
    """MSAnnika中提出plink的peptide pair
    序列+交联位点
    两条序列之间没有排序
    N端修饰位点为1
    """

    seq1 = ser['Sequence A']
    seq2 = ser['Sequence B']

    link_pos1 = ser['Crosslinker Position A'] #从1开始
    link_pos2 = ser['Crosslinker Position B']

    pep_str = seq1 + '(' + str(link_pos1) + ')-' + seq2 + '(' + str(link_pos2) + ')'

    return pep_str

def msannika2plink_csm_mod(ser):
    """MSAnnika中提出plink的modifications
    两条序列之间没有排序
    C1(Carbamidomethyl);C8(Carbamidomethyl);K10(DSBU)
    Nterm(DSBU);M8(Oxidation)
    K3(DSSO);Cterm(Oxidation) GQKPLM # 会有标位N端和C端修饰
    """

    seq1 = ser['Sequence A']
    seq2 = ser['Sequence B']
    mods1 = ser['Modifications A']
    mods2 = ser['Modifications B']

    def get_site2mod(mods, seq):
        site2mod = {}
        for mod_str in mods.split(';'):
            mod_name = mod_str.split('(')[1][:-1]
            if mod_name in MAP_TRANS_MODIFICAIONT:
                site_str = mod_str.split('(')[0]
                if site_str == 'Cterm':
                    mod_site = len(seq)
                elif site_str == 'Nterm':
                    mod_site = 1
                else:
                    mod_site = int(mod_str.split('(')[0][1:])
                site2mod[mod_site] = MAP_TRANS_MODIFICAIONT[mod_name]
        return site2mod
        
    site2mod1 = get_site2mod(mods1,seq1)
    site2mod2 = get_site2mod(mods2,seq2)

    mod_ls = []
    for k, v in site2mod1.items():
        mod_ls.append(v + '('+ str(k) + ')')
    for k, v in site2mod2.items():
        mod_ls.append(v + '('+ str(k+len(seq1)+3) + ')')
    
    if mod_ls:
        return ';'.join(mod_ls)
    else:
        return 'null'

def msannika2plink_proteintype(ser):
    """获得intra/inter"""

    fdr_group = ser['Crosslink Type'].strip()
    fdr_group2protein_type = {'Intra':'Intra-Protein', 'Inter':'Inter-Protein'}

    return fdr_group2protein_type[fdr_group]

def get_fpath_csm(dpath):
    for fname in list(os.walk(dpath))[0][2]:
        if fname[-9:] == '_CSMs.txt':
            return os.path.join(dpath, fname)
def get_fpath_crosslinks(dpath):
    for fname in list(os.walk(dpath))[0][2]:
        if fname[-15:] == '_Crosslinks.txt':
            return os.path.join(dpath, fname)
def get_fpath_xlsummary(dpath):
    for fname in list(os.walk(dpath))[0][2]:
        if fname[-21:] == '_CrosslinkSummary.txt':
            return os.path.join(dpath, fname)

def get_fdr_num(df, level, fdr):
    """获得FDR过滤后的鉴定数目"""
    ser = df[df['File name']=='Total']
    
    if fdr in ['0.05', '0.01']:
        if fdr == '0.05':
            col = '# Medium Confidence '
        if fdr == '0.01':
            col = '# High Confidence '
        if level == 'rp':
            col += 'Cross-links'
        if level == 'csm':
            col += 'CSMs'
        return list(ser[col])[0]
    else:
        if level == 'rp':
            return list(ser['# Cross-links'])[0]
        elif level == 'csm':
            return list(ser['# CSMs'])[0]

def filter_csm(df, num):
    """经过FDR过滤，去掉非TT"""
    df = df.sort_values(['Combined Score'], ascending=False)
    df = df.iloc[:num]
    df = df[df['Alpha T/D'] == 'T']
    df = df[df['Beta T/D'] == 'T']
    df = df.reset_index(drop=True)
    return df

def Annika_to_pLink(plink_template_path, inpath, fasta_path, fdr='0.01'):
    """将MSAnnika的交联CSM结果转换为pLink格式
    经过FDR过滤
    去掉非TT
    """
    if (pathlib.Path(inpath)/'csm.pLink.csv').exists():
        return

    df_plink_template = pd.read_csv(plink_template_path, sep=',')

    fpath_csm = get_fpath_csm(inpath)
    fpath_xlsummary = get_fpath_xlsummary(inpath)
    df_csm = pd.read_csv(fpath_csm, sep=delimiter)
    df_xlsummary = pd.read_csv(fpath_xlsummary, sep=delimiter)

    # 经过FDR过滤，去掉非TT
    df_csm = filter_csm(df_csm, get_fdr_num(df_xlsummary, 'csm', fdr))

    print('======CSM Number:', df_csm.shape[0])

    csm_num = df_csm.shape[0]
    df_plink = pd.DataFrame(columns=df_plink_template.columns)

    df_plink['Order'] = list(range(1, csm_num+1))
    df_plink['Title'] = list(df_csm.apply(msannika2plink_title, axis=1))
    # Should be uCSM
    df_plink["_scan_id"] = df_plink["Title"].swifter.apply(
        lambda x: re.search("(?<=^).*?\.\d+\.\d+(?=\.\d+\.\d+\.dta\;?)", x).group())
    assert df_plink['_scan_id'].duplicated().sum() == 0, 'rCSM rather than uCSM!'
    df_plink.drop(columns=['_scan_id'], inplace=True)
    df_plink['Charge'] = df_csm['Charge']

    df_plink['Peptide'] = list(df_csm.apply(msannika2plink_csm_pep, axis=1))
    df_plink['Peptide'] = df_plink['Peptide'].map(lambda x: x.replace("I", "L"))
    df_plink['Peptide_Type'] = 'Cross-Linked'
    df_plink['Linker'] = df_csm['Crosslinker']
    df_plink['Modifications'] = list(df_csm.apply(msannika2plink_csm_mod, axis=1))
    df_plink['Modifications'] = df_plink['Modifications'].map(sort_modification)

    df_plink['Protein_Type'] = list(df_csm.apply(msannika2plink_proteintype, axis=1))
    df_plink['Precursor_Mass'] = (df_csm['m/z [Da]']-H_MASS)*df_csm['Charge']+H_MASS # pLink中为不带电质量+H
    df_plink['Peptide_Mass'] = df_csm['MH+ [Da]'] # pLink中为不带电质量+H
    df_plink['Precursor_Mass_Error(ppm)'] = df_csm['DeltaM [ppm]']# pLink与MS Annika误差是带正负的
    df_plink['Score'] = df_csm['Combined Score']

    decoder = get_decoder_I2L_rev(pathlib.Path(fasta_path))
    df_plink['Proteins'] = df_plink['Peptide'].apply(
        lambda x: find_xl_seq_in_fasta_db(x, decoder))

    df_plink = df_plink.fillna(0)

    print('======Protein_Type:\n', df_plink['Protein_Type'].value_counts())

    df_plink.to_csv(pathlib.Path(inpath)/'csm.pLink.csv', index=False)
    return df_plink


def load_pp_Annika(project_dir, fasta_path):    
    # Fields  need for evaluation:
    # Peptide, 
    # Proteins, Protein_Type

    project_dir = pathlib.Path(project_dir)
    fasta_path = pathlib.Path(fasta_path)

    if (project_dir/'pp_preprocessed.csv').exists():
        return pd.read_csv(project_dir/'pp_preprocessed.csv')

    data = pd.read_csv(get_fpath_crosslinks(project_dir), delimiter='\t').fillna('')
    data = data[data['Confidence'] == 'High']
    data = data[data['Decoy'] == False]
    data.reset_index(drop=True, inplace=True)

    # FAAYA [K] AYPQEAAEFTR -> FAAYAKAYPQEAAEFTR
    data['Peptide_1'] = data['Sequence A'].apply(
        lambda x: re.sub(' |\[|\]', '', x)
    )
    data['Peptide_2'] = data['Sequence B'].apply(
        lambda x: re.sub(' |\[|\]', '', x)
    )

    data['Peptide'] = data[['Peptide_1', 'Position A', 'Peptide_2', 'Position B']].apply(
        lambda x: f'{x[0]}({x[1]})-{x[2]}({x[3]})', axis=1
    )
    data['Peptide'] = data['Peptide'].map(lambda x: x.replace("I", "L"))

    data['Modifications'] = list(data.apply(msannika2plink_csm_mod, axis=1))
    data['Modifications'] = data['Modifications'].map(sort_modification)

    data[['Peptide', 'Modifications']] = list(data[['Peptide', 'Modifications']].apply(
        lambda x: sort_alpha_beta_order(*x), axis=1))

    data['Peptide_Type'] = 'Cross-Linked'

    sequence_decoder = get_decoder_I2L_rev(fasta_path)
    data["Proteins"] = data["Peptide"].apply(
        lambda x: find_xl_seq_in_fasta_db(x, sequence_decoder))
    data["Proteins"] = data["Proteins"].map(
        sort_site_order)
        
    data['Protein_Type'] = data['Crosslink Type'].map({
        'Inter': 'Inter-Protein',
        'Intra': 'Intra-Protein',
    })

    data.to_csv(project_dir/'pp_preprocessed.csv', index=False)
    return data


def load_rp_Annika(project_dir, fasta_path):    
    # Fields  need for evaluation:
    # Peptide, 
    # Proteins, Protein_Type

    project_dir = pathlib.Path(project_dir)
    fasta_path = pathlib.Path(fasta_path)

    if (project_dir/'rp_preprocessed.csv').exists():
        return pd.read_csv(project_dir/'rp_preprocessed.csv')

    data = pd.read_csv(get_fpath_crosslinks(project_dir), delimiter='\t').fillna('')
    data = data[data['Confidence'] == 'High']
    data = data[data['Decoy'] == False]
    data.reset_index(drop=True, inplace=True)

    # FAAYA [K] AYPQEAAEFTR -> FAAYAKAYPQEAAEFTR
    data['Peptide_1'] = data['Sequence A'].apply(
        lambda x: re.sub(' |\[|\]', '', x)
    )
    data['Peptide_2'] = data['Sequence B'].apply(
        lambda x: re.sub(' |\[|\]', '', x)
    )

    data['Peptide'] = data[['Peptide_1', 'Position A', 'Peptide_2', 'Position B']].apply(
        lambda x: f'{x[0]}({x[1]})-{x[2]}({x[3]})', axis=1
    )
    data['Peptide'] = data['Peptide'].map(lambda x: x.replace("I", "L"))

    data['Modifications'] = list(data.apply(msannika2plink_csm_mod, axis=1))
    data['Modifications'] = data['Modifications'].map(sort_modification)

    data[['Peptide', 'Modifications']] = list(data[['Peptide', 'Modifications']].apply(
        lambda x: sort_alpha_beta_order(*x), axis=1))

    data['Peptide_Type'] = 'Cross-Linked'

    sequence_decoder = get_decoder_I2L_rev(fasta_path)
    data["Proteins"] = data["Peptide"].apply(
        lambda x: find_xl_seq_in_fasta_db(x, sequence_decoder))
    data["Proteins"] = data["Proteins"].map(
        sort_site_order)
        
    data['Protein_Type'] = data['Crosslink Type'].map({
        'Inter': 'Inter-Protein',
        'Intra': 'Intra-Protein',
    })

    data.to_csv(project_dir/'rp_preprocessed.csv', index=False)
    return data
