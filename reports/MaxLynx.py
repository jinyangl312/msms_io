# -*- coding: UTF-8 -*-
# @Date     : 28rd Feb, 2022
# @Author   : Northblue

"""将maxlynx的结果转换为pLink2结果格式
    去掉非TT
    只保留交联
    表中无交联剂信息
    两条序列之间按照原MaxLynx输出顺序排序
"""

import sys
import pandas as pd
from msms_io.fasta.find_protein import get_decoder_I2L_rev
from msms_io.fasta.find_protein import find_xl_seq_in_fasta_db
from msms_io.reports.sort_order import sort_modification
import pathlib


MAP_TRANS_MODIFICAIONT = {'Oxidation (M)': 'Oxidation[M]'}
MAP_FIX_MOD = {'C': 'Carbamidomethyl[C]'} 
# , 'leikermono': 'Leiker_Mono[K]', 'dssmono': 'DSS_Mono[K]', 'dssooh': 'DSSO_H2O[K]', 'dssoh2': 'DSSO_NH2[K]', 'dsbuoh': 'DSBU_H2O[K]', 'dsbuh2': 'DSBU_NH2[K]', 'dsbunh2': 'DSBU_NH2[K]', 'dsbuloop': 'DSBU_LOOP[K]'}

H_MASS = 1.00782503214

# 带mono修饰的直接当做单肽处理

delimiter = '\t'

def maxlynx2plink_title(ser):
    """MaxLynx中提出plink的谱图title"""

    mgf_name = ser['Raw file']
    scan_id = str(ser['Scan number'])
    exp_charge = str(ser['Charge'])
    pparse_id = '0'
    extended_name = 'dta'

    return '.'.join([mgf_name, scan_id, scan_id, exp_charge, pparse_id, extended_name])

def maxlynx2plink_csm_pep(ser):
    """MaxLynx中提出plink的peptide pair
    序列+交联位点
    两条序列之间没有排序
    """

    seq1 = ser['Sequence1']
    seq2 = ser['Sequence2']

    link_pos1 = ser['Peptide index of Crosslink 1'] #从1开始
    link_pos2 = ser['Peptide index of Crosslink 2']

    pep_str = seq1 + '(' + str(link_pos1) + ')-' + seq2 + '(' + str(link_pos2) + ')'

    return pep_str

def maxlynx2plink_csm_mod(ser):
    """MaxLynx中提出plink的modifications
    两条序列之间没有排序
    """

    seq1 = ser['Sequence1']
    seqmod1 = ser['Modified sequence1']
    seqmod2 = ser['Modified sequence2']
    seqmod1 = seqmod1[1:-1]
    seqmod2 = seqmod2[1:-1]

    mp_site_modstr1 = {}
    mp_site_modstr2 = {}

    def get_site2mod(seqmod):
        site2mod = {}
        site = 0
        i = 0
        while i < len(seqmod):
            if seqmod[i] == '(':
                mod_str = []
                i += 1
                while seqmod[i] != ')':
                    mod_str.append(seqmod[i])
                    i += 1
                mod_str.append(seqmod[i]) # 指向(M)的后括号
                i += 1 # 指向(Oxidation (M))的后括号
                site2mod[site] = MAP_TRANS_MODIFICAIONT[''.join(mod_str)]
            else:
                site += 1
                # 固定修饰
                if seqmod[i] in MAP_FIX_MOD:
                    site2mod[site] = MAP_FIX_MOD[seqmod[i]]
            i += 1
        return site2mod
    
    mp_site_modstr1 = get_site2mod(seqmod1)
    mp_site_modstr2 = get_site2mod(seqmod2)

    mod_ls = []
    for k, v in mp_site_modstr1.items():
        mod_ls.append(v + '('+ str(k) + ')')
    for k, v in mp_site_modstr2.items():
        mod_ls.append(v + '('+ str(k+len(seq1)+3) + ')')
    
    if mod_ls:
        return ';'.join(mod_ls)
    else:
        return 'null'

def maxlynx2plink_proteintype(ser):
    """获得intra/inter"""

    fdr_group = ser['Crosslink product type'].strip()
    fdr_group2protein_type = {'Intra-protein link':'Intra-Protein', 'Inter-protein link':'Inter-Protein'}

    return fdr_group2protein_type[fdr_group]

def MaxLynx_to_pLink(plink_template_path, inpath, fasta_path):
    """将MaxLynx的交联CSM结果转换为pLink格式
    去掉非TT
    只保留交联
    表中无交联剂信息
    """
    if pathlib.Path(inpath).with_suffix('.pLink.csv').exists():
        return

    df_plink_tmp = pd.read_csv(plink_template_path, sep=',')
    df_maxlynx = pd.read_csv(inpath, sep=delimiter)

    print('======PSM Number:', df_maxlynx.shape[0])
    assert sum(df_maxlynx[['Raw file', 'Scan number']].duplicated()) == 0, 'Find redundant PSM!'

    # 去掉反库，只保留交联
    df_maxlynx = df_maxlynx[df_maxlynx['Decoy'] == 'forward']
    df_maxlynx = df_maxlynx[df_maxlynx['Crosslink product type'].isin(['Intra-protein link', 'Inter-protein link'])]
    df_maxlynx = df_maxlynx.reset_index(drop=True)    
    print('======CSM Number(only XL & TT):', df_maxlynx.shape[0])

    csm_num = df_maxlynx.shape[0]
    df_plink = pd.DataFrame(columns=df_plink_tmp.columns)

    df_plink['Order'] = list(range(1, csm_num+1))
    df_plink['Title'] = list(df_maxlynx.apply(maxlynx2plink_title, axis=1))
    df_plink['Charge'] = df_maxlynx['Charge']
    df_plink['Peptide'] = list(df_maxlynx.apply(maxlynx2plink_csm_pep, axis=1))
    df_plink['Peptide'] = df_plink['Peptide'].map(lambda x: x.replace("I", "L"))
    df_plink['Peptide_Type'] = 'Cross-Linked'
    df_plink['Modifications'] = list(df_maxlynx.apply(maxlynx2plink_csm_mod, axis=1))

    # 必须要有list，因为df会按照index对齐
    df_plink['Protein_Type'] = list(df_maxlynx.apply(maxlynx2plink_proteintype, axis=1))
    df_plink['Precursor_Mass'] = list((df_maxlynx['m/z']-H_MASS)*df_maxlynx['Charge']+H_MASS) # pLink中为不带电质量+H
    df_plink['Peptide_Mass'] = df_maxlynx['Mass']+H_MASS
    df_plink['Precursor_Mass_Error(ppm)'] = df_maxlynx['Mass error [ppm]']
    df_plink['Score'] = df_maxlynx['Score']

    decoder = get_decoder_I2L_rev(pathlib.Path(fasta_path))
    df_plink['Proteins'] = df_plink['Peptide'].apply(
        lambda x: find_xl_seq_in_fasta_db(x, decoder))
    
    df_plink = df_plink.fillna(0)
    df_plink.to_csv(pathlib.Path(inpath).with_suffix('.pLink.csv'), index=False)

    return df_plink
