import pandas as pd
import re
import numpy as np
import functools


def not_empty(s):
    return s and s.strip()


def gen_df_from_pN_res(PSM_res, keep_top1=True):
    lines = PSM_res.strip('\n').split('\n')
    title = lines[0].split('\t')[1]

    if len(lines) == 1:
        candidates = [[np.nan]*2]
    else:
        candidates = [x.split('\t')[1:3] for x in lines[1:]]
        if keep_top1:
            candidates = candidates[0:1]

    df = pd.DataFrame(candidates, columns=['Sequence', 'Score'])
    df['Title'] = title
    return df


def process_pFind_res_line(res_line):
    def compare_mods(x, y):
        return int(re.search("\d+(?=\,)", x).group()) - int(re.search("\d+(?=\,)", y).group())

    def sort_modifications_by_site(modification):
        if modification == '':
            return modification
        mods = re.split("\;", modification.strip(';'))
        return ";".join(sorted(mods, key=functools.cmp_to_key(compare_mods)))+";"

    res_list = res_line.strip('\n').split('\n')
    mz = float(res_list[0].split('\t')[0])
    charge = int(res_list[0].split('\t')[1])
    msms_title = res_list[1]
    PSM_list = res_list[2:]
    if len(PSM_list) == 0:
        return (msms_title, mz, charge, [])

    PSM_list = [(x.split('\t')) for x in PSM_list]
    PSM_res = []
    for PSM_line in PSM_list:
        sequence = PSM_line[1].replace('I', 'L')
        final_score = float(PSM_line[3])

        modifications_num = int(PSM_line[6])
        modifications = []
        for i in range(modifications_num):
            mod_site = int(PSM_line[7+2*i])
            mod_type = PSM_line[8+2*i].strip('#0')
            modifications.append(f'{mod_site},{mod_type};')

        cleavage_type = int(PSM_line[7+2*modifications_num])
        td_type = int(PSM_line[8+2*modifications_num])

        PSM_res.append((sequence, sort_modifications_by_site(
            ''.join(modifications)), cleavage_type, td_type))

    return (msms_title, mz, charge, PSM_res)
