import pandas as pd
import re
import numpy as np


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
