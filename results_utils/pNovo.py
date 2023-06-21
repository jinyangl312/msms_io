import pandas as pd
from .functions import *
import functools
import multiprocessing
import os


def load_pNovo_res(res_path, keep_top1=True):
    '''
    Load columns from pNovo.res as pd.DataFrame
    '''

    '''
    S1	谱图名字
    P1	序列	打分	修饰丰度	母离子误差	路径排名	主干离子打分	内部离子打分	b离子序列连续性	y离子序列连续性	酶切位点打分	碎裂离子误差的标准差	碎裂离子误差的最大偏差	电荷特征	bm25打分	Spearman相似度	Gap打分	N端前两个氨基酸的Gap打分
    '''

    with open(res_path, 'r') as f:
        all_PSMs = f.read().strip('\n').split('\n\n')
    with multiprocessing.Pool(os.cpu_count()) as pool:
        all_df = pool.map(functools.partial(
            gen_df_from_pN_res, keep_top1=keep_top1), all_PSMs)

    # all_df = list(map(functools.partial(
    #     gen_df_from_pN_res, keep_top1=keep_top1), all_PSMs))

    res = pd.concat(all_df, ignore_index=True)

    return res
