from logging import handlers
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import pathlib
import venn
plt.rcParams['font.family'] = 'Times New Roman'
# plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
# config = {
#     "font.family": 'serif',
#     # "font.size": 15,
#     # "mathtext.fontset": 'stix',
#     "font.serif": ['SimHei'],
#     'axes.unicode_minus': False,  # 处理负号，即-号
# }
# plt.rcParams.update(config)

# https://zhuanlan.zhihu.com/p/501395717


def draw_weighted_venn2(set_list, label_list, save_name=None):
    # https://pypi.org/project/matplotlib-venn/
    out = venn2(set_list, label_list)
    for text in out.set_labels:
        text.set_fontsize(15)
    for x in range(len(out.subset_labels)):
        if out.subset_labels[x] is not None:
            out.subset_labels[x].set_fontsize(16)

    if not save_name == None:
        plt.savefig(save_name, dpi=512, bbox_inches='tight')
    else:
        plt.show()
    plt.close()


def draw_weighted_venn3(set_list, label_list, save_name=None):
    # https://pypi.org/project/matplotlib-venn/
    out = venn3(set_list, label_list)
    for text in out.set_labels:
        text.set_fontsize(15)
    for x in range(len(out.subset_labels)):
        if out.subset_labels[x] is not None:
            out.subset_labels[x].set_fontsize(14)
    plt.gca.legend_ =None

    if not save_name == None:
        plt.savefig(save_name, dpi=512, bbox_inches='tight')
    else:
        plt.show()
    plt.close()


def draw_venn(set_list, label_list, save_name=None, fmt="size"):
    # https://github.com/LankyCyril/pyvenn/blob/master/pyvenn-demo.ipynb
    venn_dict = dict(zip(label_list, set_list))
    if fmt == "percentage":
        venn.venn(venn_dict, fmt="{percentage:.1f}%",
                  fontsize=15, legend_loc="best")
    else:
        venn.venn(venn_dict, fmt="{size}", fontsize=15, legend_loc="best")

    if not save_name == None:
        plt.savefig(save_name, dpi=512, bbox_inches='tight')
    else:
        plt.show()
    plt.close()
