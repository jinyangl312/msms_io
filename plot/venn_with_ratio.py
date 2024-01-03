from msms_io.plot.plot_venn import draw_venn, draw_weighted_venn2, draw_weighted_venn3
import pathlib
import pandas as pd


def draw_venn_between_dataset(level, data, dataset_list,
                              dataset_name_list,
                              folder_name, file_name, add_pQuant_ratio=False, weighted=True, file_suffix='png'):
    pathlib.Path.mkdir(pathlib.Path(
        f"./png/{folder_name}"), parents=True, exist_ok=True)

    if level == "spectrum":
        def func(x): return tuple(x[["_spectrum_id"]])
    elif level == "scan":
        def func(x): return tuple(x[["_scan_id"]])
    elif level == "PSM":
        def func(x): return tuple(
            x[["Peptide", "Modifications", "Charge", "_spectrum_id"]])
    elif level == "prec":
        def func(x): return tuple(x[["Peptide", "Modifications", "Charge"]])
    elif level == "pep":
        def func(x): return tuple(x[["Peptide", "Modifications"]])
    elif level == "seq":
        def func(x): return tuple(x[["Peptide"]])
    elif level == "rp":
        def func(x): return tuple(x[["_rp"]])
    elif level == "ppi":
        def func(x): return tuple(x[["_ppi"]])
    elif level == "protein":
        def func(x): return tuple(x[["_protein"]])
    else:
        assert False

    # Keep valid items
    dataset_1 = data[pd.notna(data[dataset_list[0]])]
    dataset_2 = data[pd.notna(data[dataset_list[1]])]
    if len(dataset_list) == 3:
        dataset_3 = data[pd.notna(data[dataset_list[2]])]
    else:
        assert len(dataset_list) == 2

    if len(dataset_list) == 3:
        if weighted:
            draw_func = draw_weighted_venn3
        else:
            draw_func = draw_venn
        if add_pQuant_ratio:
            # precursor level
            draw_func([
                set(dataset_1.apply(
                    func, axis=1).values.tolist()),
                set(dataset_2.apply(
                    func, axis=1).values.tolist()),
                set(dataset_3.apply(
                    func, axis=1).values.tolist()),
            ], [
                f"""{dataset_name_list[0]}
{'{:.2%}'.format(len(dataset_1[dataset_1["has_pQuant_evidence"] == False]) / (len(dataset_1) + 1e-6))}""",
                f"""{dataset_name_list[1]}
{'{:.2%}'.format(len(dataset_2[dataset_2["has_pQuant_evidence"] == False]) / (len(dataset_2) + 1e-6))}""",
                f"""{dataset_name_list[2]}
{'{:.2%}'.format(len(dataset_3[dataset_3["has_pQuant_evidence"] == False]) / (len(dataset_3) + 1e-6))}""",
            ], f"png/{folder_name}/{file_name}.{file_suffix}")
        else:
            draw_func([
                set(dataset_1.apply(
                    func, axis=1).values.tolist()),
                set(dataset_2.apply(
                    func, axis=1).values.tolist()),
                set(dataset_3.apply(
                    func, axis=1).values.tolist()),
            ], [
                dataset_name_list[0],
                dataset_name_list[1],
                dataset_name_list[2],
            ], f"png/{folder_name}/{file_name}.{file_suffix}")
    else:
        if weighted:
            draw_func = draw_weighted_venn2
        else:
            draw_func = draw_venn
        if add_pQuant_ratio:
            # precursor level
            draw_func([
                set(dataset_1.apply(
                    func, axis=1).values.tolist()),
                set(dataset_2.apply(
                    func, axis=1).values.tolist()),
            ], [
                f"""{dataset_name_list[0]}
{'{:.2%}'.format(len(dataset_1[dataset_1["has_pQuant_evidence"] == False]) / (len(dataset_1) + 1e-6))}""",
                f"""{dataset_name_list[1]}
{'{:.2%}'.format(len(dataset_2[dataset_2["has_pQuant_evidence"] == False]) / (len(dataset_2) + 1e-6))}""",
            ], f"png/{folder_name}/{file_name}.{file_suffix}")
        else:
            draw_func([
                set(dataset_1.apply(
                    func, axis=1).values.tolist()),
                set(dataset_2.apply(
                    func, axis=1).values.tolist()),
            ], [
                dataset_name_list[0],
                dataset_name_list[1],
            ], f"png/{folder_name}/{file_name}.{file_suffix}")
