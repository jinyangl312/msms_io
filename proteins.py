import pandas as pd
from .utils.xl_utils import rmv_xl_site, split_xl_protein
from .load_results import *
from evaluation_scripts.deduction import group_by_col


def get_protein_trace_table_from_pF(spectra, quant=None, protein_path=None, suffix=None, filtered=True):
    '''
    Load PSM from pFind and return protein trace table.
    '''
    PSM_res = get_PSM_from_pF(spectra)
    if not filtered:
        PSM_res = PSM_res[PSM_res["Q-value"].apply(lambda x: x < 1)]
        PSM_res = PSM_res[PSM_res["Target/Decoy"].apply(lambda x: x == "target")]
    PSM_res = PSM_res[PSM_res["Proteins"].apply(lambda x: len(re.split("/", x))) == 2]
    PSM_res["Proteins"] = PSM_res["Proteins"].apply(lambda x: x.strip("/| "))

    if quant != None:
        quant_res = get_PSM_from_pQ(quant)[["Name_MS2", "Ratio_Sample2/Sample1"]]
        quant_res = quant_res.rename(columns={"Name_MS2": "File_Name"})
        PSM_res = PSM_res.merge(quant_res,
            left_on=["File_Name"],
            right_on=["File_Name"],
            suffixes=('', ''),
            how="left")

        PSM_res["pQuant_evidence"] = PSM_res["Ratio_Sample2/Sample1"].map(
            lambda x: True if x > 1/1024 and x < 1024 else False)
        PSM_res["Ratio_Sample2/Sample1"] = PSM_res["Ratio_Sample2/Sample1"].apply(str)

        protein_data = group_by_col(
                PSM_res,
                by=["Proteins"],
                keep_min_cols=[],
                keep_max_cols=[
                    'pQuant_evidence',
                ],
                concat_str_cols=[
                    "File_Name",
                    "Ratio_Sample2/Sample1",
                    ]
            )[["File_Name", "pQuant_evidence", "Ratio_Sample2/Sample1"]]
        protein_data = protein_data.rename(columns={
            "File_Name": f"PSM{suffix}",
            "pQuant_evidence": f"pQuant{suffix}",
            "Ratio_Sample2/Sample1": f"pQuant_ratio{suffix}",
            })
    else:
        protein_data = group_by_col(
                PSM_res,
                by=["Proteins"],
                keep_min_cols=[],
                keep_max_cols=[],
                concat_str_cols=[
                    "File_Name",
                    ]
            )[["File_Name"]]
        protein_data = protein_data.rename(columns={
            "File_Name": f"PSM{suffix}",
            })

    if filtered and protein_path != None:
        proteins_df = get_protein_from_pF(
            res_path=protein_path,
            filter=True,
            keep_subset=True)
        proteins_set = set(proteins_df[pd.notna(proteins_df[f"ID"])]["AC"].tolist())
        proteins_set.update(
            set(proteins_df[~pd.notna(proteins_df[f"ID"])]["Score"].tolist()))
        protein_data = protein_data.reset_index()
        protein_data = protein_data[protein_data["Proteins"].apply(lambda x: x in proteins_set)].set_index("Proteins")

    return protein_data


def get_protein_trace_table_from_pL(spectra, quant=None, suffix=None, xl=False, filtered=True):
    '''
    Load PSM from pLink and return protein trace table.
    '''
    PSM_res = pd.read_csv(spectra)
    PSM_res = PSM_res[PSM_res["Proteins"].apply(lambda x: len(re.split("/", x))) == 2]

    if filtered:
        if xl:
            intra = PSM_res[PSM_res["Protein_Type"] == "Intra-Protein"]
            intra["Proteins"] = intra["Proteins"].apply(rmv_xl_site).apply(lambda x: re.split("-", x)[0])
            inter = PSM_res[PSM_res["Protein_Type"] == "Inter-Protein"]
            inter["Proteins1"] = inter["Proteins"].apply(rmv_xl_site).apply(lambda x: re.split("-", x)[0])
            inter["Proteins2"] = inter["Proteins"].apply(rmv_xl_site).apply(lambda x: re.split("-", x)[1])
            inter = inter.drop(columns=["Proteins", "Proteins2"]).rename(columns={"Proteins1": "Proteins"}).append(
                inter.drop(columns=["Proteins", "Proteins1"]).rename(columns={"Proteins2": "Proteins"}),
                ignore_index=True,
            )
            PSM_res = inter.append(intra, ignore_index=True)
        else:
            PSM_res["Proteins"] = PSM_res["Proteins"].apply(rmv_xl_site)
    else:
        PSM_res = PSM_res[PSM_res["Target_Decoy"] == 2]
        xl = PSM_res[PSM_res["Peptide_Type"] == 3]
        inter = xl[xl["Protein_Type"] == 2]
        intra = xl[xl["Protein_Type"] == 1]
        rlm = PSM_res[PSM_res["Peptide_Type"] != 3]

        intra["Proteins"] = intra["Proteins"].apply(rmv_xl_site).apply(lambda x: re.split("-", x)[0])
        inter["Proteins1"] = inter["Proteins"].apply(rmv_xl_site).apply(lambda x: re.split("-", x)[0])
        inter["Proteins2"] = inter["Proteins"].apply(rmv_xl_site).apply(lambda x: re.split("-", x)[1])
        inter = inter.drop(columns=["Proteins", "Proteins2"]).rename(columns={"Proteins1": "Proteins"}).append(
            inter.drop(columns=["Proteins", "Proteins1"]).rename(columns={"Proteins2": "Proteins"}),
            ignore_index=True,
        )
        rlm["Proteins"] = rlm["Proteins"].apply(rmv_xl_site)
        PSM_res = inter.append(
            intra, ignore_index=True).append(
                rlm, ignore_index=True)
    PSM_res["Proteins"] = PSM_res["Proteins"].apply(lambda x: x.strip("/| ").strip())

    if quant != None:
        quant_res = get_PSM_from_pQ(quant)[["Name_MS2", "Ratio_Sample2/Sample1"]]
        quant_res = quant_res.rename(columns={"Name_MS2": "Title"})
        PSM_res = PSM_res.merge(quant_res,
            left_on=["Title"],
            right_on=["Title"],
            suffixes=('', ''),
            how="left")

        PSM_res["pQuant_evidence"] = PSM_res["Ratio_Sample2/Sample1"].map(
            lambda x: True if x > 1/1024 and x < 1024 else False)
        PSM_res["Ratio_Sample2/Sample1"] = PSM_res["Ratio_Sample2/Sample1"].apply(str)

        protein_data = group_by_col(
                PSM_res,
                by=["Proteins"],
                keep_min_cols=[],
                keep_max_cols=[
                    'pQuant_evidence',
                ],
                concat_str_cols=[
                    "Title",
                    "Ratio_Sample2/Sample1",
                    ]
            )[["Title", "pQuant_evidence", "Ratio_Sample2/Sample1"]]
        protein_data = protein_data.rename(columns={
            "Title": f"PSM{suffix}",
            "pQuant_evidence": f"pQuant{suffix}",
            "Ratio_Sample2/Sample1": f"pQuant_ratio{suffix}",
            })
    else:
        protein_data = group_by_col(
                PSM_res,
                by=["Proteins"],
                keep_min_cols=[],
                keep_max_cols=[],
                concat_str_cols=[
                    "Title",
                    ]
            )[["Title"]]
        protein_data = protein_data.rename(columns={
            "Title": f"PSM{suffix}",
            })
    return protein_data


def get_protein_set_from_list(proteins_list):
    '''Turn pandas serier to python set'''
    proteins_list = list(map(lambda x: x.strip("/").split("/"), proteins_list))
    protein_set = set()
    for proteins_line in proteins_list:
        for protein in proteins_line:
            protein = protein.strip()
            if protein not in protein_set:
                protein_set.add(protein)
    return protein_set


def get_protein_set_from_pFind(pFind_res_path):
    '''Load and extract protein set from pFind'''
    PSM_res = get_PSM_from_pF(pFind_res_path)
    protein_set = get_protein_set_from_list(list(set(
        PSM_res["Proteins"].tolist())))
    return protein_set


def get_xl_protein_set_from_pLink(pLink_res_prefix):
    '''Load and extract xl protein set from pLink'''
    xl_results = pd.read_csv(f"{pLink_res_prefix}.filtered_cross-linked_spectra.csv")

    xl_protein_set = get_protein_set_from_list(list(set(
            xl_results["Proteins"].apply(rmv_xl_site).apply(split_xl_protein).tolist())))
    return xl_protein_set


def get_non_xl_protein_set_from_pLink(pLink_res_prefix):
    '''Load and extract non-xl protein set from pLink'''
    mono_results = pd.read_csv(f"{pLink_res_prefix}.filtered_mono-linked_spectra.csv")
    loop_results = pd.read_csv(f"{pLink_res_prefix}.filtered_loop-linked_spectra.csv")
    regular_results = pd.read_csv(f"{pLink_res_prefix}.filtered_regular_spectra.csv")

    mono_protein_set = get_protein_set_from_list(list(set(
            mono_results["Proteins"].apply(rmv_xl_site).tolist())))
    loop_protein_set = get_protein_set_from_list(list(set(
            loop_results["Proteins"].apply(rmv_xl_site).tolist())))
    regular_protein_set = get_protein_set_from_list(list(set(
            regular_results["Proteins"].tolist())))

    protein_set = set()
    protein_set.update(mono_protein_set)
    protein_set.update(loop_protein_set)
    protein_set.update(regular_protein_set)
    return protein_set


def get_regular_protein_set_from_pLink(pLink_res_prefix):
    '''Load and extract non-xl protein set from pLink'''
    regular_results = pd.read_csv(f"{pLink_res_prefix}.filtered_regular_spectra.csv")

    regular_protein_set = get_protein_set_from_list(list(set(
            regular_results["Proteins"].tolist())))
    return regular_protein_set


def get_mono_protein_set_from_pLink(pLink_res_prefix):
    '''Load and extract non-xl protein set from pLink'''
    mono_results = pd.read_csv(f"{pLink_res_prefix}.filtered_mono-linked_spectra.csv")

    mono_protein_set = get_protein_set_from_list(list(set(
            mono_results["Proteins"].apply(rmv_xl_site).tolist())))
    return mono_protein_set


def get_loop_protein_set_from_pLink(pLink_res_prefix):
    '''Load and extract non-xl protein set from pLink'''
    loop_results = pd.read_csv(f"{pLink_res_prefix}.filtered_loop-linked_spectra.csv")

    loop_protein_set = get_protein_set_from_list(list(set(
            loop_results["Proteins"].apply(rmv_xl_site).tolist())))
    return loop_protein_set


def get_protein_set_from_pLink(pLink_res_prefix):
    '''Load and extract protein set from pLink'''
    xl_results = pd.read_csv(f"{pLink_res_prefix}.filtered_cross-linked_spectra.csv")
    mono_results = pd.read_csv(f"{pLink_res_prefix}.filtered_mono-linked_spectra.csv")
    loop_results = pd.read_csv(f"{pLink_res_prefix}.filtered_loop-linked_spectra.csv")
    regular_results = pd.read_csv(f"{pLink_res_prefix}.filtered_regular_spectra.csv")

    xl_protein_set = get_protein_set_from_list(list(set(
            xl_results["Proteins"].apply(rmv_xl_site).apply(split_xl_protein).tolist())))
    mono_protein_set = get_protein_set_from_list(list(set(
            mono_results["Proteins"].apply(rmv_xl_site).tolist())))
    loop_protein_set = get_protein_set_from_list(list(set(
            loop_results["Proteins"].apply(rmv_xl_site).tolist())))
    regular_protein_set = get_protein_set_from_list(list(set(
            regular_results["Proteins"].tolist())))
    
    protein_set = set()
    protein_set.update(xl_protein_set)
    protein_set.update(mono_protein_set)
    protein_set.update(loop_protein_set)
    protein_set.update(regular_protein_set)
    return protein_set