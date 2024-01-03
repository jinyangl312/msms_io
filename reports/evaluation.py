import pandas as pd
import re
import swifter

from msms_io.load_reports import load_spectra_pF, load_spectra_pL
from msms_io.tools.df_split_merge import *
from msms_io.evaluation.utils_PSM import add_pQuant_res_to_PSM, filter_PSM_by_protein_pF, add_pL_unfiltered_fields_to_PSM


def process_spectra_pL(
        spectra, is_filtered, 
        quant, pL_unfiltered_path=None, syn_path=None, keep_target=True, **args):
    if is_filtered:
        # Load filtered files
        if type(spectra) == list:
            PSM_res = pd.DataFrame()
            for _spectra in spectra:
                PSM_res = PSM_res.append(load_spectra_pL(
                    _spectra,
                    keep_target=True,
                    sort_modifications=True,
                    sort_alpha_beta=True,
                    evaluation_scaffold=True,
                    syn_path=syn_path,
                    **args,
                ), ignore_index=True)

            if pL_unfiltered_path == None:
                PSM_res['Q-value'] = -1
            else:
                PSM_res = add_pL_unfiltered_fields_to_PSM(
                    PSM_res, pL_unfiltered_path=pL_unfiltered_path)
        # Load one filtered file
        else:
            PSM_res = load_spectra_pL(
                spectra,
                keep_target=True,
                sort_modifications=True,
                sort_alpha_beta=True,
                evaluation_scaffold=True,
                syn_path=syn_path,
                **args,
            )
            if pL_unfiltered_path == None:
                PSM_res['Q-value'] = -1
            else:
                PSM_res = add_pL_unfiltered_fields_to_PSM(
                    PSM_res, pL_unfiltered_path=pL_unfiltered_path)
    # Load unfiltered file
    else:
        regular, mono, loop, xl, all = load_spectra_pL(
            spectra,
            keep_target=keep_target,
            sort_modifications=True,
            sort_alpha_beta=True,
            evaluation_scaffold=True,
            split_results=True,
            **args,
        )
        PSM_res = all

    if quant != None:
        PSM_res = add_pQuant_res_to_PSM(PSM_res, pQuant_path=quant)
    return PSM_res


def generate_res_table_pL(
        PSM_res, level, use_tuple,
        quant, suffix, mixed_peptide_type,
        keep_unique, split_column,
        rename_dict,
        syn_path=None,
        is_filtered=True,
        keep_min_cols=[],
        keep_max_cols=[], sum_cols=[],
        concat_str_cols=['_scan_charge_id']):
    # Set up header
    if level == "psm":
        by_column = (["_psm"] if use_tuple
                     else ["Peptide", "Charge", "Modifications", "_spectrum_id"]
                     ) + ["Peptide_Type"]
    elif level == "scan":
        by_column = ["_scan_id"]
    elif level == "spectrum":
        by_column = ["_spectrum_id"]
    elif level == "precursor":
        by_column = (["_precursor"] if use_tuple
                     else ["Peptide", "Charge", "Modifications"]
                     ) + ["Peptide_Type"]
    elif level == "peptide":
        by_column = (["_peptide"] if use_tuple
                     else ["Peptide", "Modifications"]
                     ) + ["Peptide_Type"]
    elif level == "sequence":
        by_column = (["_sequence"] if use_tuple
                     else ["Peptide"]
                     ) + ["Peptide_Type"]
    elif level == "rp":
        by_column = ["_rp"] + ["Protein_Type"]
    elif level == "ppi":
        by_column = ["_ppi"] + ["Protein_Type"]
    elif level == "protein":
        by_column = ["_protein"]
    else:
        assert False, 'Undefined level'

    # Keep unique if necessary
    if keep_unique:
        if level == "rp":
            PSM_res = PSM_res[PSM_res["_rp"].swifter.apply(
                lambda x: len(re.split(";", x))) == 1]
        elif level == "ppi" or level == "protein":
            PSM_res = PSM_res[PSM_res["_ppi"].swifter.apply(
                lambda x: len(re.split(";", x))) == 1]

    # split line with multiple ppis into lines
    if split_column:
        if level == "rp":
            PSM_res = split_column_into_columns(PSM_res, "_rp")
            PSM_res = PSM_res[PSM_res['_rp'].map(lambda x: 'REV_' not in x)]
        elif level == "ppi":
            PSM_res = split_column_into_columns(PSM_res, "_ppi")
            PSM_res = PSM_res[PSM_res['_ppi'].map(lambda x: 'REV_' not in x)]
        elif level == "protein":
            PSM_res = split_column_into_columns(PSM_res, "_protein")
            PSM_res = PSM_res[PSM_res['_protein'].map(
                lambda x: 'REV_' not in x)]

    def keep_columns(suffix):
        if level == "psm" or level == "spectrum":
            res_table = PSM_res.set_index(by_column)
        else:
            res_table = group_by_col(
                PSM_res,
                by=by_column,
                keep_min_cols=keep_min_cols,
                keep_max_cols=keep_max_cols,
                sum_cols=sum_cols,
                concat_str_cols=concat_str_cols,
            )
            # if quant != None:
            #     res_table = group_by_col(
            #         PSM_res,
            #         by=by_column,
            #         keep_min_cols=[],
            #         keep_max_cols=[
            #             'pQuant_evidence',
            #             'Score',
            #             'Q-value',
            #         ],
            #         concat_str_cols=[
            #             "Title",
            #             "Ratio_Sample2/Sample1",
            #         ]
            #     )
            # else:
            #     res_table = group_by_col(
            #         PSM_res,
            #         by=by_column,
            #         keep_min_cols=[],
            #         keep_max_cols=[
            #             'Score',
            #             'Q-value',
            #         ],
            #         concat_str_cols=[
            #             "Title",
            #         ]
            #     )

        res_table = res_table[rename_dict.keys()]
        res_table = res_table.rename(columns=rename_dict)

        # if quant != None:
        #     if level == "psm" and not is_filtered:
        #         res_table = res_table[["Title", "pQuant_evidence",
        #                                "Ratio_Sample2/Sample1", "Q-value", "Score", "Target_Decoy"]]
        #         res_table = res_table.rename(columns={
        #             "Title": f"title{suffix}",
        #             "pQuant_evidence": f"pQuant{suffix}",
        #             "Ratio_Sample2/Sample1": f"pQuant_ratio{suffix}",
        #             "Q-value": f"qvalue{suffix}",
        #             "Score": f"score{suffix}",
        #             "Target_Decoy": f"td{suffix}",
        #         })
        #     else:
        #         res_table = res_table[["Title", "pQuant_evidence",
        #                                "Ratio_Sample2/Sample1", "Q-value", "Score"]]
        #         res_table = res_table.rename(columns={
        #             "Title": f"title{suffix}",
        #             "pQuant_evidence": f"pQuant{suffix}",
        #             "Ratio_Sample2/Sample1": f"pQuant_ratio{suffix}",
        #             "Q-value": f"qvalue{suffix}",
        #             "Score": f"score{suffix}",
        #         })
        # else:
        #     if level == "psm" and not is_filtered:
        #         res_table = res_table[[
        #             "Title", "Q-value", "Score", "Target_Decoy"]]
        #         res_table = res_table.rename(columns={
        #             "Title": f"title{suffix}",
        #             "Q-value": f"qvalue{suffix}",
        #             "Score": f"score{suffix}",
        #             "Target_Decoy": f"td{suffix}",
        #         })
        #     else:
        #         res_table = res_table[["Title", "Q-value", "Score"]]
        #         res_table = res_table.rename(columns={
        #             "Title": f"title{suffix}",
        #             "Q-value": f"qvalue{suffix}",
        #             "Score": f"score{suffix}",
        #         })
        return res_table

    # Keep specific columns
    if mixed_peptide_type:
        return keep_columns(suffix)
    else:
        res_table = pd.DataFrame()
        for separate_suffix in [f"{suffix}_r", f"{suffix}_m", f"{suffix}_l", f"{suffix}_x"]:
            res_table = res_table.append(keep_columns(
                separate_suffix), ignore_index=True)
    return res_table


def process_spectra_pF(spectra, is_filtered, quant, protein_path, keep_target, keep_subset):
    PSM_res = load_spectra_pF(
        res_path=spectra,
        is_filtered=is_filtered,
        keep_target=keep_target,
        sort_modifications=True,
        evaluation_scaffold=True,
    )

    if quant != None:
        PSM_res = add_pQuant_res_to_PSM(PSM_res, pQuant_path=quant)

    if protein_path != None:
        PSM_res = filter_PSM_by_protein_pF(
            PSM_res,
            protein_path,
            keep_filtered_only=True,
            keep_subset=keep_subset,
        )
    return PSM_res


def generate_res_table_pF(
        PSM_res, level, use_tuple,
        quant, suffix, keep_unique,
        split_column, rename_dict,
        keep_min_cols=['Q-value'],
        keep_max_cols=['Final_Score'], sum_cols=[],
        concat_str_cols=['File_Name'],
):
    # Set up header
    if level == "psm":
        by_column = (["_psm"] if use_tuple
                     else ["Peptide", "Charge", "Modifications", "_spectrum_id"]
                     )
    elif level == "scan":
        by_column = ["_scan_id"]
    elif level == "spectrum":
        by_column = ["_spectrum_id"]
    elif level == "precursor":
        by_column = (["_precursor"] if use_tuple
                     else ["Peptide", "Charge", "Modifications"]
                     )
    elif level == "peptide":
        by_column = (["_peptide"] if use_tuple
                     else ["Peptide", "Modifications"]
                     )
    elif level == "sequence":
        by_column = (["_sequence"] if use_tuple
                     else ["Peptide"]
                     )
    elif level == "protein":
        by_column = ["_protein"]
    else:
        assert False, 'Undefined level'

    # Keep unique if necessary
    if keep_unique:
        if level == "protein":
            PSM_res = PSM_res[PSM_res["_protein"].swifter.apply(
                lambda x: len(re.split(";", x))) == 1]

    # split line with multiple ppis into lines
    if split_column:
        if level == "protein":
            PSM_res = split_column_into_columns(PSM_res, "_protein")
            PSM_res = PSM_res[PSM_res['_protein'].swifter.apply(
                lambda x: re.search("REV_", x) == None)]

    # Group to target level
    if level == "psm" or level == "spectrum":
        res_table = PSM_res.set_index(by_column)
    else:
        res_table = group_by_col(
            PSM_res,
            by=by_column,
            keep_min_cols=keep_min_cols,
            keep_max_cols=keep_max_cols,
            concat_str_cols=concat_str_cols,
            sum_cols=sum_cols,
        )

    # Keep specific columns
    res_table = res_table[rename_dict.keys()]
    res_table = res_table.rename(columns=rename_dict)
    return res_table
