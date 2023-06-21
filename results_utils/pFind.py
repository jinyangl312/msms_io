import pandas as pd
import re
from .functions import *
import functools
import multiprocessing
import os
from pathlib import Path


def convert_pF_mod_to_pL_mod(pFind_mod):
    # pF: 19,Oxidation[M];
    # pL: Oxidation[M](10)
    if len(pFind_mod.split(";")) == 1:
        return ""

    modification_list = list(map(
        lambda x: "{}({})".format(
            re.search("(?<=\,).*", x).group(),
            re.search("\d+", x).group()
        ), pFind_mod.strip(";").split(";")))

    return ";".join(modification_list)


def load_spectra_pF(res_path, sort_modifications=False,
                    evaluation_scaffold=False, is_filtered=True, keep_target=True):
    '''
    Load columns from .spectra text file from pFind results as pd.DataFrame
    '''

    '''
    The columns are like:
    'File_Name', 'Scan_No', 'Exp.MH+', 'Charge', 'Q-value', 'Sequence',
       'Calc.MH+', 'Mass_Shift(Exp.-Calc.)', 'Raw_Score', 'Final_Score',
       'Modification', 'Specificity', 'Proteins', 'Positions', 'Label',
       'Target/Decoy', 'Miss.Clv.Sites', 'Avg.Frag.Mass.Shift', 'Others'
    '''

    spectra_file = pd.read_csv(res_path, delimiter='\t')
    if not is_filtered:
        spectra_file = spectra_file[spectra_file["Q-value"].apply(
            lambda x: x < 1)]
    if keep_target:
        spectra_file = spectra_file[spectra_file["Target/Decoy"].apply(
            lambda x: x == "target")]

    # Modification in pLink is not sorted by site order. Resort for comparison
    # between results from different search engine.
    if sort_modifications:
        spectra_file["Modifications"] = spectra_file["Modification"].fillna("").apply(
            convert_pF_mod_to_pL_mod)

        def compare_mods(x, y):
            return int(re.search("(?<=\()\d+(?=\))", x).group()) - int(re.search("(?<=\()\d+(?=\))", y).group())

        def sort_modifications_by_site(modification):
            mods = re.split("\;", modification)
            return ";".join(sorted(mods, key=functools.cmp_to_key(compare_mods)))

        spectra_file["Modifications"] = [sort_modifications_by_site(
            mod) for mod in spectra_file["Modifications"]]

    if evaluation_scaffold:
        spectra_file["_spectrum_id"] = spectra_file["File_Name"]
        spectra_file["_scan_id"] = spectra_file["File_Name"].apply(
            lambda x: re.search("(?<=^).*?\.\d+\.\d+(?=\.\d+\.\d+\.dta\;?)", x).group())
        spectra_file["_scan_charge_id"] = spectra_file["File_Name"].apply(
            lambda x: re.search("(?<=^).*?\.\d+\.\d+\.\d+(?=\.\d+\.dta\;?)", x).group())
        if sort_modifications:
            spectra_file["_psm"] = spectra_file[["Sequence", "Charge", "Modifications", "_spectrum_id"]].apply(
                tuple, axis=1)
            spectra_file["_precursor"] = spectra_file[["Sequence", "Charge", "Modifications"]].apply(
                tuple, axis=1)
            spectra_file["_peptide"] = spectra_file[["Sequence", "Modifications"]].apply(
                tuple, axis=1)
        else:
            spectra_file["_psm"] = spectra_file[["Sequence", "Charge", "Modification", "_spectrum_id"]].apply(
                tuple, axis=1)
            spectra_file["_precursor"] = spectra_file[["Sequence", "Charge", "Modification"]].apply(
                tuple, axis=1)
            spectra_file["_peptide"] = spectra_file[["Sequence", "Modification"]].apply(
                tuple, axis=1)
        spectra_file["Peptide"] = spectra_file["Sequence"]  # comp with xl
        spectra_file["_sequence"] = spectra_file["Sequence"]

        spectra_file["_protein"] = spectra_file["Proteins"].apply(
            lambda x: ";".join(set(
                filter(not_empty, re.split("/", x)))))

    return spectra_file


def load_res_pF(res_path, sort_modifications=False,
                evaluation_scaffold=False, is_filtered=True, keep_target=True):

    if os.path.exists(f"{res_path}.pyres.npy"):
        return np.load(f"{res_path}.pyres.npy", allow_pickle=True)

    with open(res_path) as f:
        raw_text = f.read()

    # all_info = [process_line(x) for x in raw_text.strip('S\t').split('\nS\t')]
    with multiprocessing.Pool(int(os.cpu_count()/2)) as pool:
        all_info = pool.map(
            process_pFind_res_line, raw_text.strip('S\t').split('\nS\t'))

    np.save(f"{res_path}.pyres.npy", all_info)

    return all_info


def load_res_pF_dir(res_dir):
    res_dir = Path(res_dir)
    assert res_dir.is_dir()

    res_list = res_dir.glob('*.qry.res')
    all_info = []
    for res in res_list:
        all_info.extend(load_res_pF(res))
    return all_info


def load_protein_pF(res_path, keep_filtered_only=True, keep_subset=True):
    '''
    Load columns from .protein tsv file from pFind results as pd.DataFrame
    '''

    """
    The columns are like:
    ID	AC	Score	Q-Value	Coverage	No.Peptide	No.Sameset	No.Subset	Have_Distinct_Pep	Description
	    	ID	Sequence	Calc.MH+	Mass_Shift(Exp.-Calc.)	Raw_Score	Final_Score	Modification	Specificity	Proteins	Positions	Label	Target/Decoy	Miss.Clv.Sites	Avg.Frag.Mass.Shift	File_Name	Charge	Spec_Num
    """

    df = pd.read_csv(
        res_path,
        delimiter="\t",
        names=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S"])

    # rmv splitter lines
    df.set_index(keys="A", inplace=True)
    df = df.loc[:"----Summary----"].iloc[:-1]
    if keep_filtered_only:
        df = df.loc[:"----------------------------------------"].iloc[:-1]
    else:
        df = df.loc[:"----------------------------------------"].iloc[:-1].append(
            df.loc["----------------------------------------":].iloc[1:]
        )
    df = df.reset_index()

    if keep_subset:
        proteins = df[~df["B"].isna()][["A", "B", "C", "D", "E",
                                        "F", "G", "H", "I", "J"]]
        proteins.set_axis(
            df.iloc[0][["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]],
            axis=1, inplace=True)
        proteins = proteins.drop(index=0)
    else:
        proteins = df[~df["A"].isna()][["A", "B", "C", "D", "E",
                                        "F", "G", "H", "I", "J"]]
        proteins.set_axis(
            df.iloc[0][["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]],
            axis=1, inplace=True)
        proteins = proteins.drop(index=0)

    return proteins


def load_summary_pF(res_path):
    '''
    Read the .summary file in pFind results
    '''

    # Read data from .summary
    with open(res_path, "r") as f:
        summary = f.read()

    # TODO: restructure re expressions with (?<=)
    scans_num = re.search("Scans: \d*", summary)
    scans_num = re.search("\d*$", scans_num.group())
    peptides_num_sm = re.search("Peptides: \d*", summary)
    peptides_num_sm = re.search("\d*$", peptides_num_sm.group())
    sequence_num = re.search("Sequences: \d*", summary)
    sequence_num = re.search("\d*$", sequence_num.group())
    proteins_num = re.search("Proteins: \d*", summary)
    proteins_num = re.search("\d*$", proteins_num.group())
    proteins_group_num = re.search("Protein Groups: \d*", summary)
    proteins_group_num = re.search("\d*$", proteins_group_num.group())
    Avg_Seq_Cov = re.search("Avg_Seq_Cov: .*", summary)
    Avg_Seq_Cov = re.search("\d*\.\d*", Avg_Seq_Cov.group())
    ID_rate = re.search("Overall.*\%", summary)
    ID_rate = re.search("\d*\.\d*", ID_rate.group())
    specific_cleavage = re.search("Specific: .*", summary)
    specific_cleavage = re.search("\d*\.\d*", specific_cleavage.group())
    c_57_rate = re.search("Carbamidomethyl\[C\].*", summary)
    c_57_rate = re.search("\d*\.\d*", c_57_rate.group())
    Precursor_Mean = re.search("Precursor_Mean:.*", summary)
    Precursor_Mean = re.search("\d*\.\d", Precursor_Mean.group())

    return {
        "scans_num": int(scans_num.group()),
        "peptides_num_sm": int(peptides_num_sm.group()),
        "sequence_num": int(sequence_num.group()),
        "proteins_num": int(proteins_num.group()),
        "proteins_group_num": int(proteins_group_num.group()),
        "Avg_Seq_Cov": float(Avg_Seq_Cov.group()),
        "ID_rate": float(ID_rate.group()),
        "specific_cleavage": float(specific_cleavage.group()),
        "c_57_rate": float(c_57_rate.group()),
        "Precursor_Mean": float(Precursor_Mean.group())
    }, summary


if __name__ == '__main__':
    multiprocessing.freeze_support()
