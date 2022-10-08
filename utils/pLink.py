import pandas as pd
import re
import functools
from theoretical_peaks.AAMass import aamass
from .xl_utils import *


def get_unfiltered_CSM_from_pL(res_path, keep_columns=['Title', 'Modifications', 'Score', "E-value"]):
    '''
    Load all unfiltered results
    '''

    """
    The columns are like:
    'Order', 'Title', 'Charge', 'Precursor_MH', 'Peptide_Type', 'Peptide',
    'Peptide_MH', 'Modifications', 'Refined_Score', 'SVM_Score', 'Score',
    'E-value', 'Precursor_Mass_Error(Da)', 'Precursor_Mass_Error(ppm)',
    'Target_Decoy', 'Q-value', 'Proteins', 'Protein_Type', 'FileID',
    'isComplexSatisfied', 'isFilterIn'
    """

    spectra_file = pd.read_csv(res_path)
    spectra_file = spectra_file[keep_columns].fillna("")
    return spectra_file


def get_CSM_from_pL(res_path, keep_columns=None,
    rename=False, sort_modifications=False, sort_alpha_beta=False, calc_exp_mz=False,
    filter_1_scan=False):
    '''
    Load results from _spectra.csv text file from pLink results as pd.DataFrame
    '''

    """
    The columns are like:
    'Order', 'Title', 'Charge', 'Precursor_Mass', 'Peptide', 'Peptide_Type',
    'Linker', 'Peptide_Mass', 'Modifications', 'Evalue', 'Score',
    'Precursor_Mass_Error(Da)', 'Precursor_Mass_Error(ppm)', 'Proteins',
    'Protein_Type', 'FileID', 'LabelID', 'Alpha_Matched', 'Beta_Matched',
    'Alpha_Evalue', 'Beta_Evalue', 'Alpha_Seq_Coverage',
    'Beta_Seq_Coverage'
    """

    spectra_file = pd.read_csv(res_path).fillna("")
    if keep_columns != None:
        spectra_file = spectra_file[keep_columns]

    # Keep 1 PSM for 1 scan
    # from pzm
    if filter_1_scan:
        spectra_file['raw_scan'] = spectra_file['Title'].apply(lambda x:'.'.join(x.split('.')[:-4]))

        ls_rawscan = []
        ls_1scan_filter = []

        for i in range(spectra_file.shape[0]):
            ser = spectra_file.iloc[i]
            raw_scan = ser['raw_scan']
            if raw_scan in ls_rawscan:
                ls_1scan_filter.append(0)
            else:
                ls_rawscan.append(raw_scan)
                ls_1scan_filter.append(1)

        spectra_file['1scan_filter'] = ls_1scan_filter
        print('Droping duplicated PSMs in the same scan: ', len(spectra_file[spectra_file['1scan_filter']==0]))
        spectra_file = spectra_file[spectra_file['1scan_filter']==1]
        spectra_file.drop(columns=['raw_scan', '1scan_filter'], inplace=True)


    # Modification in pLink is not sorted by site order. Resort for comparison
    # between results from different search engine.
    if sort_modifications:
        def compare_mods(x, y):
            return int(re.search("(?<=\()\d+(?=\))", x).group()) - int(re.search("(?<=\()\d+(?=\))", y).group())

        def sort_modifications_by_site(modification):
            mods = re.split("\;", modification)
            return ";".join(sorted(mods, key=functools.cmp_to_key(compare_mods)))

        spectra_file["Modifications"] = [sort_modifications_by_site(mod) for mod in spectra_file["Modifications"]]
    
    # alpha peptide in pLink is the peptide with higher quality. Resort by mass
    # for comparison between results from different search engine.
    if sort_alpha_beta:
        spectra_file[["Peptide", "Modifications"]] = [sort_alpha_beta_order(row["Peptide"], row["Modifications"])
            for _, row in spectra_file.iterrows()]

    # "Precursor_Mass" in pLink is the mass of precursor. Computer exp m/z
    # for comparison between results from different search engine.
    if calc_exp_mz:
        spectra_file[["exp m/z"]] = list(map(
            lambda x, y: (x - aamass.mass_proton) / y + aamass.mass_proton,
            spectra_file["Peptide_Mass"], spectra_file["Charge"]))

    if rename:
        spectra_file = spectra_file.rename(columns={
            "Title": "title",
            "Peptide": "sequence",
            "Modifications": "modinfo",
            "Linker": "linker",
            "Proteins": "proteins"
            })
        spectra_file = spectra_file.set_index('title')   
    return spectra_file


def get_peptides_from_pL(res_path):
    '''
    Load results from _peptides.csv text file from pLink results as pd.DataFrame
    '''

    """
    The columns are like:
    Peptide_Order	Peptide	Peptide_Mass	Modifications	Proteins	Protein_Type					
	Spectrum_Order	Title	Charge	Precursor_Mass	Evalue	Score	Precursor_Mass_Error(Da)	Precursor_Mass_Error(ppm)	Alpha_Evalue	Beta_Evalue
    """

    df = pd.read_csv(
        res_path,
        names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"])
    
    peptides_file = df[~df["A"].isna()][["A", "B", "C", "D", "E", "F"]]
    peptides_file.set_axis(
        df.iloc[0][["A", "B", "C", "D", "E", "F"]],
        axis=1, inplace=True)
    peptides_file = peptides_file.drop(index=0).fillna("") 
    
    return peptides_file


def get_sites_from_pL(res_path):
    '''
    Load results from _sites.csv text file from pLink results as pd.DataFrame
    '''

    """
    The columns are like:
    Protein_Order	Protein	Unique_Peptide_Number	Spectrum_Number	Protein_Type									
	Spectrum_Order	Title	Charge	Precursor_Mass	Peptide	Peptide_Mass	Modifications	Evalue	Score	Precursor_Mass_Error(Da)	Precursor_Mass_Error(ppm)	Alpha_Evalue	Beta_Evalue
    """
    
    df = pd.read_csv(
        res_path,
        names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"])
    
    sites_file = df[~df["A"].isna()][["A", "B", "C", "D", "E"]]
    sites_file.set_axis(
        df.iloc[0][["A", "B", "C", "D", "E"]],
        axis=1, inplace=True)
    sites_file = sites_file.drop(index=0).fillna("") 
    
    return sites_file


def get_summary_from_pL(res_path):
    '''
    Read the .summary file in pLink results
    '''
    
    # Read data from .summary
    with open(res_path, "r") as f:        
        full_summary = f.read()
    
    # Spectrum View
    summary = re.split("Spectrum View", full_summary)[1]
    all_spec_sum = re.search("(?<=Spectra: )\d+", summary).group()
    valid_spec_sum = re.search("(?<=Spectra \(above threshold\): )\d+", summary).group()
    invalid_spec_sum = re.search("(?<=Spectra \(below threshold and decoy\): )\d+", summary).group()
    unknown_spec_sum = re.search("(?<=Spectra \(unknown\): )\d+", summary).group()
    xl_spec_sum = re.search("(?<=Cross-Linked Spectra \(above threshold\): )\d+", summary).group()
    loop_spec_sum = re.search("(?<=Loop-Linked Spectra \(above threshold\): )\d+", summary).group()
    mono_spec_sum = re.search("(?<=Mono-Linked Spectra \(above threshold\): )\d+", summary).group()
    regular_spec_sum = re.search("(?<=Regular Spectra \(above threshold\): )\d+", summary).group()
    
    # Peptide View
    summary = re.split("Peptide View \(above threshold\)", full_summary)[1]
    xl_pep_sum = re.search("(?<=Cross-Linked Peptides: )\d+", summary).group()
    loop_pep_sum = re.search("(?<=Loop-Linked Peptides: )\d+", summary).group()
    mono_pep_sum = re.search("(?<=Mono-Linked Peptides: )\d+", summary).group()
    regular_pep_sum = re.search("(?<=Regular Peptides: )\d+", summary).group()
    
    # Linked Sites View
    summary = re.split("Linked Sites View \(above threshold\)", full_summary)[1]
    xl_site_sum = re.search("(?<=Cross-Linked Sites: )\d+", summary).group()
    loop_site_sum = re.search("(?<=Loop-Linked Sites: )\d+", summary).group()

    # Speed
    time_cost = re.search("(?<=Time Cost: )\d+\.*\d*(?=s\n)", summary).group()

    return {
        "all_spec_sum": int(all_spec_sum),
        "valid_spec_sum": int(valid_spec_sum),
        "invalid_spec_sum": int(invalid_spec_sum),
        "unknown_spec_sum": int(unknown_spec_sum),
        "xl_spec_sum": int(xl_spec_sum),
        "loop_spec_sum": int(loop_spec_sum),
        "mono_spec_sum": int(mono_spec_sum),
        "regular_spec_sum": int(regular_spec_sum),
        "xl_pep_sum": int(xl_pep_sum),
        "loop_pep_sum": int(loop_pep_sum),
        "mono_pep_sum": int(mono_pep_sum),
        "regular_pep_sum": int(regular_pep_sum),
        "xl_site_sum": int(xl_site_sum),
        "loop_site_sum": int(loop_site_sum),
        "time": float(time_cost),
        "full_summary": full_summary
    }
    