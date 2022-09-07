import pandas as pd
import re
import functools
from theoretical_peaks.ion_calc import calc_pepmass
from theoretical_peaks.AAMass import aamass
from .utils.xl_utils import *


def get_PSM_from_pF(res_path, keep_columns=['File_Name', 'Sequence', 'Modification', 'Proteins'], rename=True):
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
    spectra_file = spectra_file[keep_columns].fillna("")
    if rename:
        spectra_file = spectra_file.rename(columns={
            "File_Name": "title", 
            "Sequence": "sequence", 
            "Modification": "modinfo",
            "Proteins": "proteins"
            })
        spectra_file = spectra_file.set_index('title')
    return spectra_file


def get_summary_from_pF(res_path):
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
    c_57_rate = re.findall("\d*\.\d*", c_57_rate.group())
    Precursor_Mean = re.search("Precursor_Mean:.*", summary)
    Precursor_Mean = re.findall("\d*\.\d", Precursor_Mean.group())

    return {
        "scans_num": int(scans_num),
        "peptides_num_sm": int(peptides_num_sm),
        "sequence_num": int(sequence_num),
        "proteins_num": int(proteins_num),
        "proteins_group_num": int(proteins_group_num),
        "Avg_Seq_Cov": float(Avg_Seq_Cov),
        "ID_rate": float(ID_rate),
        "specific_cleavage": float(specific_cleavage),
        "c_57_rate": float(c_57_rate),
        "Precursor_Mean": float(Precursor_Mean)
    }


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
    rename=True, sort_modifications=False, sort_alpha_beta=False, calc_exp_mz=False,
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
    

def convert_xi_seq_mod_site_format_to_pL_mixed_format(row):    
    xi_mod_to_pL_mod_dict = {
        "cm": "Carbamidomethyl",
        "ox": "Oxidation",
    }

    sequence1 = "".join(re.findall("[A-Z]", row["PepSeq1"]))
    sequence2 = "".join(re.findall("[A-Z]", row["PepSeq2"]))

    modification_list = []
    assert row["PepSeq1"][0] >= "A" and row["PepSeq2"][0] >= "A"
    for index1, aa in enumerate(re.findall("[A-Z][a-z]*", row["PepSeq1"])):
        if len(aa) == 1:
            continue    
        modification_list.append(f"{xi_mod_to_pL_mod_dict[aa[1:]]}[{aa[0]}]({index1+1})")
    for index2, aa in enumerate(re.findall("[A-Z][a-z]*", row["PepSeq2"])):
        if len(aa) == 1:
            continue    
        modification_list.append(f"{xi_mod_to_pL_mod_dict[aa[1:]]}[{aa[0]}]({index2+1+len(sequence1)+2+1})")

    return (
        f'{sequence1}({row["LinkPos1"]})-{sequence2}({row["LinkPos2"]})',
        ";".join(modification_list)
    )


def convert_xi_protein_site_format_to_pL_mixed_format(row):    
    xi_mod_to_pL_mod_dict = {
        "cm": "Carbamidomethyl",
        "ox": "Oxidation",
    }

    sequence1 = "".join(re.findall("[A-Z]", row["Peptide1"]))
    sequence2 = "".join(re.findall("[A-Z]", row["Peptide2"]))

    modification_list = []
    assert row["Peptide1"][0] >= "A" and row["Peptide2"][0] >= "A"
    for index1, aa in enumerate(re.findall("[A-Z][a-z]*", row["Peptide1"])):
        if len(aa) == 1:
            continue    
        modification_list.append(f"{xi_mod_to_pL_mod_dict[aa[1:]]}[{aa[0]}]({index1+1})")
    for index2, aa in enumerate(re.findall("[A-Z][a-z]*", row["Peptide2"])):
        if len(aa) == 1:
            continue    
        modification_list.append(f"{xi_mod_to_pL_mod_dict[aa[1:]]}[{aa[0]}]({index2+1+len(sequence1)+2+1})")

    return (
        f'{sequence1}({row["FromSite"]})-{sequence2}({row["ToSite"]})',
        ";".join(modification_list)
    )


def get_precursor_from_xi(res_path, format_modification_linksite=True, sort_alpha_beta=True,
replace_Isoleucine=True, add_scan_charge_id=True):
    '''
    Load results from _CSM_ csv file from xi search results as pd.DataFrame
    '''

    """
    The columns are like:
    'PSMID', 'run', 'scan', 'PeakListFileName', 'ScanId',
    'exp charge', 'exp m/z', 'exp mass', 'exp fractionalmass',
    'match charge', 'match mass', 'match fractionalmass',
    'Protein1', 'Description1', 'Decoy1',
    'Protein2', 'Description2', 'Decoy2',
    'PepSeq1', 'PepSeq2', 'PepPos1', 'PepPos2', 'PeptideLength1', 'PeptideLength2',
    'LinkPos1', 'LinkPos2', 'ProteinLinkPos1', 'ProteinLinkPos2',
    'Charge', 'Crosslinker', 'CrosslinkerModMass', 'PeptidesWithDoublets',
    'PeptidesWithStubs', 'minPepCoverage', 'Score',
    'isDecoy', 'isTT', 'isTD', 'isDD', 'fdrGroup', 'fdr', 'ifdr', 'PEP',
    'Unnamed: 43', 'PeptidePairFDR', 'Protein1FDR', 'Protein2FDR',
    'LinkFDR', 'PPIFDR', 'peptide pair id', 'link id', 'ppi id', 'info'
    """

    spectra_file = pd.read_csv(res_path).fillna("")
    if format_modification_linksite:
        spectra_file[["Peptide", "Modifications"]] = [convert_xi_seq_mod_site_format_to_pL_mixed_format(row) for _, row in spectra_file.iterrows()]
    
    # pLink has converted I to L. Replace I
    # for comparison between results from different search engine.
    if replace_Isoleucine:
        spectra_file["Peptide"] = [x.replace("I", "L") for x in spectra_file["Peptide"]]

    # alpha peptide in pLink is the peptide with higher quality. Resort by mass
    # for comparison between results from different search engine.
    if sort_alpha_beta:
        spectra_file[["Peptide", "Modifications"]] = [sort_alpha_beta_order(row["Peptide"], row["Modifications"])
            for _, row in spectra_file.iterrows()]

    if add_scan_charge_id:
        spectra_file[["scan_id"]] = list(map(
            lambda x, y: f"{x}.{y}.{y}",
            spectra_file["run"],
            spectra_file["scan"]))
        spectra_file[["scan_charge_id"]] = list(map(
            lambda x, y, z: f"{x}.{y}.{y}.{z}",
            spectra_file["run"],
            spectra_file["scan"],
            spectra_file["Charge"]))
        spectra_file[["PSM_id"]] = list(map(
            lambda x, y, z: f"{x}.{y}.{y}.{z}.0.dta",
            spectra_file["run"],
            spectra_file["scan"],
            spectra_file["Charge"]))

    return spectra_file


def get_peptides_from_xi(res_path, format_modification_linksite=True, sort_alpha_beta=True,
replace_Isoleucine=True):
    '''
    Load results from _PeptidePairs_ csv file from xi search results as pd.DataFrame
    '''

    """
    The columns are like:
    PeptidePairID	PSMIDs
    Protein1	Description1	Decoy1
    Protein2	Description2	Decoy2
    Peptide1	Peptide2	Start1	Start2	FromSite
    ToSite	FromProteinSite	ToProteinSite
    psmID	Crosslinker	Score	isDecoy
    isTT	isTD	isDD	fdrGroup	fdr	ifdr	PEP

    Protein1FDR	Protein2FDR	LinkFDR	PPIFDR
    
    link id	ppi id
    ... ...
    """

    peptides_file = pd.read_csv(res_path).fillna("")
    if format_modification_linksite:
        peptides_file[["Peptide", "Modifications"]] = [convert_xi_protein_site_format_to_pL_mixed_format(row) for _, row in peptides_file.iterrows()]
    
    # pLink has converted I to L. Replace I
    # for comparison between results from different search engine.
    if replace_Isoleucine:
        peptides_file["Peptide"] = [x.replace("I", "L") for x in peptides_file["Peptide"]]

    # alpha peptide in pLink is the peptide with higher quality. Resort by mass
    # for comparison between results from different search engine.
    if sort_alpha_beta:
        peptides_file[["Peptide", "Modifications"]] = [sort_alpha_beta_order(row["Peptide"], row["Modifications"])
            for _, row in peptides_file.iterrows()]

    return peptides_file


def get_PSM_from_pQ(res_path, keep_columns=None,
    rename=True, sort_modifications=False, sort_alpha_beta=False, calc_exp_mz=False):
    '''
    Load results from _spectra.csv text file from pLink results as pd.DataFrame
    '''

    """
    The columns are like:
    Name_MS2	Sequence	Modification	Group_Joint	Score_Identification
    Intensity_Precursor	Locus_Protein	Description_Protein	Number_Samples
    Ratio_Sample2/Sample1	Score_Interference	Similarity_IsotopicDIS_Sample1
    Similarity_IsotopicDIS_Sample2	Intensity_Sample1	Intensity_Sample2
    WidthEluting_Sample1	WidthEluting_Sample2	Flag_ProteinInfer
    """

    spectra_file = pd.read_csv(res_path, delimiter="\t").fillna("")
    if keep_columns != None:
        spectra_file = spectra_file[keep_columns]

    return spectra_file
