import pandas as pd
import re


def get_seq_mod_from_pF_res(res_path, keep_columns=['File_Name', 'Sequence', 'Modification', 'Proteins'], rename=True):
    '''
    Load columns from .spectra text file from pFind results as pd.DataFrame
    '''
    
    '''
    The columns are like:
    File_Name	Scan_No	Exp.MH+	Charge	Q-value	Sequence	Calc.MH+	Mass_Shift(Exp.-Calc.)
    Raw_Score	Final_Score	Modification	Specificity	Proteins	Positions	Label
    Target/Decoy	Miss.Clv.Sites	Avg.Frag.Mass.Shift	Others
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


def get_seq_mod_from_pL_res(res_path, keep_columns=['Title', 'Linker', 'Peptide', 'Modifications', 'Proteins'], rename=True):
    '''
    Load columns from _spectra.csv text file from pLink results as pd.DataFrame
    '''

    """
    The columns are like:
    Order, Title, Charge, Precursor_Mass, Peptide, Peptide_Type, Linker, Peptide_Mass, Modifications,
    Evalue, Score, Precursor_Mass_Error(Da), Precursor_Mass_Error(ppm), Proteins, Protein_Type, 
    FileID, LabelID, Alpha_Matched, Beta_Matched, Alpha_Evalue, Beta_Evalue
    """

    spectra_file = pd.read_csv(res_path)
    spectra_file = spectra_file[keep_columns].fillna("")
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

    
def get_all_res_from_pL(res_path, keep_columns=['Title', 'Modifications', 'Score', "E-value"]):
    '''
    Load all unfiltered results
    '''

    """
    The columns are like:
    Order,Title,Charge,Precursor_MH,
    Peptide_Type,Peptide,Peptide_MH,
    Modifications,Refined_Score,SVM_Score,Score,E-value,
    Precursor_Mass_Error(Da),Precursor_Mass_Error(ppm),
    Target_Decoy,Q-value,Proteins,Protein_Type,
    FileID,isComplexSatisfied,isFilterIn
    """

    spectra_file = pd.read_csv(res_path)
    spectra_file = spectra_file[keep_columns].fillna("")
    return spectra_file


def get_summary_from_pF_res(res_path):
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

    
def get_summary_from_pL_res(res_path):
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