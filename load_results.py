import pandas as pd


def get_seq_mod_from_pF_res(res_path):
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
    spectra_file = spectra_file[['File_Name', 'Sequence', 'Modification', 'Proteins']]
    spectra_file = spectra_file.rename(columns={
        "File_Name": "title", 
        "Sequence": "sequence", 
        "Modification": "modinfo",
        "Proteins": "proteins"
        }).set_index('title').fillna("")
    return spectra_file


def get_seq_mod_from_pL_res(res_path):
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
    spectra_file = spectra_file[['Title', 'Linker', 'Peptide', 'Modifications', 'Proteins']].fillna("")
    spectra_file = spectra_file.rename(columns={
        "Title": "title",
        "Peptide": "sequence",
        "Modifications": "modinfo",
        "Linker": "linker",
        "Proteins": "proteins"
        }).set_index('title').fillna("")
    return spectra_file
