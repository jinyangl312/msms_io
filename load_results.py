import pandas as pd


def get_seq_mod_from_pF_res(res_path):
    '''
    Load columns from pFind results
    '''
    
    spectra_file = pd.read_csv(res_path, delimiter='\t')
    spectra_file = spectra_file[['File_Name', 'Sequence', 'Modification', 'Proteins']]
    '''
    File_Name	Scan_No	Exp.MH+	Charge	Q-value	Sequence	Calc.MH+	Mass_Shift(Exp.-Calc.)
    Raw_Score	Final_Score	Modification	Specificity	Proteins	Positions	Label
    Target/Decoy	Miss.Clv.Sites	Avg.Frag.Mass.Shift	Others
    '''

    spectra_file = spectra_file.rename(columns={
        "File_Name": "title", 
        "Sequence": "sequence", 
        "Modification": "modinfo",
        "Proteins": "proteins"
        }).set_index('title').fillna("")
    return spectra_file


def get_seq_mod_from_pL_res(res_path): # title, sequence, modinfo
    '''
    Load columns from pLink results
    
    '''
    spectra_file = pd.read_csv(res_path)
    spectra_file = spectra_file[['Title', 'Linker', 'Peptide', 'Modifications', 'Proteins']].fillna("")
    """
    Order, Title, Charge, Precursor_Mass, Peptide, Peptide_Type, Linker, Peptide_Mass, Modifications,
    Evalue, Score, Precursor_Mass_Error(Da), Precursor_Mass_Error(ppm), Proteins, Protein_Type, 
    FileID, LabelID, Alpha_Matched, Beta_Matched, Alpha_Evalue, Beta_Evalue
    """

    spectra_file['Linker'].replace(['Leiker_clv'], 'Leiker', inplace=True)
    spectra_file = spectra_file.rename(columns={
        "Title": "title",
        "Peptide": "sequence",
        "Modifications": "modinfo",
        "Linker": "linker",
        "Proteins": "proteins"
        }).set_index('title').fillna("")
    return spectra_file
