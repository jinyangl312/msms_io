import pandas as pd
import re

def load_spectra_from_pF(res_path, filtered=True, keep_target=True):
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
    if not filtered:
        spectra_file = spectra_file[spectra_file["Q-value"].apply(lambda x: x < 1)]
    if keep_target:
        spectra_file = spectra_file[spectra_file["Target/Decoy"].apply(lambda x: x == "target")]

    return spectra_file


def load_protein_from_pF(res_path, keep_filtered_only=True, keep_subset=True):
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
        names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S"])
    
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
        proteins = df[~df["B"].isna()][["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]]
        proteins.set_axis(
            df.iloc[0][["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]],
            axis=1, inplace=True)
        proteins = proteins.drop(index=0)
    else:
        proteins = df[~df["A"].isna()][["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]]
        proteins.set_axis(
            df.iloc[0][["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]],
            axis=1, inplace=True)
        proteins = proteins.drop(index=0)
    
    return proteins


def load_summary_from_pF(res_path):
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
