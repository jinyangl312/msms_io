import pandas as pd


def load_spectra_pQ(res_path):
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
    # rmv redudant label line in pQuant results
    spectra_file = spectra_file[spectra_file["Name_MS2"] != "Name_MS2"]
    spectra_file["Ratio_Sample2/Sample1"] = spectra_file["Ratio_Sample2/Sample1"].swifter.apply(
        float)

    return spectra_file
