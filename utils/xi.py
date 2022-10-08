import pandas as pd
from .xl_utils import *


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