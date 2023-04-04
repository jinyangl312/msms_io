import pathlib
import pandas as pd
from .xl_utils import *
from .utils import *
from pyteomics.pyteomics import fasta
from fasta_scripts.SeqAn_pybinder import fmindex_encode, FMIndexDecoder
from evaluation_scripts.utils.utils_precursor import annote_precursor_table_fasta
from fasta_scripts.find_protein import *
import os
import swifter


def convert_xi_seq_mod_site_format_to_pL_mixed_format_xl(row):
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
        modification_list.append(
            f"{xi_mod_to_pL_mod_dict[aa[1:]]}[{aa[0]}]({index1+1})")
    for index2, aa in enumerate(re.findall("[A-Z][a-z]*", row["PepSeq2"])):
        if len(aa) == 1:
            continue
        modification_list.append(
            f"{xi_mod_to_pL_mod_dict[aa[1:]]}[{aa[0]}]({index2+1+len(sequence1)+2+1})")

    return (
        f'{sequence1}({row["LinkPos1"]})-{sequence2}({row["LinkPos2"]})',
        ";".join(modification_list)
    )


def convert_xi_seq_mod_site_format_to_pL_mixed_format_linear(row):
    xi_mod_to_pL_mod_dict = {
        "cm": "Carbamidomethyl",
        "ox": "Oxidation",
    }
    xi_mono_dict = {"dssonh", "dssonh2", 'leikerclvoh'}
    xi_loop_dict = {"dssoloop"}

    sequence1 = "".join(re.findall("[A-Z]", row["PepSeq1"]))

    modification_list = []
    mono_pos = -1
    loop_pos = -1
    assert row["PepSeq1"][0] >= "A"
    for index1, aa in enumerate(re.findall("[A-Z][a-z]*", row["PepSeq1"])):
        if len(aa) == 1:
            continue
        if aa[1:] in xi_mono_dict:
            mono_pos = index1+1
        elif aa[1:] in xi_loop_dict:
            loop_pos = index1+1
        elif 'leikerclvoh' in aa[1:]:
            mono_pos = index1+1
            modification_list.append(
                f"{xi_mod_to_pL_mod_dict[re.sub('leikerclvoh', '', aa[1:])]}[{aa[0]}]({index1+1})")
        else:
            modification_list.append(
                f"{xi_mod_to_pL_mod_dict[aa[1:]]}[{aa[0]}]({index1+1})")

    if mono_pos >= 0:
        return (
            f'{sequence1}({mono_pos})',
            ";".join(modification_list)
        )

    elif loop_pos >= 0:
        return (
            f'{sequence1}({loop_pos})(-1)',
            ";".join(modification_list)
        )

    else:
        return (
            f'{sequence1}',
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
        modification_list.append(
            f"{xi_mod_to_pL_mod_dict[aa[1:]]}[{aa[0]}]({index1+1})")
    for index2, aa in enumerate(re.findall("[A-Z][a-z]*", row["Peptide2"])):
        if len(aa) == 1:
            continue
        modification_list.append(
            f"{xi_mod_to_pL_mod_dict[aa[1:]]}[{aa[0]}]({index2+1+len(sequence1)+2+1})")

    return (
        f'{sequence1}({row["FromSite"]})-{sequence2}({row["ToSite"]})',
        ";".join(modification_list)
    )


def load_precursor_xi(res_path, evaluation_scaffold=False, is_xl=True,
                      format_modification_linksite=True, sort_alpha_beta=True,
                      replace_Isoleucine=True, fasta_path=None, syn_path=None):
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
    spectra_file = spectra_file[~spectra_file['isDecoy']]

    if format_modification_linksite:
        if is_xl:
            spectra_file[["Peptide", "Modifications"]] = [
                convert_xi_seq_mod_site_format_to_pL_mixed_format_xl(row) for _, row in spectra_file.iterrows()]
        else:
            spectra_file[["Peptide", "Modifications"]] = [
                convert_xi_seq_mod_site_format_to_pL_mixed_format_linear(row) for _, row in spectra_file.iterrows()]

    # pLink has converted I to L. Replace I
    # for comparison between results from different search engine.
    if replace_Isoleucine:
        spectra_file["Peptide"] = [
            x.replace("I", "L") for x in spectra_file["Peptide"]]

    # alpha peptide in pLink is the peptide with higher quality. Resort by mass
    # for comparison between results from different search engine.
    if sort_alpha_beta and is_xl:
        spectra_file[["Peptide", "Modifications"]] = [sort_alpha_beta_order(row["Peptide"], row["Modifications"])
                                                      for _, row in spectra_file.iterrows()]

    if evaluation_scaffold:
        spectra_file['run'] = spectra_file['PeakListFileName'].swifter.apply(
            lambda x: re.search(".*?(?=\.)", x).group())
        spectra_file[["_spectrum_id"]] = list(map(
            lambda x, y, z: f"{x}.{y}.{y}.{z}.0.dta",
            spectra_file["run"],
            spectra_file["scan"],
            spectra_file["Charge"]))
        spectra_file[["_scan_id"]] = list(map(
            lambda x, y: f"{x}.{y}.{y}",
            spectra_file["run"],
            spectra_file["scan"]))
        spectra_file[["_scan_charge_id"]] = list(map(
            lambda x, y, z: f"{x}.{y}.{y}.{z}",
            spectra_file["run"],
            spectra_file["scan"],
            spectra_file["Charge"]))

        spectra_file["_psm"] = spectra_file[["Peptide", "Charge", "Modifications", "_spectrum_id"]].swifter.apply(
            tuple, axis=1)
        spectra_file["_precursor"] = spectra_file[["Peptide", "Charge", "Modifications"]].swifter.apply(
            tuple, axis=1)
        spectra_file["_peptide"] = spectra_file[["Peptide", "Modifications"]].swifter.apply(
            tuple, axis=1)
        spectra_file["_sequence"] = spectra_file["Peptide"]

        if is_xl:
            spectra_file['Peptide_Type'] = 'Cross-Linked'
            spectra_file['Protein_Type'] = spectra_file['fdrGroup'].map(lambda x:
                                                                        'Intra-Protein' if 'Self ' in x else 'Inter-Protein')
            output_prefix = f"data/fmindex_tmp/{pathlib.Path(fasta_path).stem}/fasta.rev.fm"
            if not pathlib.Path(output_prefix).parent.exists():
                pathlib.Path(output_prefix).parent.mkdir(
                    exist_ok=True, parents=True)
                with fasta.read(fasta_path) as db:
                    transferred_data = [
                        (item.description, item.sequence.replace("I", "L")) for item in db]
                fasta.write(transferred_data, fasta_path.replace(
                    ".fasta", ".I2L.fasta"), file_mode="w")
                fasta.write_decoy_db(fasta_path.replace(".fasta", ".I2L.fasta"),
                                     fasta_path.replace(
                    ".fasta", ".I2L.rev.fasta"),
                    prefix="REV_",
                    mode="reverse",
                    file_mode="w",
                )

                print("Cleaning tmp files...")
                for file in pathlib.Path().rglob(f'{output_prefix}.*'):
                    os.remove(file)
                fmindex_encode(fasta_path.replace(
                    ".fasta", ".I2L.rev.fasta"), output_prefix)

            fm_decoder = FMIndexDecoder(output_prefix)
            spectra_file["site_inferred_from_fasta"] = spectra_file["Peptide"].swifter.apply(
                lambda x: find_xl_seq_in_fasta_db(x, fm_decoder))
            spectra_file["site_inferred_from_fasta"] = spectra_file["site_inferred_from_fasta"].map(
                sort_site_order)

            spectra_file["_rp"] = spectra_file["site_inferred_from_fasta"].swifter.apply(
                lambda x: ";".join(set(
                    filter(not_empty, re.split("/", x)))))
            spectra_file["_site"] = spectra_file["_rp"].swifter.apply(
                lambda x: ";".join(set(
                    filter(not_empty, re.split(";|(?<=\))\-", x)))))
            spectra_file["_ppi"] = spectra_file["_rp"].swifter.apply(
                rmv_xl_site)
            spectra_file["_protein"] = spectra_file["_rp"].swifter.apply(
                lambda x: ";".join(set(filter(not_empty,
                                              [rmv_xl_site(i.strip()) for i in re.split(";|(?<=\))\-", x)]))))
        else:
            spectra_file['Peptide_Type'] = spectra_file['Peptide'].map(
                lambda x: "Loop-Linked" if re.search("\(.+\)\(.+\)", x)
                else ("Mono-Linked" if re.search("\(.+\)", x) else "Common"))

        # scripts for synthesis peptides
        if is_xl and syn_path != None:
            output_prefix = f"data/fmindex_tmp/{pathlib.Path(syn_path).stem}/fasta.rev.fm"
            if not pathlib.Path(output_prefix).parent.exists():
                pathlib.Path(output_prefix).parent.mkdir(
                    exist_ok=True, parents=True)
                with fasta.read(syn_path) as db:
                    transferred_data = [
                        (item.description, item.sequence.replace("I", "L")) for item in db]
                fasta.write(transferred_data, syn_path.replace(
                    ".fasta", ".I2L.fasta"), file_mode="w")
                fasta.write_decoy_db(syn_path.replace(".fasta", ".I2L.fasta"),
                                     syn_path.replace(
                    ".fasta", ".I2L.rev.fasta"),
                    prefix="REV_",
                    mode="reverse",
                    file_mode="w",
                )

                print("Cleaning tmp files...")
                for file in pathlib.Path().rglob(f'{output_prefix}.*'):
                    os.remove(file)
                fmindex_encode(syn_path.replace(
                    ".fasta", ".I2L.rev.fasta"), output_prefix)

            fm_decoder = FMIndexDecoder(output_prefix)
            spectra_file["in_syn_db"] = spectra_file["Peptide"].swifter.apply(
                lambda x: False if find_xl_seq_in_fasta_db(x, fm_decoder) == "" else True)
    return spectra_file


def load_peptides_xi(res_path, format_modification_linksite=True, sort_alpha_beta=True,
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
        peptides_file[["Peptide", "Modifications"]] = [
            convert_xi_protein_site_format_to_pL_mixed_format(row) for _, row in peptides_file.iterrows()]

    # pLink has converted I to L. Replace I
    # for comparison between results from different search engine.
    if replace_Isoleucine:
        peptides_file["Peptide"] = [
            x.replace("I", "L") for x in peptides_file["Peptide"]]

    # alpha peptide in pLink is the peptide with higher quality. Resort by mass
    # for comparison between results from different search engine.
    if sort_alpha_beta:
        peptides_file[["Peptide", "Modifications"]] = [sort_alpha_beta_order(row["Peptide"], row["Modifications"])
                                                       for _, row in peptides_file.iterrows()]

    return peptides_file
