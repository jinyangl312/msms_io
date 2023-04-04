import pandas as pd
import re
import functools
from theoretical_peaks.AAMass import aamass
from .xl_utils import *
from .utils import *
import swifter
import pathlib
from pyteomics.pyteomics import fasta
from fasta_scripts.SeqAn_pybinder import fmindex_encode, FMIndexDecoder
import os
from fasta_scripts.find_protein import *


def load_spectra_pL(res_path, evaluation_scaffold=False,
                    sort_modifications=False, sort_alpha_beta=False, calc_exp_mz=False,
                    filter_1_scan=False, keep_target=False,
                    keep_filtered=True, split_results=False,
                    syn_path=None):
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

    """
    For all.csv, the columns are like:
    'Order', 'Title', 'Charge', 'Precursor_MH', 'Peptide_Type', 'Peptide',
    'Peptide_MH', 'Modifications', 'Refined_Score', 'SVM_Score', 'Score',
    'E-value', 'Precursor_Mass_Error(Da)', 'Precursor_Mass_Error(ppm)',
    'Target_Decoy', 'Q-value', 'Proteins', 'Protein_Type', 'FileID',
    'isComplexSatisfied', 'isFilterIn'
    """

    spectra_file = pd.read_csv(res_path).fillna("")
    if len(spectra_file) == 0:
        return spectra_file

    if keep_target:
        try:
            spectra_file = spectra_file[spectra_file["Target_Decoy"] == 2]
        except:
            pass

    if "isFilterIn" in spectra_file.columns and keep_filtered:
        try:
            spectra_file = spectra_file[spectra_file["isFilterIn"] == 1]
        except:
            pass

    # Keep 1 PSM for 1 scan
    # from pzm
    if filter_1_scan:
        spectra_file['raw_scan'] = spectra_file['Title'].swifter.apply(
            lambda x: '.'.join(x.split('.')[:-4]))

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
        print('Droping duplicated PSMs in the same scan: ', len(
            spectra_file[spectra_file['1scan_filter'] == 0]))
        spectra_file = spectra_file[spectra_file['1scan_filter'] == 1]
        spectra_file.drop(columns=['raw_scan', '1scan_filter'], inplace=True)

    # Modification in pLink is not sorted by site order. Resort for comparison
    # between results from different search engine.
    if sort_modifications:
        def compare_mods(x, y):
            return int(re.search("(?<=\()\d+(?=\))", x).group()) - int(re.search("(?<=\()\d+(?=\))", y).group())

        def sort_modifications_by_site(modification):
            mods = re.split("\;", modification)
            return ";".join(sorted(mods, key=functools.cmp_to_key(compare_mods)))

        spectra_file["Modifications"] = [sort_modifications_by_site(
            mod) for mod in spectra_file["Modifications"]]

    # alpha peptide in pLink is the peptide with higher quality. Resort by mass
    # for comparison between results from different search engine.
    if sort_alpha_beta:
        spectra_file["Peptide", "Modifications"] = spectra_file[["Peptide", "Modifications"]].swifter.apply(
            lambda x: sort_alpha_beta_order(*x), axis=1)

    if evaluation_scaffold:
        spectra_file["_spectrum_id"] = spectra_file["Title"]
        spectra_file["_scan_id"] = spectra_file["Title"].swifter.apply(
            lambda x: re.search("(?<=^).*?\.\d+\.\d+(?=\.\d+\.\d+\.dta\;?)", x).group())
        spectra_file["_scan_charge_id"] = spectra_file["Title"].swifter.apply(
            lambda x: re.search("(?<=^).*?\.\d+\.\d+\.\d+(?=\.\d+\.dta\;?)", x).group())

        spectra_file["_psm"] = spectra_file[["Peptide", "Charge", "Modifications", "_spectrum_id"]].swifter.apply(
            tuple, axis=1)
        spectra_file["_precursor"] = spectra_file[["Peptide", "Charge", "Modifications"]].swifter.apply(
            tuple, axis=1)
        spectra_file["_peptide"] = spectra_file[["Peptide", "Modifications"]].swifter.apply(
            tuple, axis=1)
        spectra_file["_sequence"] = spectra_file["Peptide"]

        if spectra_file["Peptide_Type"][0] == 'Cross-Linked':
            spectra_file["Proteins"] = spectra_file["Proteins"].map(
                sort_site_order)

        spectra_file["_rp"] = spectra_file["Proteins"].swifter.apply(
            lambda x: ";".join(set(
                filter(not_empty, re.split("/", x)))))
        spectra_file["_site"] = spectra_file["_rp"].swifter.apply(
            lambda x: ";".join(set(
                filter(not_empty, re.split(";|(?<=\))\-", x)))))
        spectra_file["_ppi"] = spectra_file["_rp"].swifter.apply(rmv_xl_site)
        spectra_file["_protein"] = spectra_file["_rp"].swifter.apply(
            lambda x: ";".join(set(filter(not_empty,
                                          [rmv_xl_site(i.strip()) for i in re.split(";|(?<=\))\-", x)]))))

    # "Precursor_Mass" in pLink is the mass of precursor. Computer exp m/z
    # for comparison between results from different search engine.
    if calc_exp_mz:
        spectra_file[["exp m/z"]] = list(map(
            lambda x, y: (x - aamass.mass_proton) / y + aamass.mass_proton,
            spectra_file["Peptide_Mass"], spectra_file["Charge"]))

    # scripts for synthesis peptides
    if syn_path != None:
        output_prefix = f"data/fmindex_tmp/{pathlib.Path(syn_path).stem}/syn/fasta.fm"
        if not pathlib.Path(output_prefix).parent.exists():
            pathlib.Path(output_prefix).parent.mkdir(
                exist_ok=True, parents=True)
            with fasta.read(syn_path) as db:
                transferred_data = [
                    (item.description, item.sequence.replace("I", "L")) for item in db]
            fasta.write(transferred_data, syn_path.replace(
                ".fasta", ".I2L.fasta"), file_mode="w")

            print("Cleaning tmp files...")
            for file in pathlib.Path().rglob(f'{output_prefix}.*'):
                os.remove(file)
            fmindex_encode(syn_path.replace(
                ".fasta", ".I2L.fasta"), output_prefix)

        fm_decoder = FMIndexDecoder(output_prefix)
        if spectra_file["Peptide_Type"][0] == 'Cross-Linked':
            spectra_file["in_syn_db"] = spectra_file["Peptide"].swifter.apply(
                lambda x: is_xl_seq_in_same_syn_group(x, fm_decoder))
                
        # Trap DB!
        # if spectra_file["Peptide_Type"][0] == 'Cross-Linked':
        #     spectra_file["in_syn_db"] = spectra_file["Peptide"].swifter.apply(
        #         lambda x: False if find_xl_seq_in_fasta_db(x, fm_decoder) == "" else True)
        # else:
        #     spectra_file["in_syn_db"] = spectra_file["Peptide"].swifter.apply(
        #         lambda x: False if len(find_seq_in_fasta_db(rmv_xl_site(x), 0, fm_decoder)) == 0 else True)

    if "Peptide_Type" in spectra_file.columns and split_results:
        return spectra_file[spectra_file["Peptide_Type"] == 0], \
            spectra_file[spectra_file["Peptide_Type"] == 1], \
            spectra_file[spectra_file["Peptide_Type"] == 2], \
            spectra_file[spectra_file["Peptide_Type"] == 3], \
            spectra_file
    else:
        return spectra_file


def load_peptides_pL(res_path):
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
        names=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"])

    peptides_file = df[~df["A"].isna()][["A", "B", "C", "D", "E", "F"]]
    peptides_file.set_axis(
        df.iloc[0][["A", "B", "C", "D", "E", "F"]],
        axis=1, inplace=True)
    peptides_file = peptides_file.drop(index=0).fillna("")

    return peptides_file


def load_sites_pL(res_path):
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
        names=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"])

    sites_file = df[~df["A"].isna()][["A", "B", "C", "D", "E"]]
    sites_file.set_axis(
        df.iloc[0][["A", "B", "C", "D", "E"]],
        axis=1, inplace=True)
    sites_file = sites_file.drop(index=0).fillna("")

    return sites_file


def load_summary_pL(res_path):
    '''
    Read the .summary file in pLink results
    '''

    # Read data from .summary
    with open(res_path, "r") as f:
        full_summary = f.read()

    # Spectrum View
    summary = re.split("Spectrum View", full_summary)[1]
    all_spec_sum = re.search("(?<=Spectra: )\d+", summary).group()
    valid_spec_sum = re.search(
        "(?<=Spectra \(above threshold\): )\d+", summary).group()
    invalid_spec_sum = re.search(
        "(?<=Spectra \(below threshold and decoy\): )\d+", summary).group()
    unknown_spec_sum = re.search(
        "(?<=Spectra \(unknown\): )\d+", summary).group()
    xl_spec_sum = re.search(
        "(?<=Cross-Linked Spectra \(above threshold\): )\d+", summary).group()
    loop_spec_sum = re.search(
        "(?<=Loop-Linked Spectra \(above threshold\): )\d+", summary).group()
    mono_spec_sum = re.search(
        "(?<=Mono-Linked Spectra \(above threshold\): )\d+", summary).group()
    regular_spec_sum = re.search(
        "(?<=Regular Spectra \(above threshold\): )\d+", summary).group()

    # Peptide View
    summary = re.split("Peptide View \(above threshold\)", full_summary)[1]
    xl_pep_sum = re.search("(?<=Cross-Linked Peptides: )\d+", summary).group()
    loop_pep_sum = re.search("(?<=Loop-Linked Peptides: )\d+", summary).group()
    mono_pep_sum = re.search("(?<=Mono-Linked Peptides: )\d+", summary).group()
    regular_pep_sum = re.search("(?<=Regular Peptides: )\d+", summary).group()

    # Linked Sites View
    summary = re.split(
        "Linked Sites View \(above threshold\)", full_summary)[1]
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
