
from logging import exception
from .ion_calc import *
from .AAMass import aamass
from .xlink import xlmass
import numpy as np
import re


def format_pL_modinfo(modinfo):
    '''
    Format modinfo in pLink into a mod list
    '''

    mod_list = []
    if not modinfo == "":
        for mod_str in modinfo.split(';'):
            modname = mod_str.split('(')[0]
            site = int(mod_str.split('(')[1][:-1])
            mod_list.append((site-1, modname))
    return mod_list


def format_pL_modinfo_xl(mods_str, len_pep_1):
    '''
    Format mods_str in pLink into a mod list
    len_pep_1: length of the first peptide in xl results
    '''

    if mods_str == "":
        return [], []

    mp_mod1 = []
    mp_mod2 = []
    for mod_str in mods_str.split(';'):
        modname = mod_str.split('(')[0]
        site = int(mod_str.split('(')[1][:-1])
        if site > len_pep_1:
            mp_mod2.append((site-4-len_pep_1, modname))
        else:
            mp_mod1.append((site-1, modname))

    return mp_mod1, mp_mod2


def format_pF_modinfo(modinfo):
    items = modinfo.split(";")
    modlist = []
    for mod in items:
        if mod != '':
            site, modname = mod.split(",")
            site = int(site)
            modlist.append((site, modname))
    return modlist


def format_linker_mass_xl(seq_length, linker_mass, linksite):
    '''Calc xlinker mass info and format it with pepmass into a list'''
    xlmassinfo = [0]*(seq_length+2)
    xlmassinfo[linksite] = linker_mass
    return xlmassinfo


def format_linker_mass(xl_seq_str, linker, mode):
    '''Get xlinker mass info and format it into a list'''

    line = re.split("\(|\)|\-", xl_seq_str)
    sequence, linksite = line[0], int(line[1])
    xlmassinfo = [0]*(len(sequence)+2)

    if mode == "mono":
        xlmassinfo[linksite] = xlmass.get_mono_mass(linker)
    elif mode == "loop":
        xlmassinfo[linksite] = xlmass.get_linker_mass(linker)
    else:
        raise exception

    return xlmassinfo


def cal_theoretical_b_y_peaks(sequence, modinfo, xlmassinfo=None):
    '''
    Return mz for b/y ions from sequence and modinfo info
    sequence: ABCDE
    modinfo: a string like "26,Carbamidomethyl[C];" or a list like [(26, Carbamidomethyl[C])]
    '''

    theoretical_peaks = {}
    bions, pepmass = calc_b_ions(sequence, modinfo, xlmassinfo)
    theoretical_peaks['b'] = bions
    theoretical_peaks['y'] = calc_y_from_b(bions, pepmass)
    #theoretical_peaks['b-ModLoss'] = calc_ion_modloss(bions, sequence, modinfo, N_term = True)
    #theoretical_peaks['y-ModLoss'] = calc_ion_modloss(theoretical_peaks['y'], sequence, modinfo, N_term = False)
    theoretical_peaks['y'].reverse()
    # theoretical_peaks['y-ModLoss'].reverse()
    # , theoretical_peaks['b-ModLoss'], theoretical_peaks['y-ModLoss']
    return theoretical_peaks['b'], theoretical_peaks['y']


def cal_theoretical_b_y_peaks_xl(sequence, modinfo, linker):
    '''
    Return mz for ab/ay/bb/by ions from sequence and modinfo info from pL
    sequence: QEDFYPFLKDNR(9)-VKYVTEGMR(2)
    modinfo: Oxidation[M](23)
    linker: DSSO
    '''

    # Split xl into two peptides
    # ['LCVLHEKTPVSEK', '7', '', 'CASIQKFGER', '6', '']
    line = re.split("\(|\)|\-", sequence)
    sequence1, sequence2, linksite1, linksite2 = line[0], line[3], int(
        line[1]), int(line[4])  # starts from 1

    # Split modinfo into two lists to the two peptides
    modinfo1, modinfo2 = format_pL_modinfo_xl(
        modinfo, len(sequence.split('-')[0].split('(')[0]))

    # Calc xl mass info and deliver to two peptides as a large modification
    if xlmass.is_cleavable(linker):
        xlmassinfo1_S = format_linker_mass_xl(
            len(sequence1), xlmass.get_short_arm_mass(linker), linksite1)
        xlmassinfo2_S = format_linker_mass_xl(
            len(sequence2), xlmass.get_short_arm_mass(linker), linksite2)
        xlmassinfo1_L = format_linker_mass_xl(
            len(sequence1), xlmass.get_long_arm_mass(linker), linksite1)
        xlmassinfo2_L = format_linker_mass_xl(
            len(sequence2), xlmass.get_long_arm_mass(linker), linksite2)
        xlmassinfo1 = format_linker_mass_xl(len(sequence1), xlmass.get_linker_mass(
            linker) + calc_pepmass(sequence2, modinfo2), linksite1)
        xlmassinfo2 = format_linker_mass_xl(len(sequence2), xlmass.get_linker_mass(
            linker) + calc_pepmass(sequence1, modinfo1), linksite2)

        # Now the xl can be regarded as two single peptides with modifications. Calculate theoretical peaks.
        return \
            cal_theoretical_b_y_peaks_cleavable_arms(sequence1, modinfo1, xlmassinfo1_S, linksite1) + \
            cal_theoretical_b_y_peaks_cleavable_arms(sequence2, modinfo2, xlmassinfo2_S, linksite2) + \
            cal_theoretical_b_y_peaks_cleavable_arms(sequence1, modinfo1, xlmassinfo1_L, linksite1) + \
            cal_theoretical_b_y_peaks_cleavable_arms(sequence2, modinfo2, xlmassinfo2_L, linksite2) + \
            cal_theoretical_b_y_peaks(sequence1, modinfo1, xlmassinfo1) + \
            cal_theoretical_b_y_peaks(sequence2, modinfo2, xlmassinfo2)
    else:
        xlmassinfo1 = format_linker_mass_xl(len(sequence1), xlmass.get_linker_mass(
            linker) + calc_pepmass(sequence2, modinfo2), linksite1)
        xlmassinfo2 = format_linker_mass_xl(len(sequence2), xlmass.get_linker_mass(
            linker) + calc_pepmass(sequence1, modinfo1), linksite2)

        return cal_theoretical_b_y_peaks(sequence1, modinfo1, xlmassinfo1) + \
            cal_theoretical_b_y_peaks(sequence2, modinfo2, xlmassinfo2)


def cal_theoretical_b_y_peaks_cleavable_arms(sequence, modinfo, xlmassinfo, linksite):
    '''
    Return mz for b/y ions from sequence and modinfo info.
    Changed mz of b ions left to the linked site and y ions right to the linked site to empty string.
    sequence: ABCDE(2)
    modinfo: a string like "26,Carbamidomethyl[C];" or a list like [(26, Carbamidomethyl[C])]
    '''

    theoretical_peaks = {}
    bions, pepmass = calc_b_ions(sequence, modinfo, xlmassinfo)
    theoretical_peaks['b'] = bions

    theoretical_peaks['y'] = calc_y_from_b(bions, pepmass)
    #theoretical_peaks['b-ModLoss'] = calc_ion_modloss(bions, sequence, modinfo, N_term = True)
    #theoretical_peaks['y-ModLoss'] = calc_ion_modloss(theoretical_peaks['y'], sequence, modinfo, N_term = False)
    theoretical_peaks['y'].reverse()
    # theoretical_peaks['y-ModLoss'].reverse()

    for i in range(1, linksite):
        theoretical_peaks['b'][i-1] = ""
    for i in range(1, len(sequence) + 1 - linksite):
        theoretical_peaks['y'][i-1] = ""

    # , theoretical_peaks['b-ModLoss'], theoretical_peaks['y-ModLoss']
    return theoretical_peaks['b'], theoretical_peaks['y']


def cal_theoretical_b_y_peaks_loop(sequence, modinfo, xlmassinfo=None):
    '''
    Return mz for b/y ions from sequence and modinfo info.
    Changed mz between loop-linked sites to empty string.
    sequence: ABCDE(1)(2)
    modinfo: a string like "26,Carbamidomethyl[C];" or a list like [(26, Carbamidomethyl[C])]
    '''

    line = re.split("\(|\)|\-", sequence)

    theoretical_peaks = {}
    bions, pepmass = calc_b_ions(line[0], modinfo, xlmassinfo)
    theoretical_peaks['b'] = bions

    theoretical_peaks['y'] = calc_y_from_b(bions, pepmass)
    #theoretical_peaks['b-ModLoss'] = calc_ion_modloss(bions, sequence, modinfo, N_term = True)
    #theoretical_peaks['y-ModLoss'] = calc_ion_modloss(theoretical_peaks['y'], sequence, modinfo, N_term = False)

    for i in range(int(line[1]), int(line[3])):
        theoretical_peaks['b'][i-1] = ""
        theoretical_peaks['y'][i-1] = ""

    theoretical_peaks['y'].reverse()
    # theoretical_peaks['y-ModLoss'].reverse()
    # , theoretical_peaks['b-ModLoss'], theoretical_peaks['y-ModLoss']
    return theoretical_peaks['b'], theoretical_peaks['y']


def get_theoretical_peaks_pL_xl(seq_mod_info, is_cleavable):
    '''
    Calc mz for ab/ay/bb/by ions from sequence and modinfo extracted from pLink
    '''

    if is_cleavable:
        seq_mod_info[['abS', 'ayS', 'bbS', 'byS', 'abL', 'ayL', 'bbL', 'byL', 'ab', 'ay', 'bb', 'by']] = np.array(list(map(
            lambda x, y, z: cal_theoretical_b_y_peaks_xl(x, y, z),
            seq_mod_info['Peptide'], seq_mod_info['Modifications'], seq_mod_info['Linker']
        )), dtype=object)
    else:
        seq_mod_info[['ab', 'ay', 'bb', 'by']] = np.array(list(map(
            lambda x, y, z: cal_theoretical_b_y_peaks_xl(x, y, z),
            seq_mod_info['Peptide'], seq_mod_info['Modifications'], seq_mod_info['Linker']
        )), dtype=object)

    return seq_mod_info


def get_theoretical_peaks_pL_regular(seq_mod_info):
    '''
    Calc mz for ab/ay/bb/by ions from sequence and modinfo extracted from pLink
    '''

    seq_mod_info[["b", "y"]] = np.array(list(map(
        lambda x, y: cal_theoretical_b_y_peaks(x, format_pL_modinfo(y)),
        seq_mod_info['Peptide'], seq_mod_info['Modifications']
    )), dtype=object)
    return seq_mod_info


def get_theoretical_peaks_pL_mono(seq_mod_info):
    '''
    Calc mz for ab/ay/bb/by ions from sequence and modinfo extracted from pLink
    '''

    seq_mod_info["linker_dict"] = np.array(list(map(
        lambda x, y: format_linker_mass(x, y, "mono"),
        seq_mod_info['Peptide'], seq_mod_info['Linker']
    )), dtype=object)

    seq_mod_info[["b", "y"]] = np.array(list(map(
        lambda x, y, z: cal_theoretical_b_y_peaks(
            x.split('(')[0], format_pL_modinfo(y), z),
        seq_mod_info['Peptide'], seq_mod_info['Modifications'], seq_mod_info["linker_dict"]
    )), dtype=object)
    return seq_mod_info


def get_theoretical_peaks_pL_loop(seq_mod_info):
    '''
    Calc mz for ab/ay/bb/by ions from sequence and modinfo extracted from pLink
    '''

    seq_mod_info["linker_dict"] = np.array(list(map(
        lambda x, y: format_linker_mass(x, y, "loop"),
        seq_mod_info['Peptide'], seq_mod_info['Linker']
    )), dtype=object)

    seq_mod_info[["b", "y"]] = np.array(list(map(
        lambda x, y, z: cal_theoretical_b_y_peaks_loop(
            x, format_pL_modinfo(y), z),
        seq_mod_info['Peptide'], seq_mod_info['Modifications'], seq_mod_info["linker_dict"]
    )), dtype=object)
    return seq_mod_info


def get_theoretical_peaks_pF(seq_mod_info):
    '''
    Calc mz for b/y ions from sequence and modinfo extracted from pFind
    '''

    seq_mod_info[["b", "y"]] = np.array(list(map(
        lambda x, y: cal_theoretical_b_y_peaks(x, format_pF_modinfo(y)),
        seq_mod_info['Sequence'], seq_mod_info['Modification']
    )), dtype=object)
    return seq_mod_info


def get_theo_peaks_array_from_precursor(sequence, modification, charge, consider_mod_loss=False):
    '''
    Calculate theoretical peaks from precursor.
    Used in search engine.
    '''

    # Generate 1-D array from seq_mod_info and sort them
    theo_peaks_array = []
    for ions_prefix, ions_mass_list in zip(
            ["b", "y"], cal_theoretical_b_y_peaks(sequence, modification)):
        for charge in range(1, charge+1):
            for length, mz in enumerate(ions_mass_list):
                if mz == "":
                    continue
                theo_peaks_array.append(((mz) / charge + aamass.mass_proton,
                                         f"{ions_prefix}{length+1}+{charge}"))
                if consider_mod_loss:
                    theo_peaks_array.append(((mz - aamass.mass_NH3) / charge + aamass.mass_proton,
                                             f"{ions_prefix}{length+1}-NH3+{charge}"))
                    theo_peaks_array.append(((mz - aamass.mass_H2O) / charge + aamass.mass_proton,
                                             f"{ions_prefix}{length+1}-H2O+{charge}"))

    theo_peaks_array.sort(key=lambda x: x[0])
    return theo_peaks_array


def get_theo_peaks_array(seq_mod_line, title, ions_prefix_list):
    '''
    Calculate theoretical peaks from theoretical ions.
    For example, given mass for b, calculate mz for b with different length and charge state.
    '''

    max_charge = int(re.findall(".*?\.", title)[3].split('.')[0])

    # Generate 1-D array from seq_mod_info and sort them
    theo_peaks_array = []
    for ions_prefix in ions_prefix_list:
        for charge in range(1, max_charge+1):
            for length, mz in enumerate(seq_mod_line[ions_prefix]):
                if mz == "":
                    continue
                theo_peaks_array.append(((mz) / charge + aamass.mass_proton,
                                         f"{ions_prefix}{length+1}+{charge}"))
                '''
                theo_peaks_array.append(((mz - aamass.mass_NH3) / charge + aamass.mass_proton, 
                    f"{ions_prefix}{length+1}-NH3+{charge}"))
                theo_peaks_array.append(((mz - aamass.mass_H2O) / charge + aamass.mass_proton,
                    f"{ions_prefix}{length+1}-H2O+{charge}"))
                '''

    theo_peaks_array.sort(key=lambda x: x[0])
    return theo_peaks_array


def get_theo_peaks_array_zero(seq_mod_line, title, ions_prefix_list):
    '''
    Calculate theoretical peaks from theoretical ions.
    For example, given mass for b, calculate mz for b with different length and charge state.
    '''

    # Generate 1-D array from seq_mod_info and sort them
    theo_peaks_array = []
    for ions_prefix in ions_prefix_list:
        for length, mz in enumerate(seq_mod_line[ions_prefix]):
            if mz == "":
                continue
            theo_peaks_array.append((mz,
                                     f"{ions_prefix}{length+1}"))
            '''
            theo_peaks_array.append(((mz - aamass.mass_NH3) / charge + aamass.mass_proton, 
                f"{ions_prefix}{length+1}-NH3+{charge}"))
            theo_peaks_array.append(((mz - aamass.mass_H2O) / charge + aamass.mass_proton,
                f"{ions_prefix}{length+1}-H2O+{charge}"))
            '''

    theo_peaks_array.sort(key=lambda x: x[0])
    return theo_peaks_array
