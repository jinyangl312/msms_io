
'''
Created on 2013-8-16

@author: RunData
'''
from .modification import get_modification, keep_one_neutral_loss


class AAMass(object):
    '''
    Get mass of amino acids.
    Forked from pDeep.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.mass_H = 1.0078250321
        self.mass_O = 15.9949146221
        self.mass_proton = 1.007276
        self.mass_N = 14.0030740052
        self.mass_C = 12.00
        self.mass_isotope = 1.003

        self.mass_H2O = self.mass_H * 2 + self.mass_O
        self.mass_CO = self.mass_C + self.mass_O
        self.mass_CO2 = self.mass_C + self.mass_O * 2
        self.mass_NH = self.mass_N + self.mass_H
        self.mass_NH3 = self.mass_N + self.mass_H * 3
        self.mass_HO = self.mass_H + self.mass_O

        self.aa_mass_dict = {}
        self.aa_mass_dict['A'] = 71.037114
        # self.aa_mass_dict['B'] = 0. # B  aspartate/asparagine? from IUPAC
        self.aa_mass_dict['C'] = 103.009185
        self.aa_mass_dict['D'] = 115.026943
        self.aa_mass_dict['E'] = 129.042593
        self.aa_mass_dict['F'] = 147.068414
        self.aa_mass_dict['G'] = 57.021464
        self.aa_mass_dict['H'] = 137.058912
        self.aa_mass_dict['I'] = 113.084064
        # self.aa_mass_dict['J'] = 114.042927 # ? or None?
        self.aa_mass_dict['K'] = 128.094963
        self.aa_mass_dict['L'] = 113.084064
        self.aa_mass_dict['M'] = 131.040485
        self.aa_mass_dict['N'] = 114.042927
        #self.aa_mass_dict['O'] = 0. #
        self.aa_mass_dict['P'] = 97.052764
        self.aa_mass_dict['Q'] = 128.058578
        self.aa_mass_dict['R'] = 156.101111
        self.aa_mass_dict['S'] = 87.032028
        self.aa_mass_dict['T'] = 101.047679
        # self.aa_mass_dict['U'] = 150.95363 # ?
        self.aa_mass_dict['V'] = 99.068414
        self.aa_mass_dict['W'] = 186.079313
        self.aa_mass_dict['X'] = 0.  # any?
        self.aa_mass_dict['Y'] = 163.06332
        # self.aa_mass_dict['Z'] = 0. # glutamate/glutamine?

        self.glyco_mass_dict = {}
        self.glyco_mass_dict["Xyl"] = 132.0422587452
        self.glyco_mass_dict["Hex"] = 162.0528234315
        self.glyco_mass_dict["dHex"] = 146.0579088094
        self.glyco_mass_dict["HexNAc"] = 203.07937253300003
        self.glyco_mass_dict["NeuAc"] = 291.09541652769997
        self.glyco_mass_dict["NeuGc"] = 307.09033114979997

        self.mod_mass_dict = keep_one_neutral_loss(get_modification())

    def fix_C57(self):
        self.aa_mass_dict['C'] += 57.021464


class AAMass_aa_ini(AAMass):
    '''
    Get mass of amino acids from aa.ini.
    Use mono mass.
    Used in search engine.
    '''

    def __init__(self):
        super().__init__()
        self.get_from_aa_ini()

    def get_from_aa_ini(self):
        from .aa import get_aa
        for aa_pFind, aamass_pFind in get_aa("mono").items():
            self.aa_mass_dict[aa_pFind] = aamass_pFind


aamass = AAMass()
aamass_aa_ini = AAMass_aa_ini()
