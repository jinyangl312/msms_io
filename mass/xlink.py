# _*_ coding: utf-8 _*_
# @Time    :   2022/03/01 20:00:19
# @FileName:   xlink.py
# @Author  :   jyl
import re
import os


def get_linker_dict():
    """
    Interface for the usage of xlink.ini
    """

    linker_dict = {}  # {mod_name: elements}
    with open(os.path.dirname(__file__)+'/xlink.ini') as f:
        # Header
        line = f.readline()  # [xlink]
        line = f.readline()
        number = int(re.search('\d+', line).group())

        for _ in range(number):
            line = f.readline()
            line = re.split("=|\n", f.readline())
            linker_dict[line[0]] = line[1].split(' ')
    return linker_dict


class LinkerMass(object):
    """
    Access linker mass info
    """

    def __init__(self):
        self.linker_dict = get_linker_dict()

    def get_mono_mass(self, linker):
        return float(self.linker_dict[linker][4])

    def get_linker_mass(self, linker):
        return float(self.linker_dict[linker][2])

    def is_cleavable(self, linker):
        if linker not in self.linker_dict:
            return False
        return self.linker_dict[linker][8] != '0'

    def get_short_arm_mass(self, linker):
        if self.is_cleavable(linker):
            return float(self.linker_dict[linker][10])
        return None

    def get_long_arm_mass(self, linker):
        if self.is_cleavable(linker):
            return float(self.linker_dict[linker][9])
        return None

    def get_long_short_arm_mass_deviation(self, linker):
        if self.is_cleavable(linker):
            return float(self.linker_dict[linker][9]) - float(self.linker_dict[linker][10])
        return None


xlmass = LinkerMass()
