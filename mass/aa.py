# _*_ coding: utf-8 _*_
# @Time    :   2022/05/10 17:50:07
# @FileName:   aa.py
# @Author  :   jyl
import re
import os
from .element import get_element


def get_aa(mode="mono"):
    """
    Interface for the usage of aa.ini
    """

    element_dict = get_element()
    aa_dict = {}  # {aa_name: mono_mass}
    with open(os.path.dirname(__file__)+'/aa.ini') as f:
        # header
        line = f.readline()
        number = int(re.search('\d+', line).group())
        line = f.readline()

        if mode == "mono":
            for _ in range(number):
                line = re.split("=|\n", f.readline())[1]
                line = re.split("\|", line)  # C(3)H(5)N(1)O(1)S(0)
                composition = re.split(
                    "\(|\)", line[1].strip(")"))  # C3H5N1O1S0

                # Calculate mass of each amino acid by monoisotope element mass
                aa_mass = 0.0
                for j in range(int(len(composition)/2)):
                    element_mono_mass = element_dict[composition[j*2]][0][0]
                    aa_mass += element_mono_mass * int(composition[j*2+1])
                aa_dict[line[0]] = aa_mass

        elif mode == "avg":
            for _ in range(number):
                line = re.split("=|\n", f.readline())[1]
                line = re.split("\|", line)
                composition = re.split("\(|\)", line[1].strip(")"))

                # Calculate mass of each amino acid by average element mass
                aa_mass = 0.0
                for j in range(int(len(composition)/2)):
                    for element_mass, element_abundance in zip(*element_dict[composition[j*2]]):
                        element_avg_mass = element_mass * element_abundance
                        aa_mass += element_avg_mass * int(composition[j*2+1])
                aa_dict[line[0]] = aa_mass
        return aa_dict
