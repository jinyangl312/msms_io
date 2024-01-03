# _*_ coding: utf-8 _*_
# @Time    :   2022/05/10 17:49:34
# @FileName:   element.py
# @Author  :   jyl
import re
import os


def get_element():
    """
    Interface for the usage of element.ini
    """

    element_dict = {}  # {element_name: elements}
    with open(os.path.dirname(__file__)+'/element.ini') as f:
        line = f.readline()
        number = int(re.search('\d+', line).group())
        for i in range(number):
            line = re.split("=|\n", f.readline())[1]
            line = re.split("\|", line)
            element_dict[line[0]] = ([float(s) for s in re.findall(
                "\d+\.?\d*", line[1])], [float(s) for s in re.findall("\d+\.?\d*", line[2])])
    return element_dict
