# -*- coding: utf-8 -*-
"""
仅window下可用,需要调用安装Chemdraw时安装的python
example:  C:\Python32\python.exe .\name2smi.py "heptyl"
"""

import sys
from ChemScript16 import StructureData

if __name__ == '__main__':
    name = sys.argv[1]
    name = '('+ name +')chloride-37Cl'  #对取代基后添加[37Cl]定位
    m = StructureData()
    name = name.replace("_"," ")
    print(name)
    m.ReadData(name)
    smiles = m.Smiles
    print(smiles)