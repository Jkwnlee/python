#!/bin/python

#Using code JH_lib

#Purpose: Moving a poscar

import JH_lib as jh
import os
import numpy as np
from optparse import OptionParser
#1. Read.
import sys
############################################################
__version__ = "1.0"
############################################################

print """[CODE] run this commend with target POSCAR file name and distance which you want:
[CODE] for example: python /team/ptcad/jhlee/b_codework/py_move-poscar/move_poscar.py POSCAR 1 1 1
[CODE]              means: moving atoms in given POSCAR 1 Angstrom in x-axis and 1 A in y-axis and 1 A z-axis
[CODE]              for the direction, please check after test this code.
"""


if __name__ == '__main__':
    from time import time
    t0 = time()
    filename = sys.argv[1]

    unitcell, compound, position = jh.r_cryst_vasp(filename)
    new_position = jh.pst_poscar_move(position, float(sys.argv[2]), \
                                    float(sys.argv[3]),float(sys.argv[4]))
                                    
    if len(new_position[0]) == 7: 
      select =  True
    else:
      select =  False
      
    jh.w_poscar(position, compound = compound, filename = 'new.vasp', unitcell = unitcell, Selective = select)

    t1 = time()
