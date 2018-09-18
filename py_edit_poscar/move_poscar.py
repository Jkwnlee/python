#!/bin/python

#Using code JH_lib

#Purpose: Moving a poscar

import JH_lib as jh
import os
#1. Read.
import sys

print """[CODE] run this commend with target POSCAR file name and distance which you want:
[CODE] for example: python /team/ptcad/jhlee/b_codework/py_move-poscar/move_poscar.py POSCAR 1 1 1
[CODE]              means: moving atoms in given POSCAR 1 Angstrom in x-axis and 1 A in y-axis and 1 A z-axis
[CODE]              for the direction, please check after test this code.
"""

## Lines for JH in SK hynix (180326)
filename = sys.argv[1]
unitcell, compound, position = jh.r_cryst_vasp(filename)
new_position = jh.pst_poscar_move(position, float(sys.argv[2]), float(sys.argv[3]),float(sys.argv[4]))
jh.w_poscar(position, compound = compound, filename = 'new.vasp', unitcell = unitcell, Selective = False)

