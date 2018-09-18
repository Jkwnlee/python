#!/bin/python

#Using code JH_lib

#Purpose: Moving a poscar

import JH_lib as jh
import os
#1. Read.
import sys

## Lines for JH in SK hynix (180326)
filename = sys.argv[1]
unitcell, compound, position = jh.r_cryst_vasp(filename)

print """
[CODE] run this commend with target POSCAR file name and expanded cell which you want:
[CODE] for example: python /team/ptcad/jhlee/b_codework/py_edit_poscar/change_unitcell.py POSCAR 0 0 15
[CODE]              means: change cell with giving 0, 0, and 15 angstrom in x-axis and y-axis and keeping z-axis
[CODE]              for the direction, please check after test this code.
[CODE]              output: %s
""" %(filename+'_new.vasp')

new_unitcell = []
new_unitcell.append([float(unitcell[0][0]) + float(sys.argv[2]), unitcell[0][1], unitcell[0][2]])
new_unitcell.append([unitcell[1][0], float(unitcell[1][1]) + float(sys.argv[3]), unitcell[1][2]])
new_unitcell.append([unitcell[2][0], unitcell[2][1], float(unitcell[2][2]) + float(sys.argv[4])])

jh.w_poscar(position, compound = compound, filename = filename+'_new.vasp', unitcell = new_unitcell, Selective = False)
