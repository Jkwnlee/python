#!/bin/python

#Using code JH_lib

#Purpose: Moving a poscar

import JH_lib as jh
import os
#1. Read.
import sys


filename = sys.argv[1]
epsilon = [float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) ]

print """
[CODE] run this commend with target POSCAR file name and strain which you want:
[CODE] for example: python /team/ptcad/jhlee/b_codework/py_edit_poscar/expand_poscar.py POSCAR 0.1 0.1 0
[CODE]              means: expand cell with given strain 0.1 and 0.1 in x-axis and y-axis and keeping z-axis
[CODE]              for some unit cell not keeping cubic/tetragonal like system (non-vertical), check the result
[CODE]              output: %s
""" %(filename+'_new.vasp')

## Lines for JH in SK hynix (180326)

unitcell, compound, position = jh.r_cryst_vasp(filename)
#print epsilon, unitcell
for a in range(3):
 for b in range(3):
#    print a,b
   unitcell[a][b] = float(unitcell[a][b]) * (1 + epsilon[a])

for a in range(len(position)):
  position[a][1] = float(position[a][1]) * (1 + epsilon[0])
  position[a][2] = float(position[a][2]) * (1 + epsilon[1])
  position[a][3] = float(position[a][3]) * (1 + epsilon[2])
  
  
jh.w_poscar(position = position, \
            compound = compound, \
            filename = filename+'_new.vasp', \
            unitcell = unitcell, \
            Selective = False)
