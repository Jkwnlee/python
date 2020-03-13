#!/bin/python

#Using code JH_lib

#Purpose: Moving a poscar


import JH_lib as jh
import os
import sys

############################################################
__version__ = 1.1
############################################################

def expand_cell(filename,expantion,output_name, Selec):
  print "[CODE]              expand the given structure  ", filename
  print """[CODE]              output: %s with selective dynamics %s
""" %(output_name,out)
  unitcell, compound, position = jh.r_cryst_vasp(filename)
  unitcell, compound, position = jh.pst_cell_expansion(unitcell, position, expantion )
  jh.w_poscar(position = position, \
              compound = compound, \
              filename = output_name,\
              unitcell = unitcell, \
              Selective= Selec)


############################################################
############################################################

if __name__ == '__main__':
  from time import time
#    opts, args = command_line_arg()
  t0 = time()
  
  filename = sys.argv[1];     out = "off";     Selec = False
  expantion =[int(sys.argv[2]), int(sys.argv[3]),int(sys.argv[4])]

  if len(sys.argv) >= 5:
    output_name=filename+'_new.vasp'
  if  len(sys.argv) >= 6:
    output_name = sys.argv[5]
  if  len(sys.argv) >= 7:
    Selec = True
    out = "On"


  if  len(sys.argv) < 5:
    print """
[CODE] run this commend with target POSCAR file name and expanded cell which you want:
[CODE] for example: python /team/ptcad/jhlee/b_codework/py_edit_poscar/expand_poscar.py POSCAR 2 2 1
[CODE]              means: expand cell in given POSCAR 2 times in x-axis and y-axis and keeping z-axis
[CODE]              for the direction, please check after test this code.
[CODE]              output: %s with selective dynamics %s
""" %(output_name,out)

  else:
    expand_cell(filename,expantion,output_name,Selec)

  t1 = time()
#  print '\nCode calc completed! Time Used: %.2f [sec]\n' % (t1 - t0)
 

