#!/bin/python

#Using code JH_lib

#Purpose: choose a poscar

import JH_lib as jh
import os
import numpy as np
from optparse import OptionParser
#1. Read.
import sys
############################################################
__version__ = "1.1"
############################################################

def command_line_arg():
    usage = "usage: %prog [options] arg1 arg2"  
    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-s", '--Selective', 
            action='store_true', dest='Selec',
            default=False,
            help='Turn on the Selective Dynamics to T T T(Default: False)')


    return  par.parse_args( )

############################################################

def choose_cell(inputs):
  print "doping the given structure  ", inputs[1]
  filename    = inputs[1]
  target_ele  = inputs[2:]
#  num_ta_ele  = int(inputs[3])
#  change_ele  = inputs[4]
  output_name = filename+'_targets.vasp'
  unitcell, compound, position = jh.r_cryst_vasp(filename)
  print "Among the elements in the target file %s containing," %filename, compound[0], ", you choose:", target_ele
  #print epsilon, unitcell
  new_position = []
  temp_position=[]
  for a in range(len(position)):
    for b in target_ele:
      if position[a][0] == b:
        temp_position = position[a]
        temp_position[0] = b
        new_position.append(temp_position)
#        print temp_position[0]
#  print temp_position,num_ta_ele,target_ele
#  for a in new_position:
#  for a in new_position: print a
  if len(new_position[0]) == 7 and opts.Selec == False: 
      select =  False
  else:
      select =  opts.Selec
      
  new_compound, new_position1 = jh.component_from_positions(new_position)
  jh.w_poscar(position = new_position1,\
              compound = new_compound, \
              filename = output_name,  \
              unitcell = unitcell,     \
              Selective = select)
          

############################################################
############################################################

if __name__ == '__main__':
  from time import time
  opts, args = command_line_arg()
  t0 = time()
  
  filename='filename'
  if len(sys.argv) > 2: 
    choose_cell(sys.argv)
  else:
    print """
[CODE] run this commend with target POSCAR file name and elements which you want:
[CODE] for example: python /team/ptcad/jhlee/b_codework/py_edit_poscar/expand_poscar.py POSCAR Ti N
[CODE]              means: find Ti and N, and build new POSCAR with selected element in the given unitcell
[CODE]              output: %s
""" %(filename+'_targets.vasp')



  t1 = time()
#  print '\nCode calc completed! Time Used: %.2f [sec]\n' % (t1 - t0)
  
  
