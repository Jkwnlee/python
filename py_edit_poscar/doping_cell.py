#!/bin/python

#Using code JH_lib

#Purpose: Moving a poscar


############################################################
############################################################

import JH_lib as jh
import os
import sys

def doping_cell(inputs):
  print "doping the given structure  ", inputs[1]
  filename    = inputs[1]
  target_ele  = inputs[2]
  num_ta_ele  = int(inputs[3])
  change_ele  = inputs[4]
  output_name = filename+'_new.vasp'
  unitcell, compound, position = jh.r_cryst_vasp(filename)
  #print epsilon, unitcell
  i = 1
  new_position = []
  temp_position=[]
  for a in range(len(position)):
    if position[a][0] == target_ele and i == num_ta_ele: 
      temp_position = position[a]
      temp_position[0] = change_ele
#      print temp_position[0], i
    else:
      new_position.append(position[a])
    i = i + 1
  new_position.append(temp_position)
#  for a in new_position:
#@    print a
  new_compound, new_position1 = jh.component_from_positions(new_position)
#  print position1, new_compound
  jh.w_poscar(position = new_position1, \
              compound = new_compound, \
              filename = filename+'_new.vasp', \
              unitcell = unitcell, \
              Selective = False)
          

############################################################
############################################################

if __name__ == '__main__':
  from time import time
#    opts, args = command_line_arg()
  t0 = time()
  
  filename='filename'
  if len(sys.argv) == 5: 
    doping_cell(sys.argv)
  else:
    print """
[CODE] run this commend with target POSCAR file name and strain which you want:
[CODE] for example: python /team/ptcad/jhlee/b_codework/py_edit_poscar/expand_poscar.py POSCAR Hf 3 Si
[CODE]              means: Change 3rd Hf atom to Si (the number of element is checked by VESTA)
[CODE]              output: %s
""" %(filename+'_new.vasp')



  t1 = time()
#  print '\nCode calc completed! Time Used: %.2f [sec]\n' % (t1 - t0)
  
  
