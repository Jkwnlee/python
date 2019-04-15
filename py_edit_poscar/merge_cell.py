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

def command_line_arg():
    usage = "usage: %prog [options] arg1 arg2"  
    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-s", '--Selective', 
            action='store_true', dest='Selec',
            default=False,
            help='Turn on the Selective Dynamics to T T T(Default: False)')


    return  par.parse_args( )

############################################################

if __name__ == '__main__':
  from time import time
  opts, args = command_line_arg()
  t0 = time()
  filename1 = sys.argv[1]
  filename2 = sys.argv[2]
  output_name = filename1+'_'+filename2+'_merged.vasp'
  ## Lines for JH in SK hynix (180326)
  unitcell1, compound1, position1 = jh.r_cryst_vasp(filename1)
  unitcell2, compound2, position2 = jh.r_cryst_vasp(filename2)
  if unitcell1 == unitcell2: 
    #pass
    print "same!"
  else: 
    print sys.stderr
       
  new_position = position1 + position2
  if len(new_position[0]) == 7 and opts.Selec == False: 
      select =  False
  else:
      select =  opts.Selec
      
  new_compound, new_position1 = jh.component_from_positions(new_position)
  jh.w_poscar(position = new_position1,\
              compound = new_compound, \
              filename = output_name,  \
              unitcell = unitcell1,     \
              Selective = select)

  t1 = time()


