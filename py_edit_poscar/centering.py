#!/bin/python

#Using code JH_lib

#Purpose: Moving a poscar

import JH_lib as jh
import os
import numpy as np
from optparse import OptionParser
#1. Read.
import sys


if __name__ == '__main__':
    from time import time
#    filename = 'CONTCAR'
    t0 = time()
    filename = sys.argv[1]
    ## Lines for JH in SK hynix (180326)
    unitcell, compound, position = jh.r_cryst_vasp(filename)
    aver=0; i = 0
    for a in position: 
        aver += a[3]
        i += 1
#    print aver/float(i)
    new_position = jh.pst_poscar_move(position, 0.0, 0.0, aver/float(i))
                                    
#    if len(new_position[0]) == 7: 
    select =  True
#    else:
#      select =  False
      
    jh.w_poscar(position, compound = compound, filename = 'new.vasp', unitcell = unitcell, Selective = select)

    t1 = time()


