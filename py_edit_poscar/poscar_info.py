#!/bin/python
#Using code JH_lib
#Purpose: choose a poscar
import JH_lib as jh
import os
import math
import numpy as np
from optparse import OptionParser

#1. Read.
import sys
############################################################
__version__ = "1.0"
#Print the information as table format
# Volume
# Species : Number of atoms 
############################################################

def command_line_arg():
    usage = """[CODE] usage: %prog [options] arg**
[CODE] run this commend with target POSCAR file name
[CODE] for ex1: python /team/ptcad/jhlee/b_codework/py_edit_poscar/expand_poscar.py -p POSCAR
[CODE]          To print out the volume, number of atoms for each kind
"""
#[CODE] output: given file name + _target.vasp
    par = OptionParser(usage=usage, version= __version__)
     
    par.add_option("-s", '--Selective', 
                  action='store_false', dest='selec',
                  default=True,
                  help='Turn on the Selective Dynamics to T T T(Default: True)')

    par.add_option("-p", '--poscar', 
                  action='store',type='string', dest='poscar',
                  default='./POSCAR',
                  help='Choose POSCAR : D-./POSCAR)')

    par.add_option("-c", '--condition', nargs=2,
                  action='store',type='float', dest='condi',
                  default=False,
                  help='Choose elements within the range: D-False)')

    par.add_option("-e", '--element', 
                  action='append', type='string', dest='t_ele',
                  default=[],
                  help='Choose elements within the POSCAR : [])')

    par.add_option("-o", '--output',
                  action='store',type='string', dest='out',
                  default=False,
                  help='Output file name: filename)')


    return   par.parse_args( )

############################################################

def main(opts):
    #print "[CODE] Start to read the given structure: ", opts.poscar

    unitcell, compound, position = jh.r_cryst_vasp(opts.poscar)
    a = np.array(unitcell[0])
    b = np.array(unitcell[1])
    c = np.array(unitcell[2])
    volume = np.dot(np.cross(a,b),c)
    vec_a= np.sum(a**2) ** (0.5)
    vec_b= np.sum(b**2) ** (0.5)
    vec_c= np.sum(c**2) ** (0.5)
    alpha= np.degrees(np.arccos(np.dot(b,c)/vec_b/vec_c))
    beta = np.degrees(np.arccos(np.dot(a,c)/vec_a/vec_c))
    gamma= np.degrees(np.arccos(np.dot(a,b)/vec_a/vec_b))
    print('a %8.8s : %20.10f' %('alpha' ,alpha))
    print('a %8.8s : %20.10f' %('beta'  ,beta))
    print('a %8.8s : %20.10f' %('gamma' ,gamma))
    print('l %8.8s : %20.10f' %('vec_a' ,vec_a))
    print('l %8.8s : %20.10f' %('vec_b' ,vec_b))
    print('l %8.8s : %20.10f' %('vec_c' ,vec_c))
    print('V %8.8s : %20.10f' %('volume',volume))
    for a,b in zip(compound[0],compound[1]):
        print('C %8.8s : %20.i' %(a,int(b)))
                    

############################################################
############################################################

if __name__ == '__main__':
    from time import time
    opts, args = command_line_arg()
    t0 = time()
    main(opts)
    t1 = time()
    
