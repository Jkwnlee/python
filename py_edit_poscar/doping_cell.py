#!/bin/python

#Using code JH_lib

#Purpose: Doping a poscar


import JH_lib as jh
import os
import sys
import re
import numpy as np
from optparse import OptionParser

############################################################
__version__=1.1
############################################################

def doping_cell(opts):
    print "[CODE] doping the given structure  ", opts.poscar
    outname = opts.poscar+'_new.vasp'
    unitcell, compound, position = jh.r_cryst_vasp(opts.poscar)
    new_position = []
    if len(opts.select) > 0:
        i = 1
        print '[CODE] You choose a mode to doping only -An-Atom %s%i to %s'\
            %(opts.select[0][0], int(opts.select[0][1]), opts.select[0][2])
        for a in range(len(position)):
            temp_position=[]
            if position[a][0] == opts.select[0][0]:
                if i == int(opts.select[0][1]): 
                    temp_position = position[a]
                    temp_position[0] = opts.select[0][2]
                    new_position.append(temp_position)
                else: 
                    new_position.append(position[a])
                i = i + 1
            else: new_position.append(position[a])
    elif opts.target_ele and opts.dopant_ele:
        wmin, wmax = opts.window
        print '[CODE] You choose a mode to doping Atoms (%s) in specific region from %f to %f along z-axis to %s'\
        %(opts.target_ele, wmin, wmax, opts.dopant_ele)
        for a in range(len(position)):
            temp_position=[]
            if position[a][0] == opts.target_ele and position[a][3]<wmax and position[a][3] > wmin:
                temp_position = position[a]
                temp_position[0] = opts.dopant_ele
                new_position.append(temp_position)
            else: 
                new_position.append(position[a])
    else: 
        print '[CODE] YOU NEED TO CHECK help (-h) option'


    new_compound, new_position1 = jh.component_from_positions(new_position)
    jh.w_poscar(position = new_position1, \
              compound = new_compound, \
              filename = outname, \
              unitcell = unitcell, \
              Selective = False)



############################################################
############################################################
def command_line_arg():
    usage = """
[CODE] run this commend with target POSCAR file name and strain which you want:
[CODE] for example: python /team/ptcad/jhlee/b_codework/py_edit_poscar/expand_poscar.py -i POSCAR -s Hf 3 Si
[CODE]              means: Change 3rd Hf atom to Si (the number of element is checked by VESTA)
    usage: %prog [options] arg1 arg2"""

    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-i", '--input', 
            action='store', type="string", dest='poscar',
            default='POSCAR',
            help='location of the POSCAR')

    par.add_option("-a", '--target', 
            action='store', type="string", dest='target_ele',
            default=False,
            help='specify which atoms to substitute')

    par.add_option("-d", '--dopant', 
            action='store', type="string", dest='dopant_ele',
            default=False,
            help='specify which atoms to dope')

    par.add_option("-w", '--window', nargs=2,
            action='store', type="float", dest='window',
            default=(0,0),
            help='Minimul and maximum z-value where the target element in')

    par.add_option('-s','--select',nargs=3, 
            action='append',  dest='select',
            default=[],
            help='''choose specific atom to others: 
--select (-s) (target-element) (target-element number) (element-to-dope)
            ''')

    return  par.parse_args( )

############################################################

if __name__ == '__main__':
    from time import time
    opts, args = command_line_arg()

    t0 = time()
    doping_cell(opts)

    t1 = time()
    print '\nCode calc completed! Time Used: %.2f [sec]\n' % (t1 - t0)
  
  
