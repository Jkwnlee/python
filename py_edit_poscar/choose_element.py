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
__version__ = "1.2"
#ADD direct coordinate option (choose atoms between given coordinate)
############################################################

def command_line_arg():
    usage = """[CODE] usage: %prog [options] arg1 arg2
[CODE] run this commend with target POSCAR file name and elements which you want:
[CODE] for ex1: python /team/ptcad/jhlee/b_codework/py_edit_poscar/expand_poscar.py -p POSCAR -e Ti N
[CODE]          TO find Ti and N, and build new POSCAR with selected element in the given unitcell
[CODE] for ex2: python /team/ptcad/jhlee/b_codework/py_edit_poscar/expand_poscar.py -p POSCAR -c 0.45 0.73
[CODE]          TO find element between the relative position of 0.45 and 0.73 in z-direction, 
[CODE]          and build new POSCAR with selected element in the given unitcell
[CODE] output: given file name + _target.vasp"""
 
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
    print "[CODE] Start to read the given structure: ", opts.poscar
    if opts.out: output_name=opts.out
    else:    output_name = opts.poscar+'_targets.vasp'

    new_position =[];     temp_position=[]; select = True
    unitcell, compound, position = jh.r_cryst_vasp(opts.poscar)
    if opts.condi:
        print '[CODE] the -c option is chosen (not the elements)'
        if opts.condi[1] < opts.condi[0]: 
            reverse=True
            print '[CODE] [0] %f is higher than [1] %f' %(opts.condi[0], opts.condi[1])
            print '[CODE] atoms higher than [0] %f and lower then [1] %f are chosen' %(opts.condi[0], opts.condi[1])
        else : 
            reverse=False
            print '[CODE] [1] %f is higher than [0] %f' %(opts.condi[1], opts.condi[0])
            print '[CODE] choose element higher than %f and lower than %f' %(opts.condi[0],opts.condi[1])

        for a in range(len(position)):
            if opts.condi[0] > 1 or opts.condi[1] >1:            z= 1
            else:                z=unitcell[2][2] 
            if reverse:

                if position[a][3]/z >  opts.condi[0] or \
                   position[a][3]/z <  opts.condi[1]:
             #       print a, len(new_position), position[a]
                    temp_position = position[a]
                    new_position.append(temp_position)
            else:
                if position[a][3]/z >  opts.condi[0] and \
                   position[a][3]/z <  opts.condi[1]:
             #       print position[a][3]/unitcell[2][2], opts.condi
                    temp_position = position[a]
                    new_position.append(temp_position)

    else:
        print "[CODE] Among the elements in the target file %s containing," %opts.poscar, compound[0], ", you choose:", opts.t_ele
        for a in range(len(position)):
            for b in opts.t_ele:
                if position[a][0] == b:
                    temp_position = position[a]
                    temp_position[0] = b
                    new_position.append(temp_position)
    if len(position[0]) == 7 and opts.selec == False: 
        select =    False
    else:
        select =    opts.selec
    new_compound, new_position1 = jh.component_from_positions(new_position)

    if len(new_position) == 0 :
        print "error!!!"
    else:
        jh.w_poscar(position = new_position1,\
               compound  = new_compound, \
               filename  = output_name,    \
               unitcell  = unitcell,         \
               Selective = select)
                    

############################################################
############################################################

if __name__ == '__main__':
    from time import time
    opts, args = command_line_arg()
    t0 = time()
    main(opts)
    t1 = time()
    
