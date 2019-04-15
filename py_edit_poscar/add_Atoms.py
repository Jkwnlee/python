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

####################################################################
####################################################################

   
def command_line_arg():
    usage = "usage: %prog [options] arg1"  
    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-f", '--filename1', 
            action='store', type="string", 
            dest='filename1', default='out.vasp',
            help='output name of the code')

    par.add_option("-F", '--filename2', 
            action='store', type="string", 
            dest='filename2', default='out.vasp',
            help='output name of the code')

    par.add_option("-o", '--output', 
            action='store', type="string", 
            dest='outname', default='out.vasp',
            help='output name of the code')

    par.add_option("-t", '--two', 
            action='store', type="string", dest='twoside',
            default=False,
            help='add element both side or only one-side: Default-Oneside')

    par.add_option("-d", '--distance', 
            action='store', type="float",
            dest='dist', default=2.2,
            help='distance from the maximum/minimun of atom in the given slab, default 2.2')

    return  par.parse_args( )
    
####################################################################
####################################################################



if __name__ == '__main__':
    from time import time
    opts, args = command_line_arg()

    t0 = time()
#    filename1 = sys.argv[1]
#    filename2 = sys.argv[2]
    ## Lines for JH in SK hynix (180326)
    sub_unitcell, sub_compound, sub_position = jh.r_cryst_vasp(opts.filename1)
    add_unitcell, add_compound, add_position = jh.r_cryst_vasp(opts.filename2)

    z=[]; 
    for a in sub_position: 
        z.append(a[3])
    
    if opts.twoside:
        min_add=[]  ;  max_add=[]
    else:  max_add=[]

    for a in add_position:
         if opts.twoside:
           mi =  - float(a[3]) + ( min(z) - opts.dist )
           ma =    float(a[3]) + ( max(z) + opts.dist )
           min_add.append([a[0],a[1],a[2], mi])
           max_add.append([a[0],a[1],a[2], ma])
         else:
           max_add.append([a[0],a[1],a[2],   float(a[3]) + ( max(z) + opts.dist )])

    if opts.twoside:
      jh.add_adsorbate(sub_position, sub_unitcell, min_add, 't1.vasp')
      sub_unitcell, sub_compound, sub_position = jh.r_cryst_vasp('t1.vasp')
      jh.add_adsorbate(sub_position, sub_unitcell, max_add, opts.outname)
      os.remove('t1.vasp')
    else: 
      jh.add_adsorbate(sub_position, sub_unitcell, max_add, opts.outname)

    t1 = time()





