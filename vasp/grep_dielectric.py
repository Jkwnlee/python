#!/bin/python

### The python code for read  OUTCAR and print certain properties
### Written by Ji-Hwan Lee
### Last update: 2019. 4. 17
### Used library:os, math, matplotlib in py_vasp_post_locpot_to_1D_fun_lib
### Detail functions are in 'py_vasp_post_locpot_to_1D_fun_lib.py'
### The python code for 
### Written by Ji-Hwan Lee


import os
import math 
import sys
import numpy as np
import JH_lib as jh
from optparse import OptionParser

############################################################
__version__ = "1.0"
# 1.4
# Update for target information
# Update for optparse (--input, --draw, --output)
############################################################


def read_file_line(filename, key):
    i = 0
    k = 0 
    with open(filename, 'r') as f:
        for line in f:
            if key in line: 
                k = i
            i += 1
    return k

############################################################
############################################################

def average_xx_yy_zz(tensor):
   return (float(tensor[0][0]) +  float(tensor[1][1]) + float(tensor[2][2]))/3

############################################################
############################################################

def get_static_dielectric(filename):
    if os.path.isfile(filename) == True:
        strings  = "MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION"
        line_num = read_file_line(filename, strings)
        if line_num ==0 : return np.zeros((3, 3)),0
        else:
            static_e=[]
            with open(filename, 'r') as f:
                i = 0 
                for line in f:
                    if i >  line_num + 1 and i < line_num + 5 : 
                        static_e.append(line.split())
    #                    print line
                    i = i + 1
            return static_e, average_xx_yy_zz(static_e)
    else:
        print("[CODE] Error, There is no %s" %filename)
        return np.zeros((3, 3)),0


####################################################################
####################################################################

   
def command_line_arg():
    usage = "usage: %prog [options] arg1"  
    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-o", '--outcar', 
            action='store', type="string", dest='outcar',
            default='./OUTCAR',
            help='Path for the OUTCAR, default: .')

    par.add_option("-t", '--tensor', 
            action='store_true', dest='tensor',
            default=False,
            help='Print the Whole tensor: \n Default: print only average 3 component xx,yy,zz')

    return  par.parse_args( )
    
####################################################################
####################################################################



if __name__ == "__main__":
    opts, args = command_line_arg()
    t, e = get_static_dielectric(opts.outcar)
    if opts.tensor and e != 0: 
        for a in t: 
            for b in a: 
                print "%10.5f" %float(b),
            print
    else: print e

