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
import pandas as pd
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

def get_born_effective(opts):
    in_unitcell, in_compound, in_position = jh.r_cryst_vasp(opts.poscar)
    n_atom = len(in_position)
    d = pd.DataFrame(np.zeros((10, n_atom)), 
          columns=['atom','xx','xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'])
    j = 0
    for a, b in zip(in_compound[0], in_compound[1]):
        for i in range(int(b)):
            d.atom.iloc[j] = a+'%s' %(int(i)+1)
            j = j + 1

    if os.path.isfile(opts.outcar) == True:
        strings  = "BORN EFFECTIVE CHARGES (in e, cummulative output)"
        start_num = read_file_line(opts.outcar, strings)
        if start_num ==0 : return d
        else:
            with open(opts.outcar, 'r') as f:
                i = 0 
                for line in f:
                    if i >  start_num + 1 and i < start_num + n_atom * 4 + 2: 
                        if 'ion'  in line: 
                            index_num = int(line.split()[1])-1
                        if (i - start_num -2) % 4  == 1:
                            for col_num in range(3):
                                d.iloc[index_num, col_num +1] = float(line.split()[col_num + 1]) 
                        elif (i - start_num -3) % 4  == 1:
                            for col_num in range(3):
                                d.iloc[index_num, col_num +4] = float(line.split()[col_num + 1]) 
                        elif (i - start_num -4) % 4  == 1:
                            for col_num in range(3):
                                d.iloc[index_num, col_num +7] = float(line.split()[col_num + 1]) 
                    i = i + 1
    else:
        print("[CODE] Error, There is no %s/OUTCAR" %filename)
    return d


####################################################################
####################################################################

   
def command_line_arg():
    usage = "usage: %prog [options] arg1"  
    par = OptionParser(usage=usage, version= __version__)

    par.add_option("--outcar", '--input', 
            action='store', type="string", dest='outcar',
            default='./OUTCAR',
            help='Path for the OUTCAR, default: .')

    par.add_option("--poscar", '--structure', 
            action='store', type="string", dest='poscar',
            default='./POSCAR',
            help='Path for the POSCAR, default: .')

    par.add_option("--out", '--output', 
            action='store_true', dest='output',
            default=False,
            help='Print the BORN EFFECTIVE CHARGE tensor: \n Default: print dataframe only')

    return  par.parse_args( )
    
####################################################################
####################################################################



if __name__ == "__main__":
    opts, args = command_line_arg()
    born = get_born_effective(opts)
    if opts.output: 
        f = open('BORN', 'w')
        f.write('#Born effective charge\n')
        f.write(born.iloc[:, 1:].to_string(index=False, header = False, 
                                      float_format= "%10.6f"))
        f.write('\n')
        f.close
    else: print born

