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
import pandas as pd
import sys
import JH_lib as jh
from optparse import OptionParser

############################################################
__version__ = "1.0"
# 1.0
# Update for optparse (--input, --draw, --output)
# read REAL DIELECTRIC FUNCTION (LOPTICS = TRUE) 
############################################################


def read_file_line(filename, key):
    i = 0
    with open(filename, 'r') as f:
        for line in f:
            if key in line: 
                k = i
            i += 1
    try:
        return k
    except UnboundLocalError:
        print ("[CODE] ERROR, there is no file")
        return 0

############################################################
############################################################

def get_optical_dielectric(filename):
    start_strings  = "REAL DIELECTRIC FUNCTION (independent particle, no local field effects)"
    end_strings    = "The outermost node in the dielectric function"
    line_num_start = read_file_line(filename, start_strings)
    line_num_end   = read_file_line(filename, end_strings)
    optics_de=[]
    with open(filename, 'r') as f:
        i = 0 
        for line in f:
            if i >  line_num_start + 2 and i < line_num_end -1: 
                optics_de.append(line.split())
            if i == line_num_start + 1 : column_names=line.split()
            i = i + 1
        df = pd.DataFrame(optics_de,columns =column_names, ) 
        for a in df.columns:
            df[a] = pd.to_numeric(df[a])


        try:
            df['mean'] = df.iloc[:,[1,2,3]].mean(axis=1)
        except IndexError:
            print('[CODE] The calculation is not finished (maybe...)')
            df['mean'] = None
            return 0
 #       print df#.dtypes#iloc[:,[1,2,3]]#.mean(axis=1)
    if opts.point or opts.tensor: return df.iloc[0]
    elif opts.lines: return df
    else: return 0



####################################################################
####################################################################

   
def command_line_arg():
    usage = "usage: %prog [options] arg1"  
    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-o", '--outcar', 
            action='store', type="string", dest='outcar',
            default='./OUTCAR',
            help='Path for the OUTCAR, default: .')

    par.add_option("-p", '--point', 
            action='store_true', dest='point',
            default=False,
            help='Print the averaged electric dielectric constant at zero energy\n Default: reading data only for all lines e, xx,yy,zz, mean')


    par.add_option("-t", '--tensor', 
            action='store_true', dest='tensor',
            default=False,
            help='Print the electric dielectric constant at zero energy\n Default: reading data only for all lines e, xx,yy,zz, mean')

    par.add_option("-l", '--lines', 
            action='store_true', dest='lines',
            default=False,
            help='Print the electric dielectric constant values as functio of energy \n Default: reading data only for all lines e, xx,yy,zz, mean')



    return  par.parse_args( )
    
####################################################################
####################################################################



if __name__ == "__main__":
    opts, args = command_line_arg()
    optics_de = get_optical_dielectric(opts.outcar)
    if   opts.point : print (optics_de['mean'])
    elif opts.tensor: print (optics_de)
    elif opts.lines : print (optics_de)

