import numpy  as np
import pandas as pd
import os
import math 
import sys
import JH_lib as jh
from optparse import OptionParser
############################################################
__version__ = "1.1"
# 1.1
# Using Revised POTCAR LIST / Edit Default Path
# Dict for Compound, EnmaxList, ValencyList, AtomicMassList
# Print or not
# 1.0
# Making POTCAR LIST 
# Read POSCAR and Make POTCAR
############################################################

def command_line_arg():
    potcar_path='/team/Process_TCAD/JHLee/c_programs/vasp/POTCAR/PAW_PBE/'
    usage = """
[CODE] READ POSCAR and make new POTCAR 
   usage: %prog [options] arg1 arg2"""

    par = OptionParser(usage=usage, version= __version__)

    par.add_option("--fp", '--folderpath', 
            action='store', type="string", dest='folderpath',
            default='.',
            help='location of the working-space to generate potcar file (Default: ./)') 

    par.add_option("-s", '--structure',  
            action='store', type="string", dest='structure',
            default='./POSCAR',
            help='Define the VASP POSCAR type file in the target folder: Default - ./POSCAR')

    par.add_option("-l", '--list', 
            action='store', type="string", dest='potcar_list',
            default=potcar_path+'potcar_list',
            help='specify potcar for each elements')

    par.add_option("--pp", '--potcarpath', 
            action='store', type="string", dest='potcarpath',
            default='/team/Process_TCAD/JHLee/c_programs/vasp/POTCAR/PAW_PBE',
            help='specify the path of VASP POTCAR (Default - GGA-PBE)')
 
    par.add_option("-p", '--printout', 
            action='store_true', dest='printout',
            default=False,
            help='Print the Code Output (Default: Silent/False')
 
    return  par.parse_args( )

############################################################
def check_potcar(elements):
    f = open(opts.potcar_list, 'r')
    for line in f:
        if "'%s'" %elements in line:
            return line.split()[2:]

if __name__ == '__main__':
    opts, args = command_line_arg()
    if opts.printout: 
        from time import time
        t0 = time()

    unitcell, compound, position = jh.r_cryst_vasp(opts.structure)
    # print(compound[0],position)
    x = jh.w_poscar(position, filename = opts.structure+'.vasp', unitcell = unitcell)
    unitcell, compound, position = jh.r_cryst_vasp(opts.structure+'.vasp')
    os.system('mv %s %s/POSCAR' %( opts.structure+'.vasp' , opts.folderpath))
    
    if opts.printout: print(compound[0],position)
    commends = 'cat ' 
    MetaList = []
    for ele in compound[0]:
        [PotcarName, Enmax, Valency, AtomicMass] = check_potcar(ele)
        MetaDict = {'element':ele, 'PotcarName':PotcarName,  'Enmax':int(Enmax), 'Valency':int(Valency), 'AtomicMass':float(AtomicMass)}
        MetaList.append(MetaDict)
        commends= commends + ' %s/%s/POTCAR ' %(opts.potcarpath, PotcarName)
    commends = commends + '>  %s/POTCAR ' %opts.folderpath
    os.system(commends)
    if opts.printout: 
        print('[CODE] compare automatically generated POTCAR and POSCAR: ',compound[0])
        os.system('grep TIT %s/POTCAR' %opts.folderpath)
    if opts.printout: 
        print(MetaList)

    
    if opts.printout: 
        t1 = time()
        print('\nCode calc completed! Time Used: %.2f [sec]\n' % (t1 - t0))
