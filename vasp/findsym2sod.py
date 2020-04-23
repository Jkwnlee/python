#!/sw/util/opt/anaconda3/bin//python3

import pymatgen as pmg
import subprocess
import os
import JH_lib as jh
import numpy as np
from optparse import OptionParser

############################################################
__version__ = "1.0"
# SOD + Findsym + POSCAR
############################################################

def command_line_arg():
    usage = """[CODE] usage: %prog [options] arg1 arg2
[CODE] run this commend with target POSCAR file name and elements which you want:
[CODE] for ex1: python findsym2sod.py -p POSCAR -e Ti N
[CODE]          TO find Ti and N, and build new POSCAR with selected element in the given unitcell
[CODE] for ex2: python /team/ptcad/jhlee/b_codework/py_edit_poscar/expand_poscar.py -p POSCAR -c 0.45 0.73
[CODE]          TO find element between the relative position of 0.45 and 0.73 in z-direction, 
[CODE]          and build new POSCAR with selected element in the given unitcell
[CODE] output: given file name + _target.vasp"""
 
    par = OptionParser(usage=usage, version= __version__)
     
    par.add_option("-p", '--poscar', 
                  action='store',type='string', dest='poscar',
                  default='./POSCAR',
                  help='Choose POSCAR : D-./POSCAR)')

    par.add_option('--path', 
                  action='store',type='string', dest='path',
                  default='.',
                  help='Choose POSCAR : D-.)')

    par.add_option("-f", '--findsym', 
                  action='store',type='string', dest='findsym',
                  default='/team/ptcad/jhlee/b_codework/py_sod_findsym/findsym/bin',
                  help='Define where the findsym bindary is')

    par.add_option("--py_jh",  
                  action='store',type='string', dest='python_jh',
                  default='/team/ptcad/jhlee/b_codework/py_vasp_post_process',
                  help="Define where the jh's python codes in")

    par.add_option("--sod",  
                  action='store',type='string', dest='sod',
                  default='/team/ptcad/jhlee/b_codework/py_sod_findsym/sod/src/sod-0.47/bin/sod_comb.sh',
                  help="Define where the sod codes (sod_comb.sh) codes")

    par.add_option("--target_element", '-t' ,
                  action='store',type='string', dest='target',
                  default=False,
                  help="Define which atom replace in Doping (target in structure)")

    par.add_option("--doping_element", '-d' ,
                  action='store',type='string', dest='doping',
                  default=False,
                  help="Define which atom is doped in structure")

    par.add_option("--n_element", '-n' ,
                  action='store',type='int', dest='Nelement',
                  default=0,
                  help="Define how many atom to dopee")

    par.add_option("-e", '--expansion',
                  action='store',type='string', dest='expand',
                  default='111',
                  help='Doping with cell expansion : Default - 111 , no expansion')

    par.add_option("-o", '--output',
                  action='store',type='string', dest='out',
                  default=False,
                  help='Output file name: filename)')


    return   par.parse_args( )

############################################################


def from_operation_xyz_to_matrix(input_opteration):
    if input_opteration ==  'x':return [ 1, 0, 0]
    if input_opteration == '-x':return [-1, 0, 0]
    if input_opteration ==  'y':return [ 0, 1, 0]
    if input_opteration == '-y':return [ 0,-1, 0]
    if input_opteration ==  'z':return [ 0, 0, 1]
    if input_opteration == '-z':return [ 0, 0,-1]

def from_operation_shift_to_matrix(input_opteration):
    if input_opteration == '0' :    return 0.
    if '/' in input_opteration:
        return float(input_opteration.split('/')[0]) / float( input_opteration.split('/')[1] )
############################################################

def main(opts): 
    print( "[CODE] Start to read the given structure: ", opts.poscar )
    # Read a POSCAR and write to a CIF.
    path1=opts.path
    path2=path1+'/xx_SOD_FINDSYM'
    path3=path2+'/0_pre_SOD_FINDSYM'
    path4=path2+'/1_SOD'
    os.system('mkdir %s' %path2); 
    os.system('mkdir %s' %path3); 
    os.system('mkdir %s' %path4); 
    print("[CODE] Generate the program output in %s", path2)

    repeat  = opts.expand
    print("[CODE] Doping with expanding structure:", repeat)

    unitcell, compound, position = jh.r_cryst_vasp('%s' %opts.poscar)
    structure = pmg.Structure.from_file(opts.poscar)
    structure.to(filename = path3+"/POSCAR.cif")

    # findsym running
    commends ='cd ' + path3 + '; bash '+ opts.findsym +'/findsym.sh %s %s %s' \
    %("./POSCAR.cif", './POSCAR.in', "./POSCAR.iso")
    os.system(commends)
    target=[]
    i = 0
    start = False
    f = open (path3+ "/POSCAR.iso", 'r') 
    for line in f:            
        if start:
    #         print(line)
            if len(line.split()) > 1:
                target.append(line.split()[1].split(',') )
        if 'Space Group' in line:
            head =line
            
        if '_space_group_symop_operation_xyz' in line:
            start = True
        
        if start:
            if 'loop_' in line:
                start = False
    # print(targlet)
    print_out  = False
    f2 = open(path4+ "/SGO", 'w')            
    f2.write('%s' %head)
    if print_out:print('%s' %head,end='')
    for temp in range(len(target)):
        output_m =[]
        output_shift = []
        for temp1 in range(3):
            if len(target[temp][temp1].split('+') ) == 1:
                output_m.append( from_operation_xyz_to_matrix(target[temp][temp1]) )
                output_shift.append('0')
            if len(target[temp][temp1].split('+') ) == 2:
                output_m.append( from_operation_xyz_to_matrix(target[temp][temp1].split('+')[0] ) )
                output_shift.append( target[temp][temp1].split('+')[1])
        f2.write('%1.0i\n' %(temp +1 ))
        if print_out:print( '%1.0i' %(temp +1 ))
        # output_m = [[1, 0, 0], [2, 1, 0], [3, 0, 1]]
        for temp1 in range(3):
            for temp2 in range(3):
                f2.write("%4.1i" %output_m[temp1][temp2])#,end='')
                if print_out:print("%4.1i" %output_m[temp1][temp2],end='')
            f2.write(  "%6.1f\n" %( from_operation_shift_to_matrix(output_shift[temp1]) ))
            if print_out:print(  "%6.1f" %( from_operation_shift_to_matrix(output_shift[temp1]) ))
    f2.write('0\n')
    if print_out:print(0)

    f2.close()



    commends ='python2 %s/poscar_info.py -p %s > %s/out.txt' %(opts.python_jh,opts.poscar,path3)
    os.system(commends)
    element=[]; Nelement  = []
    with open('%s/out.txt'%path3, 'r') as f:
        for line in f:
            if 'C  ' in line:
                element.append(line.split()[1])
                Nelement.append(line.split()[3])
            elif 'l  ' in line:
                if   'vec_a' in line:  len_a = float(line.split()[3])
                elif 'vec_b' in line:  len_b = float(line.split()[3])
                elif 'vec_c' in line:  len_c = float(line.split()[3])
            elif 'a  ' in line:
                if   'alpha' in line:  alpha = float(line.split()[3])
                elif  'beta' in line:  beta  = float(line.split()[3])
                elif 'gamma' in line:  gamma = float(line.split()[3])
            elif 'V  ' in line: 
                if      'vo' in line:  volume= line.split()[3]
    os.system('rm %s/out.txt' %path3)
    # print(element)     
                
    f3 = open('%s/INSOD'%path4, 'w')
    f3.write('''#Title
    Python-Based-Structure-Generate

    # a,b,c,alpha,beta,gamma
    %9.8f  %9.8f  %9.8f  %6.3f  %6.3f  %6.3f\n'''%(len_a,len_b,len_c,alpha,beta,gamma))

    f3.write('''\n# nsp: Number of species
    %i\n'''%len(element))

    f3.write('\n# symbol(1:nsp): Atom symbols\n')
    for ele in element:f3.write('%5s' %ele)

    f3.write('\n\n# natsp0(1:nsp): Number of atoms for each species in the assymetric unit\n')
    for Nele in Nelement:f3.write('%5i' %int(Nele))

    f3.write('\n\n# coords0(1:nat0,1:3): Coordinates of each atom (one line per atom)\n')
    for pos in position:
        #print ('%8.7f %8.7f %8.7f' %(pos[1],pos[2],pos[3]))
        f3.write('%8.7f %8.7f %8.7f\n' %(pos[1],pos[2],pos[3]))

    f3.write('\n# na,nb,nc (supercell definition)\n')
    for r in repeat: f3.write('%3i' %int(r))

    f3.write('\n\n# sptarget: Species to be substituted\n')
    i = 1
    for ele in element:
        if ele == opts.target:
            f3.write('%3i' %i)
        i = i + 1

    f3.write('\n\n# nsubs: Number of substitutions in the supercell\n%3i\n'%opts.Nelement)


    f3.write('''\n# newsymbol(1:2): Symbol of atom to be inserted in the selected position, 
    #                 symbol to be inserted in the rest of the positions for the same species.
    %5s %5s\n''' %(opts.doping, opts.target))
    # Ca Mg 

    f3.write('''
    # FILER, MAPPER
    # # FILER:   0 (no calc files generated), 1 (GULP), 2 (METADISE), 11 (VASP)
    # # MAPPER:  0 (no mapping, use currect structure), >0 (map to structure in MAPTO file)
    # # (each position in old structure is mapped to MAPPER positions in new structure)
    11 0

    ''')
    f3.close()
    os.system('cd ' + path4 + '; bash ' + opts.sod)



############################################################
############################################################

if __name__ == '__main__':
    from time import time
    opts, args = command_line_arg()
    t0 = time()
    main(opts)
    t1 = time()
    
