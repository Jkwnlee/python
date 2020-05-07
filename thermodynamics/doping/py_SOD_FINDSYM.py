import pymatgen as pmg
import subprocess
import os
import numpy as np
from optparse import OptionParser

__version__ = '1.0',
__auther__  = "Ji-Hwan Lee (jihwan.lee@sk.com)",
__summary__ = "Merge Two GNU CODE for doping (SOD and Findsym)"
############################################################

def from_operation_xyz_to_matrix(input_opteration):
    if input_opteration ==  'x':return [ 1, 0, 0]
    if input_opteration == '-x':return [-1, 0, 0]
    if input_opteration ==  'y':return [ 0, 1, 0]
    if input_opteration == '-y':return [ 0,-1, 0]
    if input_opteration ==  'z':return [ 0, 0, 1]
    if input_opteration == '-z':return [ 0, 0,-1]

############################################################

def from_operation_shift_to_matrix(input_opteration):
    if input_opteration == '0' :    return 0.
    if '/' in input_opteration:
        return float(input_opteration.split('/')[0]) / float( input_opteration.split('/')[1] )
    
############################################################
def command_line_arg():
    usage = """
[CODE] Read Log-file from Ginestra License Manager and calculate/display the usage 
    usage: %prog [options] arg1 arg2"""

    par = OptionParser(usage=usage, version= __version__)

    par.add_option("--fp", '--folderpath', 
            action='store', type="string", dest='folderpath',
            default='./',
            help='location of the working-space, containing crystal structure file (Default: ./, ./CONTCAR)')

    par.add_option("-s", '--structure',  
            action='store', type="string", dest='structure',
            default='CONTCAR',
            help='Define the VASP POSCAR type file in the target folder: Default - CONTCAR')

    par.add_option("-a", '--target', 
            action='store', type="string", dest='target_ele',
            default=False,
            help='specify which atoms to substitute')

    par.add_option("-d", '--dopant', 
            action='store', type="string", dest='dopant_ele',
            default=False,
            help='specify which atoms to dope')
 
    par.add_option("-n", '--number',
            action='store', type="float", dest='number',
            default=False,
            help='Number of element to doping in unitcell (considering symmetry)')
    
 
    par.add_option("-r", '--repeat',
            action='store', type="string", dest='repeat',
            default='111',
            help='Cell expansion of unitcell (Default: 111 (meaning 1x1x1)')
    
    par.add_option("--findsympath", 
            action='store', type="string", dest='findsympath',
            default='/team/ptcad/jhlee/b_codework/py_sod_findsym/findsym/bin',
            help='location of findsym binary folder: Default - ~findsym/bin)')

    par.add_option("--sodpath", 
            action='store', type="string", dest='sodpath',
            default='/team/ptcad/jhlee/b_codework/py_sod_findsym/sod/src/sod-0.47/bin/sod_comb.sh',
            help='location of sod_comb.sh folder: Default - ~sod-0.47/bin/sod_comb.sh)')
    
    
    return  par.parse_args( )

############################################################

if __name__ == '__main__':
    from time import time
    opts, args = command_line_arg()

    t0 = time()

    # Read a POSCAR and write to a CIF.
    path2=opts.folderpath+'/xx_SOD_FINDSYM'
    path3=path2+'/0_pre_SOD_FINDSYM'
    path4=path2+'/1_SOD'
    os.system('mkdir %s' %path2); 
    os.system('mkdir %s' %path3); 
    os.system('mkdir %s' %path4); 
    findsympath=opts.findsympath
    
    tar_ele = opts.target_ele
    sub_ele = opts.dopant_ele
    sub_Nele= opts.number
    repeat  = opts.repeat

    structure = pmg.Structure.from_file('%s/%s' %(opts.folderpath,opts.structure))
    structure.to(filename = path3+"/POSCAR.cif")

    # findsym running
    commends ='cd ' + path3 + '; bash '+findsympath+'/findsym.sh %s %s %s' %("./POSCAR.cif", './POSCAR.in', "./POSCAR.iso")
    os.system(commends)
    target=[]
    i = 0
    start =  False
    Natom =  False; element=[]; Nelement  = []
    rposition = False; position = []
    f = open (path3+ "/POSCAR.iso", 'r') 
    for line in f:            
        if start:
            if len(line.split()) > 1:
                target.append(line.split()[1].split(',') )
        if Natom:
            for xx in line.split():
                element.append(xx.split('*')[1])
            Natom = False
            Nelement = np.zeros(np.shape(element))
            
        if rposition:
            if len(line.split()) <2 : rposition = False
            else: 
                position.append(line.split()[1:7])
                for i in range(len(element)):
                    if line.split()[1] == element[i]:
                        Nelement[i] += 1
            
        if 'Space Group' in line: head =line
        elif '_space_group_symop_operation_xyz' in line: start = True
        elif '_cell_length_a'    in line: len_a = float(line.split()[1])
        elif '_cell_length_b'    in line: len_b = float(line.split()[1])
        elif '_cell_length_c'    in line: len_c = float(line.split()[1])
        elif '_cell_angle_alpha' in line: alpha = float(line.split()[1])
        elif '_cell_angle_beta'  in line:  beta = float(line.split()[1])
        elif '_cell_angle_gamma' in line: gamma = float(line.split()[1])
        elif '_cell_volume'      in line: vol   = float(line.split()[1])
        elif 'Type of each atom:' in line: Natom=True
        elif '_atom_site_fract_symmform'in line: rposition = True 

        if start:
            if 'loop_' in line:
                start = False
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
        for temp1 in range(3):
            for temp2 in range(3):
                f2.write("%4.1i" %output_m[temp1][temp2])#,end='')
                if print_out:print("%4.1i" %output_m[temp1][temp2],end='')
            f2.write(  "%6.1f\n" %( from_operation_shift_to_matrix(output_shift[temp1]) ))
            if print_out:print(  "%6.1f" %( from_operation_shift_to_matrix(output_shift[temp1]) ))
    f2.write('0\n')
    if print_out:print(0)

    f2.close()


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
        f3.write('%8.7f %8.7f %8.7f\n' %(float(pos[3]),float(pos[4]),float(pos[5])))

    f3.write('\n# na,nb,nc (supercell definition)\n')
    for r in repeat: f3.write('%3i' %int(r))

    f3.write('\n\n# sptarget: Species to be substituted\n')
    i = 1
    for ele in element:
        if ele == tar_ele:
            f3.write('%3i' %i)
        i = i + 1

    f3.write('\n\n# nsubs: Number of substitutions in the supercell\n%3i\n'%sub_Nele)


    f3.write('''\n# newsymbol(1:2): Symbol of atom to be inserted in the selected position, 
#                 symbol to be inserted in the rest of the positions for the same species.
%5s %5s\n''' %(sub_ele,tar_ele))
    # Ca Mg 

    f3.write('''
# FILER, MAPPER
# # FILER:   0 (no calc files generated), 1 (GULP), 2 (METADISE), 11 (VASP)
# # MAPPER:  0 (no mapping, use currect structure), >0 (map to structure in MAPTO file)
# # (each position in old structure is mapped to MAPPER positions in new structure)
11 0

''')
    f3.close()

    os.system('cd ' + path4 + '; bash ' + opts.sodpath)
    
    t1 = time()
    print('\nCode calc completed! Time Used: %.2f [sec]\n' % (t1 - t0))
