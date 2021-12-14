### Written by Ji-Hwan Lee
### Last update: 2019. 12. 26
### Used library:os, math, time, sys
### Last Update : read non-orthogonal structure in r_cryst_vasp
### The python code for 
### Written by Ji-Hwan Lee


import os
import math
import sys
import time
import numpy as np
import pandas as pd
############################################################
__version__ = "3.0"
#Tuning for python3
#POSITION to DataFrame
############################################################

## ------------------------------------------------------------------------------
"""
LIST-A [FOR STRUCTURES]
1. r_cryst_vasp
: reading VASP crystal structure file type 
(.vasp or called 'POSCAR, CONTCAR')
    1-1. post processing (pst)
        1-1-1. pst_poscar_move
        1-1-2. pst_poscar_ratation_atoms[not yet]
        1-1-3. pst_rot_compar_R[not yet]
        1-1-4. pst_cell_expansion
    1-2. molecule building
        1-2-1. mlc_ch4
        1-2-2. mlc_ch3

2. add-ons (side-tools)
    2-1. poscar_head
         : Not arranged
    2-2. component_from_positions


3. judge_pri_conv
: For some structure, we prefer to have a lattice constant in conventional cell than primitive cell
  In order to do that, we use this function to get conventional cell lattice constant


4. write structure
    4.1 w_poscar(position, compound = None, filename = None, unitcell = None, Selective = False):
    : To write VASP crystal structure file type, this function need to include
    (1) filename
    (2) unitcell=[a, b, c] axis within xyz-coordinate unit of \AA
    (3) compound=[[kind of compound],[number of atoms for each]]
    (4) position=[[kind of compound, position of atom 1 in xyz], [kind of compound, position of atom 2 in xyz]...]
    (5) Selective= True or False 
        This is required to apply Selective Dynamics in slab model or others

    4.2 w_lammps_str(filename, unitcell, compound, position, Selective):
    : To write LAMMPS crystal structure file type, 


5. surface building: 
    5-1. scaling_making_slab(output_name, unitcell_print, initial_position,  \
                        test_x, test_y, diff_atom_position,       \
                        num_fix_lay, num_atom, seq_num, \
                        centering):

    5-2. build_slab(system, metal_kind, lat_con, index, output_name, \
                   vaccume,num_atom,num_fix_lay, symmetric):
            if system == 'fcc': 
                5-2-1. index = '100'
                5-2-2. index = '110'
                5-2-3. index = '111'
                5-2-4. index = '210'
                5-2-5. index = '211'
                5-2-6. index = '311'
                5-2-7. index = '331'
            if system == 'bcc': 'will be updated'

    5-3. add_adsorbate(slab_position, initial_position, mlc_position, output_name)
        # rotation will be updated

6. etc
    6-1 normal_vector(array_for_three_point)
        # Find a normal vector with three point
    6-2 cross_product(vec1,vec2)
        # Calculate cross product for two vector
    6-3 dot_product(vec1,vec2)
        # Calculate dot product for one vector and one matrix (3x3)
    6-4 shift_ads_w_normal_vec(inital_point, distance, vec)
        # Shift the initial point with following vector
    6-5 [up] distance_atoms(a1,a2)
        # Calculate the distance between atoms (two vectors)
    6-6 macroscopic_average(x,y,unitcell, number_of_grid): return macroscopic_average

"""
## ------------------------------------------------------------------------------

def r_cryst_vasp(filename='POSCAR', printout=False):
    def _judge_coord(i,r):
        judge_ = 0.
        len_   = 0.
        for a in range(len(i)):
            len_   = len_ + i[a]**2
            judge_ = judge_ + i[a] * r[a+1]
        judge_ =  judge_ / (len_**0.5)
        if judge_ > 0 and judge_ < 1: return True
        else: return False

    """
    Filename: VASP, POSCAR TYPE
    Function: Read POSCAR type 
    """
    with open(filename, 'r') as f:
        poscar=[]; uc=[]; comp=[]; posi=[];scales=[]
        num_atoms=0; scale_num = 1
        i = 1;  j = 1;  k =0
        Selective=False; kn = 0
        Check_Selective = True
        for l in f:
        #         print(i,l, end='')
            if l == None: break

            if i == 2 : scale_num = float(l)
            elif i > 2 and i < 6: uc.append([float(l_i)*scale_num for l_i in l.split()])
            elif i > 5 and i < 8: comp.append(l.split())
            elif i == 8 : 
                #After Collecting Componenets, Count the total number of atoms
                if num_atoms == 0:  
                    for temp in comp[1]: 
                        num_atoms=num_atoms+int(temp)
#                    print(num_atoms)

                ## Judge the input have Selective Dynamics or Not
                if Check_Selective:
                    if l.split()[0][0] == 'S':
                        Selective= True; i = i - 1
                Check_Selective = False # Pass
                
                if Check_Selective == False:
                # After Checking Selective Dynamics, 
                    if   l.split()[0][0] == 'C': scales=[[1.,0,0],[0,1.,0],[0,0,1.]]
                    elif l.split()[0][0] == 'D' or l.split()[0][0] == 'd': scales=uc

            # Finish reading Unitcell, Composition, Dynamic Conditions
            # Start reading atomic position 
            elif i > 8 and i <= 8 + num_atoms:
                te=[]
                coo = [float(l.split()[l_i]) for l_i in range(3)]
                for svec_num2 in range(3):
                    scaled_coo = 0.0
                    for svec_num1 in range(3):
                    #        print(scales, svec_num1,svec_num2)
                        scaled_coo = scaled_coo +  scales[svec_num1][svec_num2] *coo[svec_num1]
                    te.append(scaled_coo)
                [x, y, z] = te
                ## k = k-th atom , kn = Criteria for Composition
                if int(comp[1][kn]) == 0:
                    # if there is number of atom  ==  0, skip it [bug-fix 206002]
                    kn = kn + 1
                    k = 0 # Initialize the K (reset total number of k-th element in comp[0])

                if k <= int(comp[1][kn]):
                    if Selective: posi.append([comp[0][kn], x, y, z, l.split()[3], l.split()[4], l.split()[5]])
                    else:         posi.append([comp[0][kn], x, y, z])
                k= k+1

                if k == int(comp[1][kn]):
                    kn = kn + 1
                    k = 0 # Initialize the K (reset total number of k-th element in comp[0])
            
            ##After reading whole informations
            if i == 8 + num_atoms:
                #judge the position (r) is in the unitcell (a,b,c) : 0<dot(i,r)/|i|< 1
                for a in posi:
                    for x in range(3):
                        if _judge_coord(uc[x], a): pass
                        else: #move with consideration of periodic boundary condition
                            a = a + uc[x]
                if Selective:df = pd.DataFrame( posi, columns=['element', 'x', 'y', 'z','rx','ry','rz'])
                else: df = pd.DataFrame( posi, columns=['element', 'x', 'y', 'z'])
                
                df = df.sort_values(by=['element','z'])
                if printout:
                    print("[CODE] The atomic postions are sorted by alphabet, Please revise your POTCAR")
                posi = df.values.tolist()
#                for tetete in posi: print(tetete)
                return uc, comp, posi
            elif i > 8 + num_atoms: break
            else:  i = i + 1

## ------------------------------------------------------------------------------
def pst_poscar_move(position, mv_x, mv_y,mv_z):
    new_position = position
    for temp in range(len(position)):
        new_position[temp][1]=position[temp][1] - mv_x
        new_position[temp][2]=position[temp][2] - mv_y
        new_position[temp][3]=position[temp][3] - mv_z
    return new_position

## ------------------------------------------------------------------------------


def pst_cell_expansion(unitcell, position, dimension):
    """
    1-1-4. pst_cell_expansion

    To make a supercell, expand the cell and add atoms properly
    output datas are new_unitcell, new_compound, new_position
    """
    new_position=[]; new_unitcell=[[0,0,0],[0,0,0],[0,0,0]]; new_compound=[]; temp=[]
    for temp1 in range(len(position)):
        ## expanding to a
        for temp2 in range(dimension[0]): 
            ## expanding to b
            for temp3 in range(dimension[1]):
                ## expanding to c
                for temp4 in range(dimension[2]):
                    if len(position[temp1]) == 4: 
                        new_position.append(\
                            [position[temp1][0],\
                             position[temp1][1] + temp2 * float(unitcell[0][0]) + temp3 * float(unitcell[1][0]) + temp4 * float(unitcell[2][0]) , \
                             position[temp1][2] + temp2 * float(unitcell[0][1]) + temp3 * float(unitcell[1][1]) + temp4 * float(unitcell[2][1]) , \
                             position[temp1][3] + temp2 * float(unitcell[0][2]) + temp3 * float(unitcell[1][2]) + temp4 * float(unitcell[2][2]) , ])
                    if len(position[temp1]) == 7:
                        new_position.append(\
                            [position[temp1][0],\
                             position[temp1][1] + temp2 * float(unitcell[0][0]) + temp3 * float(unitcell[1][0]) + temp4 * float(unitcell[2][0]) , \
                             position[temp1][2] + temp2 * float(unitcell[0][1]) + temp3 * float(unitcell[1][1]) + temp4 * float(unitcell[2][1]) , \
                             position[temp1][3] + temp2 * float(unitcell[0][2]) + temp3 * float(unitcell[1][2]) + temp4 * float(unitcell[2][2]) , \
                             position[temp1][4], position[temp1][5], position[temp1][6]])

    for a in unitcell: temp.append(a)
    for temp1 in range(len(temp)):
        # each axis
        for temp2 in range(len(temp[temp1])):
            #each component
            new_unitcell[temp1][temp2] = float(temp[temp1][temp2]) * dimension[temp1]

    ## Count Component from positions
    new_compound,new_position = component_from_positions(new_position)

    return new_unitcell, new_compound, new_position


## ------------------------------------------------------------------------------


def mlc_ch4(initial_position_C=None, bond_leng=None):
    """
    1-2-1. mlc_ch4: building \ce{CH4} molecule 
    : As first varibalbe, define the position of C as like [0,0,0] 
      and 
      (Optional) insert the bonding length as second variable
                 [default is  1.09601 angstrom]

    : rotation will be [TODO:updated]

    """
    compound = [['C', 'H'], [1, 4]]
    if bond_leng == None:
        bond_leng = 1.09601
    if initial_position_C == None:
        initial_position_C = [0, 0, 0]

    position=[['C',                          0. ,                       0. ,              0. , 'T', 'T', 'T' ],\
              ['H',                          0. ,                       0. , -     bond_leng , 'T', 'T', 'T' ],\
              ['H',  8.**(0.5) * bond_leng / 3. ,                       0. ,   bond_leng / 3 , 'T', 'T', 'T' ],\
              ['H', bond_leng * -2.**(0.5) / 3. ,  bond_leng * 2 / 6**(0.5),   bond_leng / 3 , 'T', 'T', 'T' ],\
              ['H', bond_leng * -2.**(0.5) / 3. , -bond_leng * 2 / 6**(0.5),   bond_leng / 3 , 'T', 'T', 'T' ] ]
    position = pst_poscar_move(position, - initial_position_C[0], - initial_position_C[1], - initial_position_C[2]) 

    return position, compound


## ------------------------------------------------------------------------------


def mlc_ch3(initial_position_C=None, bond_leng=None):
    """
    1-2-2. mlc_ch3: building \ce{CH3} molecule 
    : As first varibalbe, define the position of C as like [0,0,0] 
      and 
      (Optional) insert the bonding length as second variable
                 [default is  1.09601 angstrom]

    : rotation will be [TODO:updated]

    """
    compound = [['C', 'H'], [1, 3]]
    if bond_leng == None:
        bond_leng = 1.09601
    if initial_position_C == None:
        initial_position_C = [0, 0, 0]

    position=[['C',                          0. ,                       0. ,              0. , 'T', 'T', 'T' ],\
              ['H',  8.**(0.5) * bond_leng / 3. ,                       0. ,   bond_leng / 3 , 'T', 'T', 'T' ],\
              ['H', bond_leng * -2.**(0.5) / 3. ,  bond_leng * 2 / 6**(0.5),   bond_leng / 3 , 'T', 'T', 'T' ],\
              ['H', bond_leng * -2.**(0.5) / 3. , -bond_leng * 2 / 6**(0.5),   bond_leng / 3 , 'T', 'T', 'T' ] ]
    position = pst_poscar_move(position, - initial_position_C[0], - initial_position_C[1], - initial_position_C[2]) 

    return position, compound


## ------------------------------------------------------------------------------


def poscar_head(num_atom,num_rlx):
    scale=( num_atom - 1 ) / 2
    num_fix_lay=( num_atom - num_rlx ) / 2 


## ------------------------------------------------------------------------------


def component_from_positions(insert_position, printout=False):
    if len(insert_position[0]) > 4:Selective = True
    else: Selective = False
    
    if Selective:
    	df = pd.DataFrame( insert_position, columns=['element', 'x', 'y', 'z','rx','ry','rz'])
    else: 
    	df = pd.DataFrame( insert_position, columns=['element', 'x', 'y', 'z'])
    df = df.sort_values(by='element')
    new_compound=[]
    new_compound.append(list(df.groupby(by='element').agg('count').index))
    new_compound.append(df.groupby(by='element').agg('count').loc[:,'x'].tolist()) 


#    new_compound=[]
#    componets = []
#    num_componets = []
#    i = 0; k= 0
#    for temp1 in range(len(insert_position)):
#        if temp1 == 0:
#            componets.append(insert_position[temp1][0])
#            num_componets.append(1)
#        else:
#            if insert_position[temp1][0] == temp:
#                if componets[i] == insert_position[temp1][0]:
#                    num_componets[i] = num_componets[i] + 1
#            else:
#                i = i + 1
#                componets.append(insert_position[temp1][0])
#                num_componets.append(1)
#
#        temp = insert_position[temp1][0]
#
#    new_compound = [componets , num_componets]
    if printout:
        print("[CODE] The atomic postions are sorted by alphabet, Please revise your POTCAR")
#    print new_compound
    return new_compound,df.values.tolist()
    

## ------------------------------------------------------------------------------


def judge_pri_conv(filename, system):
    ## Function check the condition of file wheather it is primitive or not
    ##
    ## Condition1: number atom == 1
    ## Condition2: nit
    unitcell, compound, position=r_cryst_vasp('CONTCAR')
    if len(compound[1]) == 1:
        # print("this cell include only one kind of atom")
        if compound[1][0] == '1':
            # print("this cell include only one of atom")
            True
    else:
        print("this cell include only more than one kind of atoms", compound[1])
        pass
    if system == 'fcc':
        material_name=compound[0][0]
        latcon=unitcell[0][0]
        if float(latcon) == 0.:
            # print("the unit cell is not conventional")
            latcon=float(unitcell[0][1]) * 2
        else:
            pass
    return latcon


## ------------------------------------------------------------------------------


def w_poscar(position, compound = None, filename = None, unitcell = None, Selective = False, printout=False):
    """
    4.1 Write a poscar file type
    Officially, it needs 2 parameters. "compound and position"

    [default] filename is 'py_POSCAR.vasp'
    [default] unitcell is a cubic with a lattice, 1.1 times bigger than maximum in "position"
    [default] Selective == None, if you want to turn on, 
              Please type True with proper "position" including Selective information inside
              ( ex. [[Cu, 0, 0, 0, F, F, F]]  )
    """
    if compound == None:
        if printout:
            print('[CODE] There is no initial compound data:%s' %compound)
            print('[CODE] Check the output once again')
        compound,position = component_from_positions(position)
    compound,position = component_from_positions(position)
   # for i in position: print(i)
    
    if filename == None:
        filename = 'py_POSCAR.vasp'
    

    if unitcell == None:
        tempx = 10.
        tempy = 10.
        tempz = 10.
        for temp in position:
            tempx = max(abs(float(temp[1])), tempx)
            tempy = max(abs(float(temp[2])), tempy)
            tempz = max(abs(float(temp[3])), tempz)
        unitcell = [[tempx*1.1, 0., 0.], \
                    [ 0.,tempy*1.1, 0.], \
                    [ 0., 0.,tempz*1.1]]
                    
#    print len(position[0]),position[0]
    if len(position[0]) == 7:
        if position[0][6] == 'T' or position[0][6] == 'F': Selective = True


    f=open(filename, 'w')
    f.write('test\n %19.16f\n' % 1.00000000000000)
    for temp in range(len(unitcell)):
        xyz=unitcell[temp]
        if len(xyz) == 3:
            f.write(' %22.16f%22.16f%22.16f \n' % (float(xyz[0]), float(xyz[1]), float(xyz[2]))) 
        else: 
            print(" WARNING the xyz is not proper 'unitcell' SHOULD CONTAIN ONLY 3 NUMBERS")
            break

    for temp in range(len(compound)):
        for temp2 in range(len(compound[temp])):
            f.write('  %s' %compound[temp][temp2])
#            print 'test',compound[temp][temp2]
        f.write('\n')
        
    if Selective:
        f.write('Selective Dynamics\n')
        
    f.write('Cartesian\n')
    for temp in range(len(position)):
        xyz=position[temp][1:4]
        if len(xyz) == 3:
            if Selective:
                if len(position[temp]) ==  4:
                    dynamic = ['T','T','T']
                else:
                    dynamic=position[temp][4:]
                f.write(' %22.16f%22.16f%22.16f %5.1s %5.1s %5.1s \n' \
                    % (float(xyz[0]), float(xyz[1]), float(xyz[2]), dynamic[0],dynamic[1], dynamic[2])) 
            else:
                f.write(' %22.16f%22.16f%22.16f \n' % (float(xyz[0]), float(xyz[1]), float(xyz[2]))) 
        else: 
            print(" WARNING the iserted position is not proper. IT SHOULD CONTAIN ONLY 3 NUMBERS")
            print( xyz)
            break



## ------------------------------------------------------------------------------


def w_lammps_str(compound, position, filename = None, unitcell = None, Selective = None):
    """
    4.2 w_lammps_str(filename, unitcell, compound, position, Selective):
    : To write LAMMPS crystal structure file type, 
    Officially, it needs 2 parameters. "compound and position"

    [default] filename is 'py_POSCAR.vasp'
    [default] unitcell is a cubic with a lattice, 1.1 times bigger than maximum in "position"
    [default] Selective == None, if you want to turn on, 
              Please type True with proper "position" including Selective information inside
              ( ex. [[Cu, 0, 0, 0, F, F, F]]  )
    """
    if filename == None:
        filename = 'py_POSCAR.vasp'

    if unitcell == None:
        tempx = 10.
        tempy = 10.
        tempz = 10.
        for temp in position:
            tempx = max(abs(float(temp[1])), tempx)
            tempy = max(abs(float(temp[2])), tempy)
            tempz = max(abs(float(temp[3])), tempz)
        unitcell = [[tempx*1.1, 0., 0.], \
                    [ 0.,tempy*1.1, 0.], \
                    [ 0., 0.,tempz*1.1]]

    if len(position[1]) == 7:
        Selective = True


    f=open(filename, 'w')
    f.write('test\n %19.16f\n' % 1.00000000000000)
    for temp in range(len(unitcell)):
        xyz=unitcell[temp]
        if len(xyz) == 3:
            f.write(' %22.16f%22.16f%22.16f \n' % (float(xyz[0]), float(xyz[1]), float(xyz[2]))) 
        else: 
            print(" WARNING the xyz is not proper IT SHOULD CONTAIN ONLY 3 NUMBERS")
            break

    for temp in range(len(compound)):
        for temp2 in range(len(compound[temp])):
            f.write(' %4.3s' %compound[temp][temp2])
        f.write('\n')
    if Selective == None:
        pass
    else:
        f.write('Selective Dynamics\n')
    f.write('C\n')

    for temp in range(len(position)):
        xyz=position[temp][1:4]
        if len(xyz) == 3:
            if Selective:
                if len(position[temp]) == 4:
                    dynamic=['T', 'T', 'T']
                else:
                    dynamic=position[temp][4:]
                f.write(' %22.16f%22.16f%22.16f %5.1s %5.1s %5.1s \n' \
                    % (float(xyz[0]), float(xyz[1]), float(xyz[2]), dynamic[0],dynamic[1], dynamic[2])) 
            else:
                f.write(' %22.16f%22.16f%22.16f \n' % (float(xyz[0]), float(xyz[1]), float(xyz[2]))) 
        else: 
            print(" WARNING the xyz is not proper IT SHOULD CONTAIN ONLY 3 NUMBERS")
            break


## ------------------------------------------------------------------------------

def scaling_making_slab(output_name, unitcell_print, initial_position,  \
                        test_x, test_y, diff_atom_position,       \
                        num_fix_lay, num_atom, seq_num, \
                        centering):
    ## store positions into positions
    position_temp=initial_position[0][1:4]
    position_after=position_temp
    position_print=[]
    for temp in range(num_atom):
        position_print.append(['Cu',position_temp[0], position_temp[1], position_temp[2]])
        ## Define x (a) coefficient!
        position_after[0] = position_temp[0] + diff_atom_position[0]    
        if position_after[0] >= test_x: position_after[0]= position_after[0] - test_x
        ## Define y (b) coefficient!
        position_after[1] = position_temp[1] + diff_atom_position[1]    
        if position_after[1] >= test_y: position_after[1]= position_after[1] - test_y
        ## Define z (c) coefficient!
        position_after[2] = position_temp[2] + diff_atom_position[2]
        position_temp=position_after

    if centering:
        z_collect=[]
        temp_num_rlx=1 ; temp_num_unrlx=1
        for temp in range(len(position_print)):
            z_collect.append(position_print[temp][3])
        averaged=reduce(lambda x, y: x + y, z_collect) / len(z_collect)
        for temp in range(len(position_print)):
            position_print[temp][3]=position_print[temp][3]-averaged
            if ( num_atom - num_fix_lay ) / 2.  >=  temp_num_rlx :
                position_print[temp].append('T');position_print[temp].append('T');position_print[temp].append('T')
                temp_num_rlx=temp_num_rlx + 1
            else:
                if temp_num_unrlx > num_fix_lay:
                    position_print[temp].append('T');position_print[temp].append('T');position_print[temp].append('T')
                else:
                    position_print[temp].append('F');position_print[temp].append('F');position_print[temp].append('F')
                    temp_num_unrlx = temp_num_unrlx + 1
    else:
        for temp in range(num_atom):
            if temp <=  num_fix_lay :
                position_print[temp].append('F');position_print[temp].append('F');position_print[temp].append('F')
            else:
                position_print[temp].append('T');position_print[temp].append('T');position_print[temp].append('T')

    compound_print=[[position_print[0][0]],[len(position_print)]]
    w_poscar(position_print, compound_print, output_name, unitcell_print, True)


## ------------------------------------------------------------------------------


def build_slab(system, metal_kind, lat_con, index, output_name, vaccume,num_atom,num_fix_lay, symmetric):
    if system == 'fcc':
        ## CHOOSE face-centered crystal
        position=[[metal_kind, 0.0, 0.0, 0.0, 'F', 'F', 'F']]
        if index == '100':
            unitcell_print, diff_atom_position = poscar_fcc_100_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=2
        elif index == '110':
            unitcell_print, diff_atom_position = poscar_fcc_110_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=2
        elif index == '111':
            unitcell_print, diff_atom_position = poscar_fcc_111_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=3
        elif index == '210':
            unitcell_print, diff_atom_position = poscar_fcc_210_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=10
        elif index == '211':
            unitcell_print, diff_atom_position = poscar_fcc_211_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=6
        elif index == '311':
            unitcell_print, diff_atom_position = poscar_fcc_311_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=11
        elif index == '331':
            unitcell_print, diff_atom_position = poscar_fcc_331_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=19

    test_x  = unitcell_print[0][0] + unitcell_print[1][0] + unitcell_print[2][0]
    test_y  = unitcell_print[0][1] + unitcell_print[1][1] + unitcell_print[2][1]

    scaling_making_slab(output_name, unitcell_print, position, \
        test_x ,test_y ,diff_atom_position,\
        num_fix_lay,  num_atom, seq_num, symmetric)


## ------------------------------------------------------------------------------


def poscar_fcc_100_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con / 2.**(0.5) ;vec_ay = 0.0 ;vec_az = 0.0
    vec_bx = 0.0         ;vec_by = lat_con / 2.**(0.5) ;vec_bz = 0.0
    vec_cx = 0.0         ;vec_cy = 0.0
    vec_cz = lat_con / 2 * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]
    ## Position of next atom 
    len_x = vec_ax / 2.  ; len_y= vec_by / 2.  ;len_z = lat_con / 2.
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position     


## ------------------------------------------------------------------------------


def poscar_fcc_110_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con         ;vec_ay = 0.0         ;vec_az = 0.0
    vec_bx = 0.0             ;vec_by = lat_con / 2.**(0.5)  ;vec_bz = 0.0
    vec_cx = 0.0             ;vec_cy = 0.0
    vec_cz = lat_con / 2.**(0.5) / 2 * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]
    ## Position of next atom 
    len_x = vec_ax  / 2.    ;len_y = vec_by  / 2.  ;len_z = lat_con / 2.**(0.5) / 2
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position


## ------------------------------------------------------------------------------


def poscar_fcc_111_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con / 2.**(0.5)  ;vec_ay = 0.0     ;vec_az = 0.0
    vec_bx = lat_con / 2.**(0.5) / 2.
    vec_by = lat_con / 2.**(0.5) / 2. * 3. **(0.5) ;vec_bz = 0.0
    vec_cx = 0.0                  ;vec_cy = 0.0
    vec_cz = lat_con * 3. **(0.5) / 3 * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]

    ## Position of next atom 
    len_x = ( vec_bx + vec_ax) / 3. ; len_y = vec_by / 3. ;len_z = lat_con * 3. **(0.5) / 3
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position


## ------------------------------------------------------------------------------


def poscar_fcc_211_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con * 2.**(0.5) / 2. ;vec_ay = 0.0        ;vec_az = 0.0
    vec_bx = 0.0         ;vec_by = lat_con * 3.**(0.5)     ;vec_bz = 0.0
    vec_cx = 0.0         ;vec_cy = 0.0
    vec_cz = lat_con * 6. **(0.5) / 2 / 6 * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]

    ## Position of next atom 
    len_x = vec_ax / 2.  ;len_y = vec_by * 2. / 3.         ;len_z = lat_con * 6. **(0.5) / 2. / 6.
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position


## ------------------------------------------------------------------------------


def poscar_fcc_311_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con * 6.**(0.5) / 2. ;vec_ay = 0.0        ;vec_az = 0.0
    vec_bx = lat_con * 6.**(0.5) / 2. * 5 / 6. 
    vec_by = lat_con * 6.**(0.5) / 2  * 11.**(0.5) / 6     ;vec_bz = 0.0
    vec_cx = 0.0         ;vec_cy = 0.0
    vec_cz = lat_con / 11. **(0.5) * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]

    ## Position of next atom  (CONTROVERSY for 0.2727)
    len_x = (vec_ax + vec_bx ) * 0.2727  ;len_y = vec_by * 0.2727  ;len_z = lat_con / 11. **(0.5)
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position


## ------------------------------------------------------------------------------


def poscar_fcc_210_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con * 6.**(0.5) / 2. ;vec_ay = 0.0        ;vec_az = 0.0
    vec_bx = vec_ax * 2. / 3.  ;vec_by = vec_ax * 5.**(0.5) / 3  ;vec_bz = 0.0
    vec_cx = 0.0         ;vec_cy = 0.0
    vec_cz = lat_con * 5.**(0.5) / 10. * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]

    ## Position of next atom  (CONTROVERSY for 0.2727)
    len_x = (vec_ax + vec_bx ) * 0.7  ;len_y = vec_by * 0.7  ;len_z = lat_con * 5. **(0.5) / 10
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position


## ------------------------------------------------------------------------------


def poscar_fcc_331_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con * 10.**(0.5) / 2. ;vec_ay = 0.0        ;vec_az = 0.0
    vec_bx = lat_con * 10.**(0.5) / 2. * 0.9
    vec_by = lat_con * 10.**(0.5) / 2. * (0.19)**(0.5)      ;vec_bz = 0.0
    vec_cx = 0.0         ;vec_cy = 0.0
    vec_cz = lat_con / 19.**(0.5)  * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]

    ## Position of next atom  (CONTROVERSY for 0.2727)
    len_x = (vec_ax + vec_bx ) * (1 - 10.**(0.5) / 10 )
    len_y = vec_by * (1 - 10.**(0.5) / 10 )  ;len_z = lat_con / 19. **(0.5)
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position


## ------------------------------------------------------------------------------


def add_adsorbate(slab_position, slab_unitcell, adsorbate_position, output_name):
    new_position = slab_position + adsorbate_position
    new_compound,new_position = component_from_positions(new_position)
    # print new_position, new_compound
    w_poscar(new_position, compound = new_compound, filename = output_name, unitcell = slab_unitcell, Selective = True)


## ------------------------------------------------------------------------------


def print_error(judge = True, message_from=False):
    if judge:
        print("""
     _____ ____  ____   ___  ____
    | ____|  _ \|  _ \ / _ \|  _ \\
    |  _| | |_) | |_) | | | | |_) |
    | |___|  _ <|  _ <| |_| |  _ <
    |_____|_| \_\_| \_\\\___/|_| \_\\
    """)
    if message_from != False:
        print("Error Message From the function of ", message_from)


## ------------------------------------------------------------------------------


def normal_vector(array_for_three_point):
    """
    Calculate the normal vector with a inserted three points
    """

    if len(array_for_three_point) != 3:
        return print_error(judge = True, message_from='normal_vector')
    point0 = array_for_three_point[0]; point1 = array_for_three_point[1];point2 = array_for_three_point[2]
    vec1=[]
    vec2=[]
    center_of_points=[]
    for a in range(3):
        vec1.append(point0[a] - point1[a])
        vec2.append(point0[a] - point2[a])
        center_of_points.append( (point0[a]+ point1[a] + point2[a]) / 3. )
    # print vec1
    # print vec2
    temp= []
    if cross_product(vec1, vec2)[2] < 0:
        temp = cross_product(vec2, vec1)
        for a in range(3):
            temp[a] = temp[a]/math.sqrt(temp[0]**2 + temp[1]**2 + temp[2]**2)        
        
    else:
        temp = cross_product(vec1, vec2)
        for a in range(3):
            temp[a] = temp[a]/math.sqrt(temp[0]**2 + temp[1]**2 + temp[2]**2)
    return temp, center_of_points
    # if cross_vec12


## ------------------------------------------------------------------------------


def cross_product(vec1, vec2):
    """
    Calculate the corss product from vec1 to vec2
    """
    if len(vec1) != 3 or len(vec2) != 3:
        return print_error(judge = True, message_from='cross_product')
    else:
        return [vec1[1] * vec2[2]  - vec1[2] * vec2[1], 
                vec1[2] * vec2[0]  - vec1[0] * vec2[2],
                vec1[0] * vec2[1]  - vec1[1] * vec2[0]]


## ------------------------------------------------------------------------------


def shift_ads_w_normal_vec(inital_point,distance, vec):
    """
    Shift the initial point with following vector
    """
    if len(inital_point) !=3 or len(vec) !=3:
        return print_error(judge = True, message_from='shift_ads_w_normal_vec')
        
    len_vec = math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
    temp_vec = []; next_position=[]
    for temp in range(3):
        temp_vec.append(vec[temp]/len_vec * distance)
        next_position.append(inital_point[temp] + temp_vec[temp])
    # print(' %22.16f%22.16f%22.16f T  T  T\n' % (float(next_position[0]), float(next_position[1]), float(next_position[2]))) 

    return temp_vec,next_position


## ------------------------------------------------------------------------------


def Rot(theta,r_axis):
    ## r_axis: The following three basic rotation matrices rotate vectors by
    ## an angle "theta" about the x-, y-, or z-axis, in three dimensions, using the right-hand rule

    ## CHEK THETA IS FLOAT
    co=math.cos(theta / 180. * math.pi)
    si=math.sin(theta / 180. * math.pi)

    if r_axis == 'x' or 'a':
        v = [1, 0, 0]
    elif r_axis == 'y' or 'b':
        v = [0, 1, 0]
    elif r_axis == 'z' or 'c':
        v = [0, 0, 1]
    elif len(r_axis) == 3:
        pass
    else:
        return print_error(judge = True)
    
    return [[co+(v[0]**2)*(1-co),v[0]*v[1]*(1-co)-v[2]*si,v[0]*v[2]*(1-co)+v[1]*si], 
            [v[1]*v[0]*(1-co)+v[2]*si,co+(v[1]**2)*(1-co),v[1]*v[2]*(1-co)-v[0]*si], 
            [v[2]*v[0]*(1-co)-v[1]*si,v[2]*v[1]*(1-co)+v[0]*si,co+(v[2]**2)*(1-co)]] 


## ------------------------------------------------------------------------------


def dot_product(x, y):
    # x: (3x3) matrix
    # y: point
    return [ x[0][0] * y[0] + x[0][1] * y[1] + x[0][2] * y[2], 
        x[1][0] * y[0] + x[1][1] * y[1] + x[1][2] * y[2],
        x[2][0] * y[0] + x[2][1] * y[1] + x[2][2] * y[2] ]


## ------------------------------------------------------------------------------


def rotate_mlc(input_data, angle, R_axis, output_name=False):
    in_compound = []
    if type(input_data) is str:
        in_unitcell, in_compound, in_position = r_cryst_vasp(input_data)
        input_name = input_data
    elif type(input_data) is list:
        in_position = input_data
        input_name = 'None'
    else:
        print_error(judge = True, message_from='rotate_mlc')

    x=0; y=0; z=0; i = 0
    for a in in_position:
        x = x + float(a[1])
        y = y + float(a[2])
        z = z + float(a[3])
        i = i + 1
    center=[x/i, y/i, z/i]

    new_position=[]
    R_matrix = Rot(angle, R_axis)
    position = pst_poscar_move(in_position, center[0], center[1], center[2])
    for a in in_position:
        temp=[]
        temp=dot_product(R_matrix , a[1:])
        new_position.append([a[0],temp[0], temp[1], temp[2]])


    new_position = pst_poscar_move(new_position, -center[0], -center[1], -center[2])
    if in_compound == None:
        new_compound,new_position = component_from_positions(new_position)
    else:
        new_compound = in_compound


    if output_name == False:
        pass
    elif output_name ==  True:
        output_name1 = "R_"+str(angle)+'_'+input_name
        w_poscar(new_position, new_compound, output_name1)
    else:
        w_poscar(new_position, new_compound, output_name)

    
    return new_compound, new_position


## ------------------------------------------------------------------------------


def w_lammps_str(compound, position, filename = None, unitcell = None, sequence_of_compound = None):
    """
    4.2 w_lammps_str(filename, unitcell, compound, position, Selective):
    : To write LAMMPS crystal structure file type, 
    Officially, it needs 2 parameters. "compound and position"

    [default] filename is 'py_POSCAR.vasp'
    [default] unitcell is a cubic with a lattice, 1.1 times bigger than maximum in "position"
    [default] Selective == None, if you want to turn on, 
              Please type True with proper "position" including Selective information inside
              ( ex. [[Cu, 0, 0, 0, F, F, F]]  )
    """
    if filename == None:
        filename = 'py_structure.data'

    if unitcell == None:
        tempx = 10.
        tempy = 10.
        tempz = 10.
        for temp in position:
            tempx = max(abs(float(temp[1])), tempx)
            tempy = max(abs(float(temp[2])), tempy)
            tempz = max(abs(float(temp[3])), tempz)
        unitcell = [[tempx*1.1, 0., 0.], \
                    [ 0.,tempy*1.1, 0.], \
                    [ 0., 0.,tempz*1.1]]

    if sequence_of_compound == None:
        sequence_of_compound = compound[0]  ## TODO later
    n_atom=0
    for a in compound[1]:
        n_atom=n_atom + int(a)

    f=open(filename, 'w')
    f.write('# LAMMPS data file created from JH\'s python\n')
    f.write('%5.0f atoms\n' %(n_atom ) )
    f.write('%3.0f atom types\n' %(len(compound[0]) ) )
    f.write('%10.6f %8.6f  xlo xhi\n' %( 0, float(unitcell[0][0]) ) )
    f.write('%10.6f %8.6f  ylo yhi\n' %( 0, float(unitcell[1][1]) ) )
    f.write('%10.6f %8.6f  zlo zhi\n' %( 0, float(unitcell[2][2]) ) )
    f.write('%10.6f %8.6f %8.6f  xy xz yz\n\nMasses\n\n' %( 0,0,0  ) )

    compound = mass_of_element(compound)
    i = 1
    for a in sequence_of_compound[0]:
        for b in range(len(compound[0])):
            if a == compound[0][b]:
                f.write('%1.0f %8.6f\n' %(i, float(compound[2][b])))
                i = i + 1
    f.write('\nAtoms\n\n')

    i = 1
    for a in range(len(sequence_of_compound[0])):    
        for temp in position:
            if sequence_of_compound[0][a] == temp[0]:
                f.write('%1.0f %1.0f %2.2s %8.6f %8.6f %8.6f\n' \
                         %(i, a+1, sequence_of_compound[1][a], float(temp[1]), float(temp[2]), float(temp[3]) ) )
                i = i + 1


## ------------------------------------------------------------------------------

def mass_of_element(compound):
    lib=[   [  'H',  1.00], [ 'He',  4.00], [ 'Li',  6.94], [ 'Be',  9.01], [  'B', 10.81], [  'C', 12.01], [  'N', 14.00],  \
            [  'O', 15.99], [  'F', 18.99], [ 'Ne', 20.17], [ 'Na', 22.98], [ 'Mg', 24.30], [ 'Al', 26.98], [ 'Si', 28.08],  \
            [  'P', 30.97], [  'S', 32.06], [ 'Cl', 35.45], [  'K', 39.09], [ 'Ar', 39.94], [ 'Ca', 40.08], [ 'Sc', 44.95],  \
            [ 'Ti', 47.90], [  'V', 50.94], [ 'Cr', 51.99], [ 'Mn', 54.93], [ 'Fe', 55.84], [ 'Ni', 58.70], [ 'Co', 58.93],  \
            [ 'Cu', 63.54], [ 'Zn', 65.38], [ 'Ga', 69.72], [ 'Ge', 72.59], [ 'As', 74.92], [ 'Se', 78.96], [ 'Br', 79.90],  \
            [ 'Kr', 83.80], [ 'Rb', 85.46], [ 'Sr', 87.62], [  'Y', 88.90], [ 'Zr', 91.22], [ 'Nb', 92.90], [ 'Mo', 95.94],  \
            [ 'Tc', 98.00], [ 'Ru',101.07], [ 'Rh',102.90], [ 'Pd',106.40], [ 'Ag',107.86], [ 'Cd',112.41], [ 'In',114.82],  \
            [ 'Sn',118.69], [ 'Sb',121.75], [  'I',126.90], [ 'Te',127.60], [ 'Xe',131.30], [ 'Cs',132.90], [ 'Ba',137.33],  \
            [ 'La',138.90], [ 'Ce',140.12], [ 'Pr',140.90], [ 'Nd',144.24], [ 'Pm',145.00], [ 'Sm',150.40], [ 'Eu',151.96],  \
            [ 'Gd',157.25], [ 'Tb',158.92], [ 'Dy',162.50], [ 'Ho',164.93], [ 'Er',167.26], [ 'Tm',168.93], [ 'Yb',173.04],  \
            [ 'Lu',174.96], [ 'Hf',178.49], [ 'Ta',180.94], [  'W',183.85], [ 'Re',186.20], [ 'Os',190.20], [ 'Ir',192.22],  \
            [ 'Pt',195.09], [ 'Au',196.96], [ 'Hg',200.59], [ 'Tl',204.37], [ 'Pb',207.20], [ 'Bi',208.98], [ 'Po',209.00],  \
            [ 'At',210.00], [ 'Rn',222.00], [ 'Fr',223.00], [ 'Ra',226.02], [ 'Ac',227.02], [ 'Pa',231.03], [ 'Th',232.03],  \
            [ 'Np',237.04], [  'U',238.02], [ 'Pu',242.00], [ 'Am',243.00], [ 'Bk',247.00], [ 'Cm',247.00], [ 'No',250.00],  \
            [ 'Cf',251.00], [ 'Es',252.00], [ 'Hs',255.00], [ 'Mt',256.00], [ 'Fm',257.00], [ 'Md',258.00], [ 'Lr',260.00],  \
            [ 'Rf',261.00], [ 'Bh',262.00], [ 'Db',262.00], [ 'Sg',263.00], ['Uun',269.00], ['Uuu',272.00], ['Uub',277.00]
          ]
    temp=[]
    for a in compound[0]:
        for b in lib:
            if a == b[0]:
                temp.append(b[1])

    compound.append(temp)
#     return(compound)

#     for a in compound[0]:
#         for b in lib:
#             if a == b[1]:
#                 temp.append(b[0])

#     compound.append(temp)
#     return(compound)


## ------------------------------------------------------------------------------


def mlc_grycerol(initial_position_C=None, bond_leng=None):
    """
    1-2-3. mlc_grycerol: add 95th configuration of grycerol
    Generated by UC Kim at 2017 06 22
    """
    compound = [['C','O', 'H'], [3, 3, 8]]

    position= [['C', 1.318079978, -0.5108237259999999, -0.022517889999999596 , 'T', 'T', 'T' ],\
                 ['C', 2.2400230170000004, 0.5552262070000005, 0.5413368339999991 , 'T', 'T', 'T' ],\
                 ['C', 3.6814138290000002, 0.08655339500000103, 0.6524235009999995 , 'T', 'T', 'T' ],\
                 ['O', 0.0, 0.0, 0.0 , 'T', 'T', 'T' ],\
                 ['O', 4.230886698000001, -0.3631360829999988, -0.5680285400000002 , 'T', 'T', 'T' ],\
                 ['O', 2.2279642519999996, 1.7212823040000007, -0.2598315480000011 , 'T', 'T', 'T' ],\
                 ['H', -0.5653873090000001, -0.5642719559999989, -0.517293811 , 'T', 'T', 'T' ],\
                 ['H', 1.3221055269999997, 1.9483476880000001, -0.44187158400000115 , 'T', 'T', 'T' ],\
                 ['H', 4.234495759, 0.36649555000000156, -1.1786125600000013 , 'T', 'T', 'T' ],\
                 ['H', 1.6184166070000003, -0.7554851469999999, -1.0342277590000002 , 'T', 'T', 'T' ],\
                 ['H', 1.3887509700000003, -1.40960589, 0.583291053 , 'T', 'T', 'T' ],\
                 ['H', 1.8963165580000005, 0.8011206990000002, 1.5444263809999992 , 'T', 'T', 'T' ],\
                 ['H', 4.276929796000001, 0.8955478670000012, 1.0631129139999995 , 'T', 'T', 'T' ],\
                 ['H', 3.737431765, -0.7485926149999997, 1.3380751009999994 , 'T', 'T', 'T' ]]

    position = pst_poscar_move(position, - initial_position_C[0], - initial_position_C[1], - initial_position_C[2]) 
    return position


## ------------------------------------------------------------------------------


def mlc_O3(initial_position=None, bond_leng=None):
    """
    1-2-4. mlc_O3: add Ozone element
    Generated by UC Kim at 2018 07 17
    """
    compound = [['O'], [3]]      
    position= [['O', 0.0000000000000000,    0.0000000000000000,    0.0000000000000000, 'T', 'T', 'T' ],\
               ['O', 1.1027226066840949,   -0.0707783483217419,    0.6550693560555967, 'T', 'T', 'T' ],\
               ['O',-1.1019250332549309,   -1.1019250332549309,    0.6565392326959074, 'T', 'T', 'T' ]]

    position = pst_poscar_move(position, - initial_position[0], - initial_position[1], - initial_position[2]) 
    return position


## ------------------------------------------------------------------------------


def distance_atoms(a1, a2):
    if len(a1) == 3 and len(a2) ==3:
        pass
        return (  (a1[0]-a2[0])**2  +  (a1[1]-a2[1])**2 + (a1[2]-a2[2])**2 ) ** 0.5
    else:
        print_error(judge = True, message_from= 'distance_atoms')
        print( a1, a2)

## ------------------------------------------------------------------------------

def macroscopic_average(x,y,lp): 
    import numpy as np
    """
    y   : potential_function(x)
    lp  : length of phase to find the macroscopic average (target phase)
    nog : number of grid - to define 
    """
    m_aver = np.zeros(shape=(len(y)))            
    dx = max(x) / len(x)  # dx
    period_points = int( lp / dx)               # unit of integration
    # Period points must be even
    if period_points % 2 != 0 or period_points == 0: 
        period_points = period_points + 1
    len_y= len(y)
    for i in range(len_y):
        start = i - int(period_points / 2)
        end = i + int(period_points / 2)
        if start < 0:
            start = start + len_y
            m_aver[i] = m_aver[i] + sum(y[0:end]) + sum(y[start:len_y])
#            print("%10.9s %10.9s %10.9s %10.9s" %(x[start], x[end], str(period_points / 2), str(len(y))))
            m_aver[i] = m_aver[i] / period_points
        elif end >= len_y:
            end = end - len_y
            m_aver[i] = m_aver[i] + sum(y[start:len_y]) + sum(y[0:end])
            m_aver[i] = m_aver[i] / period_points
#            print("%10.9s %10.9s %10.9s %10.9s" %(x[start], x[end], str(period_points / 2), str(len(y))))
        else:
#            print("%10.9s %10.9s %10.9s %10.9s" %(x[start], x[end], str(period_points / 2), str(len(y))))
            m_aver[i] = m_aver[i] + sum(y[start:end]) / period_points

    macroscopic_average = np.average(m_aver)
#    print("Average of the average = ", numpy.average(m_aver))
    return macroscopic_average, m_aver





# Test for Lamms
# target = str(sys.argv[1])
# sequence_of_compound = [['Zr', 'O', 'Y'],['+4','-2','+3']]
# unitcell, compound, position = r_cryst_vasp(target)
# w_lammps_str(compound, position, "structure.data", unitcell, sequence_of_compound)


## Test for normal vector 
# temp, center=normal_vector([[6.0598862556776316,0.0000000000000000,4.3528408899082995],[ 6.0042276279323659,2.5098483562000000,4.2973650164883086],[4.2408565966276779,1.2419402764485341,5.5002446586139300]])
# shift_ads_w_normal_vec(center,1.09601,temp)



# ## Test for build_slab
# unitcell, compound, position=r_cryst_vasp('CONTCAR')
# num_atom=13 ## will be no_layer
# vaccume=18
# num_fix_lay=5
# output_name='cu331_13al.vasp'
# metal_kind=compound[0][0]
# lat_con=float(judge_pri_conv('CONTCAR', 'fcc'))
# ## build_slab(system, metal_kind, lat_con, index, output_name, vaccume,num_atom,num_fix_lay, symmetric)
# build_slab('fcc', compound[0][0], lat_con, '331', output_name, vaccume, num_atom, num_fix_lay, True)

# index == '100'
# index == '110'
# index == '111'
# index == '210'
# index == '211'
# index == '311'
# index == '331'

# ## TEST FOR mlc_ch4
# position,compound = mlc_ch4([0,0,0], 1.09601)
# print position
# w_poscar(position)


# ## TEST for pst_cell_expansion
# unitcell, compound, position=r_cryst_vasp('cu_111_4al.vasp')
# unitcell, compound, position=pst_cell_expansion(unitcell, position,[3,3,1])
# w_poscar(compound, position, 'test33.vasp', unitcell)

## TEST for ads
# unitcell, compound, position = r_cryst_vasp('cu_111_4al.vasp')
# slab_unitcell, compound, slab_position = pst_cell_expansion(unitcell, position,[3,3,1])
# # unitcell, compound, position = r_cryst_vasp('cu_311_7al.vasp')
# # slab_unitcell, compound, slab_position = pst_cell_expansion(unitcell, position,[2,3,1])
# # adsorbate_position, compound = mlc_ch4([6.1413235439999996,    1.8516787189999999, 10.3293388149376888], 1.09601)
# adsorbate_position, compound = mlc_ch3([1.2549236450848389/2. ,2.1735915129064742/2. , 8.2422636020000006], 1.09601)

# add_adsorbate(slab_position, slab_unitcell, adsorbate_position, 'cu111_b.vasp')

# unitcell, compound, position = r_cryst_vasp('test33.vasp')
# slab_unitcell, compound, slab_position = pst_cell_expansion(unitcell, position,[1,1,1])
# adsorbate_position, compound = mlc_ch4([0.0000000000000000,0.0000000000000000,8.2074228061374646], 1.09601)

# add_adsorbate(slab_position, slab_unitcell, adsorbate_position, 'cu111_b.vasp')


## ------------------------------------------------------------------------------
## Example for rotation_mlc
# for angle in range(10):
#     temp = 180./float(angle+1)
#     rotate_mlc('grycerol.vasp',  temp, 'z')



# ###########################################################################
# ###########################################################################
# ####   ____   ___    _   _  ___ _____   _____ ___  _   _  ____ _   _   ####
# ####  |  _ \ / _ \  | \ | |/ _ \_   _| |_   _/ _ \| | | |/ ___| | | |  ####
# ####  | | | | | | | |  \| | | | || |     | || | | | | | | |   | |_| |  ####
# ####  | |_| | |_| | | |\  | |_| || |     | || |_| | |_| | |___|  _  |  ####
# ####  |____/ \___/  |_| \_|\___/ |_|     |_| \___/ \___/ \____|_| |_|  ####
# ####                                                                   ####
# ###########################################################################
# ###########################################################################

