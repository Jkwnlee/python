### The python code for atomatic Wulff construction
### Written by Ji-Hwan Lee
### Last update: 2018. 07. 24
### Used library:os, math, time, sys

### The python code for 
### Written by Ji-Hwan Lee


import os
import math 
import sys
import time


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


"""
## ------------------------------------------------------------------------------


def r_cryst_vasp(filename):
    """
    Filename: VASP, POSCAR TYPE
    Function: Read POSCAR type 
    """
    poscar=[]; unitcell=[]; compound=[]; position=[];scales=[];num_atoms=0
    with open(filename, 'r') as f:
        i = 1; j = 1; k =0; Selective=False; kn = 0
        for line in f:
            if line == None: break
            if i > 2 and i <6:
                unitcell.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
            elif i > 5 and i < 8:
                compound.append(line.split())
            elif i == 8 :
                if j == 1 :
                    if line.split()[0][0] == 'S':
                        # print 'Selective Dynamics are applied'
                        Selective= True
                        i = i - 1
                    j = 2
                if j == 2:
                    if line.split()[0][0] == 'C':
                        scales=[[1.,0,0],[0,1.,0],[0,0,1.]]
                    elif line.split()[0][0] == 'D':
                        scales=[[ float(unitcell[0][0]), float(unitcell[0][1]), float(unitcell[0][2])],\
                                [ float(unitcell[1][0]), float(unitcell[1][1]), float(unitcell[1][2])],\
                                [ float(unitcell[2][0]), float(unitcell[2][1]), float(unitcell[2][2])]]
                if num_atoms == 0:
                    for temp in compound[1]: num_atoms=num_atoms+int(temp)

            elif i > 8 :
                if i <= 9 + num_atoms:
                    x = scales[0][0] * float(line.split()[0]) + scales[1][0] * float(line.split()[1]) + scales[2][0] * float(line.split()[2])
                    y = scales[0][1] * float(line.split()[0]) + scales[1][1] * float(line.split()[1]) + scales[2][1] * float(line.split()[2])
                    z = scales[0][2] * float(line.split()[0]) + scales[1][2] * float(line.split()[1]) + scales[2][2] * float(line.split()[2])
                    if k <= int(compound[1][kn]):
                        if Selective:
                            position.append([compound[0][kn], x, y, z, line.split()[3], line.split()[4], line.split()[5]])
                        else:
                            position.append([compound[0][kn], x, y, z])
                    k= k+1

                    if k == int(compound[1][kn]):
                        kn = kn + 1
                        k = 0


            if i == 8 + num_atoms:
                return unitcell, compound, position
            else:
                i = i + 1


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


def component_from_positions(insert_position):
    new_compound=[]
    componets = []
    num_componets = []
    i = 0; k= 0
    insert_position = sorted(insert_position)
    for temp1 in range(len(insert_position)):
        if temp1 == 0:
            componets.append(insert_position[temp1][0])
            num_componets.append(1)
        else:
            if insert_position[temp1][0] == temp:
                if componets[i] == insert_position[temp1][0]:
                    num_componets[i] = num_componets[i] + 1
            else:
                i = i + 1
                componets.append(insert_position[temp1][0])
                num_componets.append(1)

        temp = insert_position[temp1][0]


    new_compound = [componets , num_componets]
    print "[CODE] The atomic postions are sorted by alphabet, Please revise your POTCAR"
    print new_compound
    return new_compound,insert_position
    

## ------------------------------------------------------------------------------


def judge_pri_conv(filename, system):
    ## Function check the condition of file wheather it is primitive or not
    ##
    ## Condition1: number atom == 1
    ## Condition2: nit
    unitcell, compound, position=r_cryst_vasp('CONTCAR')
    if len(compound[1]) == 1:
        # print "this cell include only one kind of atom"
        if compound[1][0] == '1':
            # print "this cell include only one of atom"
            True
    else:
        print "this cell include only more than one kind of atoms", compound[1]
        pass
    if system == 'fcc':
        material_name=compound[0][0]
        latcon=unitcell[0][0]
        if float(latcon) == 0.:
            # print "the unit cell is not conventional"
            latcon=float(unitcell[0][1]) * 2
        else:
            pass
    return latcon


## ------------------------------------------------------------------------------


def w_poscar(position, compound = None, filename = None, unitcell = None, Selective = False):
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
        compound,position = component_from_positions(position)

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
            print " WARNING the xyz is not proper 'unitcell' SHOULD CONTAIN ONLY 3 NUMBERS"
            break

    for temp in range(len(compound)):
        for temp2 in range(len(compound[temp])):
            f.write(' %4.3s' %compound[temp][temp2])
        f.write('\n')
        
    if Selective:
        f.write('Selective Dynamics\n')
        
    f.write('C\n')

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
            print " WARNING the iserted position is not proper. IT SHOULD CONTAIN ONLY 3 NUMBERS"
            print xyz
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
            print " WARNING the xyz is not proper IT SHOULD CONTAIN ONLY 3 NUMBERS"
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
            print " WARNING the xyz is not proper IT SHOULD CONTAIN ONLY 3 NUMBERS"
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


def print_error(True,a=False):
    print """
 _____ ____  ____   ___  ____
| ____|  _ \|  _ \ / _ \|  _ \\
|  _| | |_) | |_) | | | | |_) |
| |___|  _ <|  _ <| |_| |  _ <
|_____|_| \_\_| \_\\\___/|_| \_\\
"""
    if a != False:
        print "Error Message From the function of ", a


## ------------------------------------------------------------------------------


def normal_vector(array_for_three_point):
    """
    Calculate the normal vector with a inserted three points
    """

    if len(array_for_three_point) != 3:
        return print_error(True, 'normal_vector')
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
        return print_error(True, 'cross_product')
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
        return print_error(True, 'shift_ads_w_normal_vec')
        
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
        return print_error(True)
    
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
        print print_error(True, 'rotate_mlc')

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
    lib=[['1.0079','H'], ['4.0026','He'], ['6.941','Li'], ['9.0122','Be'],\
        ['10.811','B'], ['12.0107','C'], ['14.0067','N'], ['15.9994','O'],\
        ['18.9984','F'], ['20.1797','Ne'], ['22.9897','Na'], ['24.305','Mg'],\
        ['26.9815','Al'], ['28.0855','Si'], ['30.9738','P'], ['32.065','S'],\
        ['35.453','Cl'], ['39.0983','K'], ['39.948','Ar'], ['40.078','Ca'],\
        ['44.9559','Sc'], ['47.867','Ti'], ['50.9415','V'], ['51.9961','Cr'],\
        ['54.938','Mn'], ['55.845','Fe'], ['58.6934','Ni'], ['58.9332','Co'],\
        ['63.546','Cu'], ['65.39','Zn'], ['69.723','Ga'], ['72.64','Ge'],\
        ['74.9216','As'], ['78.96','Se'], ['79.904','Br'], ['83.8','Kr'],\
        ['85.4678','Rb'], ['87.62','Sr'], ['88.9059','Y'], ['91.224','Zr'],\
        ['92.9064','Nb'], ['95.94','Mo'], ['98','Tc'], ['101.07','Ru'],\
        ['102.9055','Rh'], ['106.42','Pd'], ['107.8682','Ag'], ['112.411','Cd'],\
        ['114.818','In'], ['118.71','Sn'], ['121.76','Sb'], ['126.9045','I'],\
        ['127.6','Te'], ['131.293','Xe'], ['132.9055','Cs'], ['137.327','Ba'],\
        ['138.9055','La'], ['140.116','Ce'], ['140.9077','Pr'], ['144.24','Nd'],\
        ['145','Pm'], ['150.36','Sm'], ['151.964','Eu'], ['157.25','Gd'],\
        ['158.9253','Tb'], ['162.5','Dy'], ['164.9303','Ho'], ['167.259','Er'],\
        ['168.9342','Tm'], ['173.04','Yb'], ['174.967','Lu'], ['178.49','Hf'],\
        ['180.9479','Ta'], ['183.84','W'], ['186.207','Re'], ['190.23','Os'],\
        ['192.217','Ir'], ['195.078','Pt'], ['196.9665','Au'], ['200.59','Hg'],\
        ['204.3833','Tl'], ['207.2','Pb'], ['208.9804','Bi'], ['209','Po'],\
        ['210','At'], ['222','Rn'], ['223','Fr'], ['226','Ra'], ['227','Ac'],\
        ['231.0359','Pa'], ['232.0381','Th'], ['237','Np'], ['238.0289','U'],\
        ['243','Am'], ['244','Pu'], ['247','Cm'], ['247','Bk'], ['251','Cf'],\
        ['252','Es'], ['257','Fm'], ['258','Md'], ['259','No'], ['261','Rf'],\
        ['262','Lr'], ['262','Db'], ['264','Bh'], ['266','Sg'], ['268','Mt'],\
        ['272','Rg'], ['277','Hs'], ['Darmstadtium','110'], ['Ununbium','112'],\
        ['Ununtrium','113'], ['Ununquadium','114'], ['Ununpentium','115'], ['Ununhexium','116'],\
        ['Ununseptium','117'], ['Ununoctium','118']]
    temp=[]
    for a in compound[0]:
        for b in lib:
            if a == b[1]:
                temp.append(b[0])

    compound.append(temp)
    return(compound)


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
        print_error(True, 'distance_atoms')
        print a1, a2



## ------------------------------------------------------------------------------


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
