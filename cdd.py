#!/bin/python

### The python code for atomatic Wulff construction
### Written by Ji-Hwan Lee
### Last update:20 jan 2017
### Used library:os, math, ase
### Purpose: make wulff construction  --> go to 
### step 1: Revise surface_slab_fcc_bcc into python
### setp 2: 

### The python code for 
### Written by Ji-Hwan Lee
### Please use this code and cite our surface energy library
### S.-H. Yoo, J.-H. Lee, Y.-K. Jung, and A. Soon, Phys. Rev. B 93, 035434 (2016).

import os
import math 
import sys
import time


def r_cryst_vasp(filename):
    """
    Filename: VASP, POSCAR TYPE
    Function: Read POSCAR type 
    """
    poscar=[]; unitcell=[]; compound=[]; position=[];scales=[];num_atoms=0
    with open(filename, 'r') as f:
        i = 1; j = 1; k =0
        for line in f:
            if line == None: break
            if i > 2 and i <6:
                unitcell.append(line.split())
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
                        scales=[1.,1.,1.]
                    elif line.split()[0][0] == 'D':
                        scales=[float(unitcell[0][0]) + float(unitcell[1][0]) + float(unitcell[2][0]),\
                                float(unitcell[0][1]) + float(unitcell[1][1]) + float(unitcell[2][1]),\
                                float(unitcell[0][2]) + float(unitcell[1][2]) + float(unitcell[2][2])]          
                if num_atoms == 0:
                    for temp in compound[1]: num_atoms=num_atoms+int(temp)

            elif i > 8 :
                if i <= 8 + num_atoms:
                    x = scales[0]  * float(line.split()[0])
                    y = scales[1]  * float(line.split()[1])
                    z = scales[2]  * float(line.split()[2])
                    if k < int(compound[1][0]):
                        if Selective:
                            position.append([compound[0][0], x, y, z, line.split()[3], line.split()[4], line.split()[5]])
                        else:
                            position.append([compound[0][0], x, y, z])
                    elif int(compound[1][0]) <= k <= int(compound[1][1]):
                        if Selective:
                            position.append([compound[0][1], x, y, z, line.split()[3], line.split()[4], line.split()[5]])
                        else:
                            position.append([compound[0][1], x, y, z])
                    k= k+1

            if i == 8 + num_atoms:
                return unitcell, compound, position
            else:
                i = i + 1
    # poscar=[unitcell, compound, position]


def read_CHGCHR_line(filename):
	with open(filename, 'r') as f:
		density=[]; i = 1; n_point_total= 1; target=50; n_point_plane=0; line_density=[]
		for line in f:
			if 'augmentation' in line:
				break
			if target  ==  i:
				# print i
				n_point_total = float(line.split()[0]) * float(line.split()[1]) * float(line.split()[2])
				n_point_plane = float(line.split()[0]) * float(line.split()[1]) 
				# print n_point_total
			if line == ' \n':
				target = i + 1
				# print i
			if  i > target:
				# print i
				for a in range(len(line.split())):
					# print float(line.split()[a])
					density.append( float(line.split()[a]) )
					if len(density) == n_point_plane:
						line_density.append(sum(density) / n_point_plane)
						# print sum(density) / n_point_plane
						density=[]
			i = i + 1

	return line_density

"""
This python code is written by Ji-Hwan Lee and Jongmin Yun for post-process of
Charge density (electron density) difference calculation in 1 Dimention.

1. You need to put following input in a same folder
:   POSCAR, CHGCAR_A, CHGCAR_B, CHGCAR_Total
* POSCAR: need to get volume and height
* CHGCAR_Total: VASP output of CHGCAR of total system (A+B)
* CHGCAR_A: VASP output of CHGCAR of a part of system (A)
* CHGCAR_B: VASP output of CHGCAR of a part of system (B)

2. Working following    
* Read the electron density from each CHGCAR per cube 
(unit: e/A^3/N, where N is total number of cube)
* Average the electron density per xy layer 
* convert unit into  e/A^3 or e as function of z

"""

if os.path.isfile('POSCAR') == False:
    print('[ERROR]  There is no input files for structure: POSCAR')
    print('[ERROR]  Stop the calculation')
    sys.exit()
elif os.path.isfile('CHGCAR_A') == False:
    print('[ERROR]  There is no input files for structure: CHGCAR_A')
    print('[ERROR]  Please change name of CHGCAR of a section of total system into: CHGCAR_A')
    print('[ERROR]  Stop the calculation')
    sys.exit()
elif os.path.isfile('CHGCAR_B') == False:
    print('[ERROR]  There is no input files for structure: CHGCAR_B')
    print('[ERROR]  Please change name of CHGCAR of a section of total system into: CHGCAR_B')
    print('[ERROR]  Stop the calculation')
    sys.exit()
elif os.path.isfile('CHGCAR_Total')== False:
    print('[ERROR]  There is no input files for structure: CHGCAR_Total')
    print('[ERROR]  Please change name of CHGCAR of a total system into: CHGCAR_Total')
    print('[ERROR]  Stop the calculation')
    sys.exit()


unitcell, compound, position = r_cryst_vasp('POSCAR')
Volume = float(unitcell[0][0]) * float(unitcell[1][1]) * float(unitcell[2][2])
sn=read_CHGCHR_line('CHGCAR_A')
cu=read_CHGCHR_line('CHGCAR_B')
total=read_CHGCHR_line('CHGCAR_Total')
result=[]
 

f=open('output.txt','w' )
f.write(' %10.9s %30.20s %30.20s \n' %('height(A)', 'cdd (e/A^3)', 'cdd (e)'))

for a in range(len(total)):
	f.write(' %10.9f  %30.20f %30.20f \n' \
        %( float(unitcell[2][2])/len(total) * a +float(unitcell[2][2])/len(total)/2 ,  \
           (total[a] - sn[a] - cu[a]) / Volume  ,  \
           ( total[a] - sn[a] - cu[a])  ))

