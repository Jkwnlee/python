#!/bin/python

### The python code for atomatic Wulff construction
### Written by Ji-Hwan Lee
### Last update:2018. 07. 06
### Used library:os, math


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


def write_CHGCHR_line(unitcell,chg_density,plot=False, output_name='output_chg.txt'):
    result=[]
    Volume = float(unitcell[0][0]) * float(unitcell[1][1]) * float(unitcell[2][2])
    f=open(output_name,'w' )
    f.write(' %10.9s %30.20s %30.20s \n' %('height(A)', 'chg_d (e/A^3)', 'chg_d (e)'))
    for a in range(len(local_potental)):
        f.write(' %10.9f  %30.20f %30.20f \n' \
                   %( float(unitcell[2][2])/len(chg_density[a]) * a + float(unitcell[2][2])/len(total)/2 ,  \
                     chg_density[a]/ Volume  ,  \
                     chg_density[a])
    f.close()

    if plot:
        import matplotlib.pyplot as plt

        f = open(output_name,'r')
        x=[]; y=[]; i = 1
        for a in f:
            if i == 1: i = i + 1
            else:
                x.append(a.split()[0])
                y.append(a.split()[1])
        plt.plot(x,y)
        plt.savefig(output_name.replace('txt', 'png'))



"""
This python code is written by Ji-Hwan Lee and Jongmin Yun for post-process of
Charge density (electron density) difference calculation in 1 Dimention.

1. You need to put following input in a same folder (in input POSCAR or not)
:   POSCAR, CHGCAR
* POSCAR: need to get volume and height
* CHGCAR: VASP output of CHGCAR of system

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
elif os.path.isfile('CHGCAR') == False:
    print('[ERROR]  There is no input files for structure: CHGCAR')
    print('[ERROR]  Stop the calculation')
    sys.exit()

if len(sys.argv) == 2:
    poscar = './input/POSCAR'
    chgcar = './input/CHGCAR'
    print '''\
[CODE] Automatically choose the input file: POSCAR and CONTCAR like,\n
[CODE] POSCAR : %s
[CODE] CHGCAR : %s
[CODE] 
[CODE] if you want to choose certain file,
[CODE] $ python  py_vasp_post_chgcar_to_1D.py  ./input/POSCAR ./input/CHGCAR ''' %(poscar, chgcar)
elif:
    poscar = sys.argv[1]
    chgcar = sys.argv[2]
    print '''\
[CODE] You choose the input file: POSCAR and CHGCAR like,\n
[CODE] POSCAR : %s
[CODE] CONTCAR: %s ''' %(poscar, chgcar)
                  


def main():
	unitcell, compound, position = r_cryst_vasp(poscar)
	total = read_CHGCHR_line(chgcar)
	write_CHGCHR_line(  unitcell        =unitcell  ,\
	                    chg_density     =total     ,\
	                    plot            =True      ,\
	                    output_name     ='total.txt')

if __name__ == "__main__":
	main()

		
