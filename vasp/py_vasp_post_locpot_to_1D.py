#!/bin/python

### The python code for read LOCPOT and convert into 1D along Z-direction
### Written by Ji-Hwan Lee
### Last update: 2018. 7. 5
### Used library:os, math, matplotlib in py_vasp_post_locpot_to_1D_fun_lib
### Detail functions are in 'py_vasp_post_locpot_to_1D_fun_lib.py'
### The python code for 
### Written by Ji-Hwan Lee


import os
import math 
import sys
import JH_lib as jh

def read_LOCPOT_line(filename):
	with open(filename, 'r') as f:
		density=[]; i = 1; n_point_total= 1; target=50; n_point_plane=0; line_density=[]
		for line in f:
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


def write_LOCPOT_line(unitcell, local_potental, plot=False,  output_name='output_lp'):
    f=open(output_name+'.txt','w' )
    f.write('%22.16s%22.16s \n' %('distance','local_potental'))
    for a in range(len(local_potental)):
        f.write('%22.16f%22.16f \n' %(float(unitcell[2][2])*a/len(local_potental),  local_potental[a] ))
    f.close()

    if plot:
        import matplotlib.pyplot as plt

        f = open(output_name+'.txt','r')
        x=[]; y=[]; i = 1
        for a in f:
            if i == 1: i = i + 1
            else:
                x.append(float(a.split()[0]))
                y.append(float(a.split()[1]))
        plt.plot(x,y)
        plt.savefig(output_name+'.png')
#        plt.savefig(output_name.replace('txt', 'png'))


def check_the_input_path(target):
    if os.path.isfile(target) == False:
        print('[ERROR]  There is no input files for structure: %s' %target)
        print('[ERROR]  Stop the calculation')
        sys.exit()

def main():
    if len(sys.argv) == 1:
        poscar = './input/POSCAR'
        locpot = './input/LOCPOT'
        ploting= False
        outname=sys.argv[len(sys.argv)-1]
        print '''\
        [CODE] Automatically choose the input file: POSCAR and LOCPOT:
        [CODE]     POSCAR : %s
        [CODE]     LOCPOT : %s
        [CODE] 
        [CODE] if you want to choose certain file,
        [CODE] $ python  py_vasp_post_chgcar_to_1D.py  ./input/POSCAR ./input/LOCPOT ''' %(poscar, locpot)
    elif len(sys.argv) == 3:
        poscar = sys.argv[1]
        locpot = sys.argv[2]
        outname=sys.argv[len(sys.argv)-1]
        ploting= False
        print '''\
        [CODE] You choose the input file: POSCAR and LOCPOT:
        [CODE]     POSCAR : %s
        [CODE]     LOCPOT : %s 
        [CODE] 
        [CODE] if you want to draw local potentail as function of z-axis,
        [CODE] $ python  py_vasp_post_chgcar_to_1D.py  ./input/POSCAR ./input/LOCPOT True''' %(poscar, locpot)
    elif len(sys.argv) == 4:
        poscar = sys.argv[1]
        locpot = sys.argv[2]
        ploting= sys.argv[3]
        outname=sys.argv[len(sys.argv)-1]
        print '''\
        [CODE] You choose the input file: POSCAR and LOCPOT:
        [CODE]     POSCAR : %s
        [CODE]     LOCPOT : %s 
        [CODE] 
        [CODE] You allow the drawing LOCPOT (taking a time a bit)''' %(poscar, locpot)
    elif len(sys.argv) < 2:
        outname='out_locpot_1D'

    check_the_input_path(poscar)
    check_the_input_path(locpot)


    print '        [CODE] Start the reading 3D LOCPOT and convert to 1D'
    unitcell, compound, position = jh.r_cryst_vasp(poscar)
    total = read_LOCPOT_line(locpot)
    print '        [CODE] Finish converting. Start to write'
    write_LOCPOT_line(  unitcell        =unitcell  ,\
	                local_potental  =total     ,\
	                plot            =ploting   ,\
	                output_name     =outname   )


if __name__ == "__main__":
    main()
