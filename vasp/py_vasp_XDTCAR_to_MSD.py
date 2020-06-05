import os
import math 
import sys
import time

from optparse import OptionParser
############################################################
__version__ = "1.1"
# 1.1
# ADD DATAFRAME
# Image number / counting time
# 1.0
# Merging two codes into one
# Add function to plot data
############################################################
####################################################################
####################################################################

   
def command_line_arg():
    usage = """
Usage: %prog [options] arg1

This python code is written by Ji-Hwan Lee and Jongmin Yun for post-process of
Charge density (electron density) difference calculation in 1 Dimention.

1. This code is merged two script from 'http://markummitchell.github.io/engauge-digitizer/'
"""
    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-x", '--xdatcar', 
            action='store', type="string", dest='xdatcar',
            default='XDATCAR',
            help='read XDATCAR with MD calulation, default: XDATCAR')

    par.add_option("--o1", '--outxyz', 
            action='store', type="string", dest='out_xyz',
            default='XDATCAR_fract.xyz',
            help='Write XDATCAR to fractional coordinates, default: XDATCAR_fract.xyz')

    par.add_option("--o2", '--outmsd', 
            action='store', type="string", dest='out_msd',
            default='XDATCAR_msd.xyz',
            help='Write XDATCAR to fractional coordinates, default: XDATCAR_msd.out')

    par.add_option("-p", '--plots', 
            action='store_true', dest='plotting',
            default=False,
            help='Plot the 1D graph or not: Default-False')

    return  par.parse_args( )
    
####################################################################
####################################################################



import numpy as np
import pandas as pd
from copy import deepcopy

def xdatcar_to_xyz(opts):
    xdatcar   = open(opts.xdatcar, 'r')
    xyz_fract = open(opts.out_xyz, 'w')
    lat_rec   = open(opts.xdatcar+'_lattice.vectors', 'w')


    system = xdatcar.readline()
    scale = float(xdatcar.readline().rstrip('\n'))

    #get lattice vectors
    a1 = np.array([ float(s)*scale for s in xdatcar.readline().rstrip('\n').split() ])
    a2 = np.array([ float(s)*scale for s in xdatcar.readline().rstrip('\n').split() ])
    a3 = np.array([ float(s)*scale for s in xdatcar.readline().rstrip('\n').split() ])

    #Save scaled lattice vectors
    lat_rec.write(str(a1[0])+' '+str(a1[1])+' '+str(a1[2])+'\n')
    lat_rec.write(str(a2[0])+' '+str(a2[1])+' '+str(a2[2])+'\n')
    lat_rec.write(str(a3[0])+' '+str(a3[1])+' '+str(a3[2]))
    lat_rec.close()

    #Read xdatcar
    element_names   = xdatcar.readline().rstrip('\n').split()
    element_dict    = {}
    element_numbers = xdatcar.readline().rstrip('\n').split()
    images_number   = -2
    i = 0
    N = 0
    for t in range(len(element_names)):
        element_dict[element_names[t]] = int(element_numbers[i])
        N += int(element_numbers[i])
        i += 1

    while True:
        images_number = images_number + 1
        line = xdatcar.readline()
        if len(line) == 0:
            break
        xyz_fract.write(str(N)+"\ncomment\n")
        for el in element_names:
            for i in range(element_dict[el]):
                p = xdatcar.readline().rstrip('\n').split() 
                coords = np.array([ float(s) for s in p ])
                cartesian_coords = coords[0]*a1+coords[1]*a2+coords[2]*a3 
                xyz_fract.write(el+ " " + str(coords[0])+ " " + str(coords[1]) + " " + str(coords[2]) +"\n")
    xdatcar.close()
    xyz_fract.close()
    
    return images_number

def read_lat_vec(opts):
    lat_file = open(opts.xdatcar+'_lattice.vectors','r')
    line = []
    for i in range(3):
        line.append([float(x) for x in lat_file.readline().rstrip('\n').split()])
        print(line[i])
    lattice = np.array([line[0],line[1],line[2]])
    return lattice
 # This function reads an XYZ file and a list of lattice vectors L = [x,y,z] and gives MSD + unwrapped coordinates

def MSD(opts,image_number):
    a = read_lat_vec(opts)
    l = [np.sqrt(np.dot(a[0],a[0])), np.sqrt(np.dot(a[1],a[1])), np.sqrt(np.dot(a[2],a[2]))] #basis vector lengths

    file      = open(opts.out_xyz, 'r')

    origin_list = [] # Stores the origin as [element,[coords]]
    prev_list = [] # Stores the wrapped previous step
    unwrapped_list = [] # Stores the instantenous unwrapped

    msd = [] #Stores atom-wise MSD  Stores msd as [msd]
    msd_dict ={} #Stores element-wise MSD
    msd_lattice = []
    msd_dict_lattice ={}

    element_list = [] # element list
    element_dict = {} # number of elements stored

    content = file.readline()
    N = int(content)

    for i in range(N):
        msd.append(np.float64('0.0'))
        msd_lattice.append([0.0, 0.0, 0.0 ])

    file.readline()
    step = 0

    while True:
        step += 1
        # Get and store the origin coordinates in origin_dict at first step
        if step == 1:
            for i in range(N):
                t = file.readline().rstrip('\n').split()
                element = t[0]
                if element not in element_list:
                    element_list.append(element)
                if element not in element_dict:
                    element_dict[element] = 1.0
                else:
                    element_dict[element] += 1.0
                coords = np.array( [ float(s) for s in t[1:] ] )
                origin_list.append([element,coords])
            # Copy the first set of coordinates as prev_dict and unwrapped
            unwrapped_list = deepcopy(origin_list)
            prev_list = deepcopy(origin_list)
            recorder_column=["step"]            
            for element in element_list:
                recorder_column.append(element+"_sum")  
                recorder_column.append(element+"_x")
                recorder_column.append(element+"_y")
                recorder_column.append(element+"_z")

            msd_df = pd.DataFrame(columns=recorder_column)
            
        # Read wrapped coordinates into wrapped_dict
        content = file.readline()
        if len(content) == 0:
            print ("\n---End of file---\n")
            break
        N = int(content)
        file.readline()
        wrapped_list = [] # Erease the previous set of coordinates
        for i in range(N):
            t = file.readline().rstrip('\n').split()
            element = t[0]
            coords = np.array( [ float(s) for s in t[1:] ] )
            wrapped_list.append([element,coords])

        # Unwrap coodinates and get MSD
        for atom in range(N):

            msd[atom] = 0.0

            # decompose wrapped atom coordinates to onto lattice vectors:
            w1, w2, w3 = wrapped_list[atom][1][0], wrapped_list[atom][1][1], wrapped_list[atom][1][2]

            # decompose prev atom coordinates to onto lattice vectors:
            p1, p2, p3 = prev_list[atom][1][0], prev_list[atom][1][1], prev_list[atom][1][2]

            #get distance between periodic images and use the smallest one
            if np.fabs(w1 - p1) > 0.5:  u1 = w1 - p1 - np.sign(w1 - p1)
            else:                       u1 = w1 - p1

            if np.fabs(w2 - p2) > 0.5:  u2 = w2 - p2 - np.sign(w2 - p2)
            else:                       u2 = w2 - p2

            if np.fabs(w3 - p3) > 0.5:  u3 = w3 - p3 - np.sign(w3 - p3)
            else:                       u3 = w3 - p3

            #add unwrapped displacements to unwrapped coords

            unwrapped_list[atom][1][0] += u1
            unwrapped_list[atom][1][1] += u2
            unwrapped_list[atom][1][2] += u3

            uw = unwrapped_list[atom][1][0]*a[0] + unwrapped_list[atom][1][1]*a[1] +unwrapped_list[atom][1][2]*a[2]
            ol =    origin_list[atom][1][0]*a[0] +    origin_list[atom][1][1]*a[1] +   origin_list[atom][1][2]*a[2]

            msd[atom] = np.linalg.norm(uw-ol)**2
            msd_lattice[atom] = [np.linalg.norm(uw[0]-ol[0])**2,np.linalg.norm(uw[1]-ol[1])**2,np.linalg.norm(uw[2]-ol[2])**2]


        prev_list = [] # Store current wrapped coordinates for the next step
        prev_list = deepcopy(wrapped_list)
        
        # record msd
        recorder_lines=[int(step)]

        for el in element_list:
             msd_dict[el] = 0.0
             msd_dict_lattice[el]=[0.,0.,0.]

        for atom in range(len(msd)):
            msd_dict[wrapped_list[atom][0]] += msd[atom]/element_dict[wrapped_list[atom][0]]
            for i in range(3):
                msd_dict_lattice[wrapped_list[atom][0]][i] += msd_lattice[atom][i]/element_dict[wrapped_list[atom][0]]

        for el in element_list:
            recorder_lines.append(float(msd_dict[el]))
            for temp in range(3): recorder_lines.append(float(msd_dict_lattice[el][temp]))
            
        record_series= pd.Series(recorder_lines, index = recorder_column)
        msd_df       = msd_df.append(record_series, ignore_index=True)
        
        if step % 10 == 0:
            step = float(step)
            sys.stdout.write("\r[" + "=" * int(step/ (image_number / 20)) \
                                   + " " * int((image_number - step)/ (image_number / 20)) + "]" \
                                   + "%9.3f% %" %float(step / (image_number / 100 )))
            sys.stdout.flush()
    file.close()
    sys.stdout.write("\n")
    return msd_df


        
####################################################################
####################################################################

if __name__ == "__main__":
    from time import time 
    t0 = time()
    opts, args = command_line_arg()
    
    image_number = xdatcar_to_xyz(opts)
    t1 = time()
    print ('\n[CODE] Converting INPUT-XDATCAR to xyz is completed! Time Used: %.2f [sec]\n' % (t1 - t0))
    msd =  MSD(opts,image_number)
    t2 = time()
    print ('\n[CODE] Calculating the Mean square displacement is completed! Time Used: %.2f [sec]\n' % (t2 - t1))

    print(image_number)


    if opts.plotting:
        set_plot(ax,opts)#,opts.outname)
        if opts.macroscopic:         set_plot(ax0, opts)  #,opts.outname+'0')
        fig.tight_layout(pad=0.50)
        fig.subplots_adjust(hspace= 0.30,\
                            left  = 0.20,\
                            bottom= 0.10,\
                            right = 0.90,\
                            top   = 0.90 \
                            )
        plt.savefig(opts.outname+'.png', dpi=opts.dpi)
        plt.show()
    
    print ('\n[CODE] calc completed! Time Used: %.2f [sec]\n' % (t2 - t0))
