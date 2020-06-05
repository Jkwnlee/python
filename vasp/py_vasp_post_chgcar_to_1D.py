#!/bin/python

### The python code for charge density difference (cdd) in 1D
### Written by Ji-Hwan Lee
### Last update:2018.7.5
### Used library:os, math
### INPUT : POSCAR, CHGCAR_A, CHGCAR_B,CHGCAR_Total
### OUTPUT: 'output.txt' including distance, cdd in unit of e/A^3, and  cdd in unit of e

import os
import math 
import sys
import time
import JH_lib as jh
from optparse import OptionParser
############################################################
__version__ = "1.0"
# 1.0
# Update from cdd
############################################################

def check_the_input_path(target):
    if os.path.isfile(target) == False:
        print('[ERROR]  There is no input files for structure: %s' %target)
        return False
    else: 
        print('[CODE] There is given File, %s .' %target) 
        return True

####################################################################
####################################################################

def read_CHGCHR_line(filename):
    density=[]; i = 1; line_density=[]; NGX = 0; NGY = 0 ; NGZ = 0;
    with open(filename, 'r') as f:
        xxx=1; target = -10000
        for line in f: 
            xxx = xxx + 1
            if  line == ' \n':        target = xxx
            elif xxx == target + 1:
                NGX, NGY, NGZ =  float(line.split()[0]), float(line.split()[1]),  float(line.split()[2])
                break
    npt = NGX * NGY * NGZ
    npp = NGX * NGY
    
    print "READING File, %20.20s , containing %10.i points  with %10.i points in a plane normal to set axis %s " %(filename, npt, npp, opts.direction)
    with open(filename, 'r') as f:
        for line in f:
            if 'augmentation' in line: break
            if i > target:     
                for a in range(len(line.split())):
                    density.append( float(line.split()[a]) )
                    if len(density) == npp: 
                        line_density.append(sum(density) / npp)
                        density=[]
            i = i + 1
    return line_density #, NGX, NGY, NGZ 



####################################################################
####################################################################

   
def command_line_arg():
    usage = """
Usage: %prog [options] arg1

This python code is written by Ji-Hwan Lee and Jongmin Yun for post-process of
Charge density (electron density) difference calculation in 1 Dimention.

1. You need to indicate following inputs with given options
:  --pab ./0_all/POSCAR  --pa  ./1_sub/POSCAR --pb ./2_on/POSCAR  -a ./1_sub/CHGCAR -b  ./2_on/CHGCAR  -c ./0_all/CHGCAR_AB 
* ./POSCAR: need to get volume and height
* ./CHGCAR: VASP output of CHGCAR of system
(You could choose LOCPOT instead of CHGCAR)

2. Working following    
* Read the electron density from each CHGCAR (LOCPOT) per cube 
(unit: e/A^3/N, where N is total number of cube)
* Average the electron density per xy layer 
* convert unit into  e/A^3 or e as function of z
"""

    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-p", '--poscar', 
            action='store', type="string", dest='poscar',
            default='./POSCAR',
            help='Path for the poscar, default: ./POSCAR')

    par.add_option("-c", '--chgcar', 
            action='store', type="string", dest='chg',
            default='./CHGCAR',
            help='Path for the CHGCAR, default: ./CHGCAR')

    par.add_option("-d", '--direction', 
            action='store', type="int", dest='direction',
            default=3,
            help='''
Direction for 1D plot: detualt-3 (z-direction)
1: x-direction
2: y-direction
3: z-direction''')
    
    par.add_option("--pv", '--plot_volume', 
            action='store', type="int", dest='plot_v',
            default=1,
            help='''
Ploting unit: detualt-1
1: (e/A^3)
2: (e) ''')

    par.add_option("--yl", '--ylabel', 
            action='store', type="string", dest='ylabel',
            default='cdd',
            help='Label of Y-axis, default: cdd')

    par.add_option("-o", '--output', 
            action='store', type="string", dest='outname',
            default='py_1D_CDD',
            help='OUTPUT TEXT, default: py_1D_CDD.txt')
             
    par.add_option("-g", '--grid', 
            action='store_false', dest='grid',
            default=True,
            help='Plot the grid in graph: Default-True')
  
    par.add_option("-f", '--figure', 
            action='store_true', dest='plotting',
            default=False,
            help='Plot the 1D graph or not: Default-False')

    par.add_option("-t", '--title', 
            action='store', type="string", dest='p_title',
            default='',
            help='title of the graph')

    par.add_option('-s', '--size', nargs=2,
            action='store', type="float", dest='figsize',
            default=(4.8, 3.0),
            help='figure size of the output plot')

    par.add_option('-x', nargs=2,
            action='store', type="float", dest='xlim',
            default=(0, 0),
            help='distance limit of the cdd plot')

    par.add_option('-y', nargs=2,
            action='store', type="float", dest='ylim',
            default=(0, 0),
            help='charge range of the cdd plot')
    
    par.add_option('--area', 
            action='store_true', dest='cal_area',
            default=False,
            help='Calculate area of charge')

    par.add_option('--label','-l', 
            action='store', type="string", dest='label',
            default='chgcar',
            help='label for the fitting system')


    par.add_option('--dpi', 
            action='store', type="int", dest='dpi',
            default=360,
            help='resolution of the output image')

    par.add_option('--macro', '--macroscopic', 
            action='store_true', dest='macroscopic',
            default=False,
            help='Calculate macroscopic average (integration 2 times)')
 
    par.add_option('--macro_dist', nargs=2,
            action='store', type="float", dest='macro_len',
            default=[4.0, 4.0],
            help='''
Define the length for  macroscopic average calculation (integration 2 times with l1, l2) 
where the default l1 = 4.0, and l2 = 4.0 Ang. ''')
    

    return  par.parse_args( )
    
####################################################################
####################################################################

def main(opts):
    if check_the_input_path(opts.outname+'.txt'):
        print("[CODE] There is output file ( %s )"  %(opts.outname+'.txt'))
        print("[CODE] Stop to re-reading/writing the potential")
        print("[CODE] If you want to run, please remove the file or change the outputname (-o)")
    elif    check_the_input_path(opts.poscar) and check_the_input_path(opts.chg):
        chg  =read_CHGCHR_line(opts.chg)
        length = len(chg)
        unitcell, compound, position = jh.r_cryst_vasp(opts.poscar)
        Volume = float(unitcell[0][0]) * float(unitcell[1][1]) * float(unitcell[2][2])
        result=[]
        f=open(opts.outname+'.txt','w' )
        f.write(' %10.9s %30.20s %30.20s\n'                        \
            %('height(A)', 'chg (e/A^3)', 'chg (e)'))
         
        for a in range(length):
            f.write(' %10.9s %30.20s %30.20s \n'                    \
                    %( float(unitcell[2][2])/length * a + float(unitcell[2][2])/length/2 ,  \
                        chg[a] / Volume  ,                       \
                        chg[a]))
    else: 
        print('[ERROR]  There is no input files')
        print('[ERROR]  Stop the calculation')
        sys.exit()
    
####################################################################
    if opts.plotting and check_the_input_path(opts.outname+'.txt'):
        import matplotlib
        import matplotlib.pyplot as plt

        #set the figure format 
        width, height = opts.figsize
        xmin, xmax = opts.xlim
        dpi = opts.dpi
        fig,(ax, ax0) = plt.subplots(2,1)
        #if opts.macroscopic:    fig,(ax, ax0) = plt.subplots(2,1)
        #else:                   fig,(ax, ax0) = plt.subplots(2,1)
        fig.set_size_inches(width, height)
        fig.tight_layout(pad=0.50)

        def set_plot(axx, opts): 
            if opts.plot_v == 2: 
                axx.set(xlabel='height(A)', ylabel=r'%s $(e)$' %opts.ylabel     )
            if opts.plot_v == 1:
                axx.set(xlabel='height(A)', ylabel=r'%s $(e/A^3)$' %opts.ylabel )        
            if opts.grid: axx.grid()
            axx.legend(fontsize='small',
                      frameon=True,
                      framealpha=0.6)
            if opts.p_title: 
                str2 = os.getcwd().split('/')  # get the title from the path
                axx.set_title('%s of %s' % (str2[len(str2)-1], opts.p_title))

            if xmax == 0:    pass
            else: axx.set_xlim([xmin,xmax])

        if opts.plot_v == 2: 
            V = 1.
        elif opts.plot_v == 1:
            unitcell, compound, position = jh.r_cryst_vasp(opts.poscar)
            V = float(unitcell[0][0]) * float(unitcell[1][1]) * float(unitcell[2][2])
 
        with open(opts.outname + '.txt' , 'r') as f:
            x=[]; y=[]; yp =[]; ym=[]; z=[]; i = 1 # ; NGZ = 0
            y = []
            y_a = []
            y_b = []
            for line in f:
                if i == 1: 
                    pass
                else:
                    x.append(float(line.split()[0]))   
                    y.append(float(line.split()[2])/V)
                i = i + 1

           
        ax.plot(x, y,     '-',  label=opts.label, color='black',  linewidth=1)

####################################################################
    if opts.macroscopic:
        l1= opts.macro_len[0] #4.2#abs(x_dipole[0] - x_dipole[1])   #length of 1st target phase
        l2= opts.macro_len[1]# 4.0#unitcell[2][2] - l1              #length of 2nd target phase#
        print 'plotting macroscopic average (integrated twice with length scale of %s %s' %(l1, l2)
        v_m_aver1, m_aver1 = jh.macroscopic_average(x,y,l1)
        v_m_aver2, m_aver2 = jh.macroscopic_average(x,m_aver1,l2)
        v_m_aver1, m_aver1 = jh.macroscopic_average(x,y,l1)
        v_m_aver2, m_aver2 = jh.macroscopic_average(x,m_aver1,l2)
       

        if opts.plotting:
            ax0.plot(x, m_aver2, '-', label='2nd integration', color= 'red', linewidth=1)


####################################################################
    if opts.plotting:
        set_plot(ax,opts)#,opts.outname)
        set_plot(ax0, opts)#,opts.outname+'0')
        fig.tight_layout(pad=0.50)
        fig.subplots_adjust(hspace= 0.30,\
                            left  = 0.10,\
                            bottom= 0.10,\
                            right = 0.90,\
                            top   = 0.90 \
                            )
        plt.savefig(opts.outname+'.png', dpi=opts.dpi)
        plt.show()

####################################################################
####################################################################

if __name__ == "__main__":

    opts, args = command_line_arg()
    main(opts)




    
