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
__version__ = "2.1"
# 1.4
# Update for target information
# Update for optparse (--input, --draw, --output)
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
* 0_all/POSCAR: need to get volume and height
* 1_sub/POSCAR: need to find interface region
* 2_on/POSCAR : need to find interface region
* 0_all/CHGCAR: VASP output of CHGCAR of total system (A+B)
* 1_sub/CHGCAR: VASP output of CHGCAR of a part of system (A)
* 2_on/CHGCAR : VASP output of CHGCAR of a part of system (B)
(You could choose LOCPOT instead of CHGCAR)

2. Working following    
* Read the electron density from each CHGCAR (LOCPOT) per cube 
(unit: e/A^3/N, where N is total number of cube)
* Average the electron density per xy layer 
* convert unit into  e/A^3 or e as function of z
"""

    par = OptionParser(usage=usage, version= __version__)

    par.add_option("--pab", '--poscar_AB', 
            action='store', type="string", dest='poscar_AB',
            default='./0_all/POSCAR',
            help='Path for the poscar, default: ./0_all/POSCAR')

    par.add_option("--pa", '--poscar_A', 
            action='store', type="string", dest='poscar_A',
            default='./1_sub/POSCAR',
            help='Path for the poscar of part A, default: ./1_sub/POSCAR')

    par.add_option("--pb", '--poscar_B', 
            action='store', type="string", dest='poscar_B',
            default='./2_on/POSCAR',
            help='Path for the poscar of part B, default: ./2_on/POSCAR')

    par.add_option("-a", '--chgcar_A', 
            action='store', type="string", dest='chg_a',
            default='./1_sub/CHGCAR',
            help='Path for the CHGCAR of part A, default: ./1_sub/CHGCAR')

    par.add_option("-b", '--chgcar_B', 
            action='store', type="string", dest='chg_b',
            default='./2_on/CHGCAR',
            help='Path for the CHGCAR of part B, default: ./2_on/CHGCAR')

    par.add_option("--ab", '--chgcar_AB', 
            action='store', type="string", dest='chg_ab',
            default='./0_all/CHGCAR',
            help='Path for the CHGCAR of part AB, default: ./0_all/CHGCAR')

 
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
            default='cdd, $\Delta \\rho$',
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

    par.add_option('-y', nargs=2,
            action='store', type="float", dest='ylim',
            default=(0, 0),
            help='distance limit of the cdd plot')

    par.add_option('-x', nargs=2,
            action='store', type="float", dest='xlim',
            default=(0, 0),
            help='distance limit of the cdd plot')
   
    par.add_option('--area', 
            action='store_true', dest='cal_area',
            default=False,
            help='Calculate area of chage')
   
    par.add_option('--plot_all', 
            action='store_true', dest='plot_all',
            default=False,
            help='Plot all of potentials in picture')


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
    elif    check_the_input_path(opts.poscar_AB) and check_the_input_path(opts.chg_ab)  and \
            check_the_input_path(opts.poscar_A)  and check_the_input_path(opts.chg_a)   and \
            check_the_input_path(opts.poscar_B)  and check_the_input_path(opts.chg_b):
        chg_ab  =read_CHGCHR_line(opts.chg_ab)
        chg_a   =read_CHGCHR_line(opts.chg_a)
        chg_b   =read_CHGCHR_line(opts.chg_b)
        length = len(chg_ab)
        unitcell, compound, position = jh.r_cryst_vasp(opts.poscar_AB)
        Volume = float(unitcell[0][0]) * float(unitcell[1][1]) * float(unitcell[2][2])
        result=[]
        f=open(opts.outname+'.txt','w' )
        f.write(' %10.9s %30.20s %30.20s %30.20s %30.20s %30.20s \n'                        \
            %('height(A)', 'cdd (e/A^3)', 'cdd (e)', 'chg_ab (e)','chg_a (e)', 'chg_b (e)'))
         
        for a in range(length):
            f.write(' %10.9s %30.20s %30.20s %30.20s %30.20s %30.20s \n'                    \
                    %( float(unitcell[2][2])/length * a + float(unitcell[2][2])/length/2 ,  \
                        (chg_ab[a] - chg_a[a] - chg_b[a]) / Volume  ,                       \
                        (chg_ab[a] - chg_a[a] - chg_b[a]),                                  \
                         chg_ab[a], chg_a[a], chg_b[a]))
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
        ymin, ymax = opts.ylim
        dpi = opts.dpi
        #fig,(ax, ax0) = plt.subplots(2,1)
        if opts.macroscopic or opts.plot_all:    fig,(ax,ax0) = plt.subplots(2,1)
        else:                   fig,(ax) = plt.subplots(1,1)
        fig.set_size_inches(width, height)
        fig.tight_layout(pad=0.50)
        option1 = {'horizontalalignment': 'right', 'verticalalignment': 'baseline', 'fontsize':10}#, 'bbox':props}
        option2 = {'horizontalalignment':  'left', 'verticalalignment': 'baseline', 'fontsize':10}#, 'bbox':props}
 
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

            if xmax == 0:  pass
            else: axx.set_xlim([xmin,xmax])
            if ymax == 0:  pass
            else: axx.set_ylim([ymin,ymax])
            

        if opts.plot_v == 2: 
            V = 1.
        elif opts.plot_v == 1:
            unitcell, compound, position = jh.r_cryst_vasp(opts.poscar_AB)
            V = float(unitcell[0][0]) * float(unitcell[1][1]) * float(unitcell[2][2])
 
        with open(opts.outname + '.txt' , 'r') as f:
            x=[]; y=[]; yp =[]; ym=[]; z=[]; i = 1 # ; NGZ = 0
            y_ab = []
            y_a = []
            y_b = []
            for line in f:
                if i == 1: 
                    pass
                else:
                    x.append(float(line.split()[0]))   
                    y.append(float(line.split()[2])/V)
                    y_ab.append( float(line.split()[3])/V)
                    y_a.append(  float(line.split()[4])/V)
                    y_b.append(  float(line.split()[5])/V)
                    if float(line.split()[2]) > 0:   yp.append(float(line.split()[2])/V); ym.append(0)
                    if float(line.split()[2]) < 0:   yp.append(0); ym.append(float(line.split()[2])/V)
                i = i + 1

           
        ax.plot(x, y,     '-',  label='Difference', color='black',  linewidth=1)
        if opts.plot_all:
            ax0.plot(x, y_ab, '-',  label='(AB)'  , color='black',   linewidth=1)
            ax0.plot(x, y_a,  '--', label='(A)'   , color='silver',    linewidth=1)
            ax0.plot(x, y_b,  '--', label='(B)'   , color='grey',     linewidth=1)
        
        ax.plot([],[],linewidth=5, label='Electron Accumulated', color='r',alpha=0.5)           
        ax.plot([],[],linewidth=5, label='Electron Loss', color='g',alpha=0.5)
        ax.fill_between(x, yp, facecolor='r', alpha=0.5)
        ax.fill_between(x, ym, facecolor='g', alpha=0.5)   

####################################################################
    if opts.cal_area or opts.macroscopic:
        """
        Define the regions for area calculation (center of interface) 
        Condition for center of interface 
            1. Elements
            2. Charge density (or potential density difference)
        """
        x_dipole=[]
        centering=sum(x) / float(len(x))
        unitcell, compound, position_b = jh.r_cryst_vasp(opts.poscar_B);     
        spcl=[];  # 2nd phase coordinate list
        #Get the max/min coordination for the defined 2nd phase
        for a in position_b:
            if a[opts.direction] > 0 : 
                spcl.append(a[opts.direction])
            else: 
                spcl.append(a[opts.direction] + unitcell[opts.direction-1][opts.direction-1]) # Applying PBC 
        
        gap = max(x[len(x)-1] - x[len(x)-2], spcl[len(spcl)-1] - spcl[len(spcl) - 2]) * 2
        # When the sign of charge/potential difference change, 
        pcc = []  # potential/charge density sign chage coordinate
        for num in range(len(y)-1):
            if  y[num] * y[num+1] < 0.0 : # and abs(x[num]- float(a)) < gap : 
                pcc.append(x[num]) # potential/charge density sign chage coordinate

        if len(pcc) > 2: 
            temp1=[]; temp2=[]; temp3=[];less_cent=[]; higher_cent=[]
            for a in pcc: # x_dipole:
                if  a  < min(spcl) : 
                    temp1.append(a)
                if  a  > max(spcl) : 
                    temp2.append(a)
            x_dipole.append(max(temp1)) ; x_dipole.append(min(temp2))

        elif ( max(pcc)-centering )  *  (min(pcc) - centering) < 0 : 
            x_dipole = pcc
            
    if opts.cal_area:
        if opts.plotting:
            for a in set(x_dipole):
               ax.plot([a,a], [min(y), max(y)],'--', color='black', linewidth=1)

        ## Calculate Area
        import numpy as np
        from numpy import trapz
        dist_r11=[]; dist_r12=[]; dist_r21=[]; dist_r22=[] ;
        chg_q11 =[] ; chg_q12=[];  chg_q21=[];  chg_q22=[] ;
        D_qor11 =[] ; D_qor12=[];  D_qor21=[];  D_qor22=[] ;
        for num in range(len(y)): 
            if x[num] < min(x_dipole) or x[num] == min(x_dipole): 
                dist_r11.append(x[num]) ; chg_q11.append(y[num]); D_qor11.append(y[num]/x[num])
            elif x[num] > min(x_dipole) and x[num] < centering :  
                dist_r12.append(x[num]) ; chg_q12.append(y[num]) ; D_qor12.append(y[num]/x[num])
            elif x[num] < max(x_dipole) and x[num] > centering or x[num] == max(x_dipole) : 
                dist_r21.append(x[num]) ; chg_q21.append(y[num]) ; D_qor21.append(y[num]/x[num])
            else: 
                dist_r22.append(x[num]) ; chg_q22.append(y[num]) ; D_qor22.append(y[num]/x[num])

        print("area1-1 = %.12f; Dipole1-1 = %.12f" %(trapz(chg_q11, dist_r11), sum(D_qor11)))
        print("area1-2 = %.12f; Dipole1-2 = %.12f" %(trapz(chg_q12, dist_r12), sum(D_qor12)))
        print("area2-1 = %.12f; Dipole2-1 = %.12f" %(trapz(chg_q21, dist_r21), sum(D_qor21)))
        print("area2-2 = %.12f; Dipole2-2 = %.12f" %(trapz(chg_q22, dist_r22), sum(D_qor22)))

        if opts.plotting:
        #    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(x_dipole[0]*0.99, max(y)*0.95, r'%.5f (1-1)$\rightarrow$' % (trapz(chg_q11, dist_r11)), **option1 )
            ax.text(x_dipole[0]*1.01, max(y)*0.90, r'$\leftarrow$(1-2) %.5f'  % (trapz(chg_q12, dist_r12)), **option2) 
            ax.text(x_dipole[1]*0.99, max(y)*0.95, r'%.5f (2-1)$\rightarrow$' % (trapz(chg_q21, dist_r21)), **option1 ) 
            ax.text(x_dipole[1]*1.01, max(y)*0.90, r'$\leftarrow$(2-2) %.5f'  % (trapz(chg_q22, dist_r22)), **option2 )


####################################################################
    if opts.macroscopic:
        print opts.macro_len
        l1= opts.macro_len[0] #4.2#abs(x_dipole[0] - x_dipole[1])   #length of 1st target phase
        l2= opts.macro_len[1]# 4.0#unitcell[2][2] - l1              #length of 2nd target phase#
        v_m_aver1, m_aver1 = jh.macroscopic_average(x,y,l1)
        v_m_aver2, m_aver2 = jh.macroscopic_average(x,m_aver1,l2)

        v_m_aver1_ab, m_aver1_ab = jh.macroscopic_average(x,y_ab,l1)
        v_m_aver2_ab, m_aver2_ab = jh.macroscopic_average(x,m_aver1_ab,l2)
       
        sum_m1 = 0 ; i = 0
        i2=[]
        for a in range(len(x)): 
            if ( x[a] - min(x_dipole) )  * ( x[a] - max(x_dipole) )  < 0 and \
               abs(x[a] - centering) < (max(x_dipole) -min(x_dipole) ) / 6:
                sum_m1 = sum_m1 + m_aver2_ab[a]
                i = i  + 1
            if max(m_aver2_ab) == m_aver2_ab[a]: i2 = a

        sum_M1 = sum_m1/i
        if (m_aver2_ab[i2] - sum_M1) * (m_aver2_ab[0] - sum_M1) < 0 :
            sum_M2 = max(m_aver2_ab)
        else: 
            sum_M2 = min(m_aver2_ab)
        print "MACROSCOPIC AVERAGE: min - %s , max-of-center-regiron: %5.4f" \
                %(min(m_aver2_ab), sum_m1/i)  
        print "MACROSCOPIC AVERAGE: diff: %5.4f" %(sum_m1/i - min(m_aver2_ab)) 

        if opts.plotting:
            if opts.plot_all: ax0.plot(x, m_aver1_ab,'--', label='1st integration', color='blue', linewidth=1)
            ax0.plot(x, m_aver2_ab, '-', label='$\Delta \widebar \\rho$', color= 'red', linewidth=1)
            ax0.plot([min(x),centering], [sum_M1, sum_M1], ':', color='red' , linewidth=1)
            ax0.plot([min(x),centering], [sum_M2, sum_M2], ':', color='red' , linewidth=1)
            arwx0 = (min(x)+centering)/6;   arwy0 = min(sum_M2, sum_M1)
            darwx = 0                   ;   darwy = abs(sum_M1 - sum_M2 ) 
            ax0.arrow(arwx0, arwy0, darwx, darwy,  color='red', head_width=min(.5,darwy / 10), head_length=darwy / 6, length_includes_head=True, linewidth=1)
            ax0.text( arwx0, arwy0 + darwy/3 , r'%.3f'  % (darwy), **option2) 


####################################################################
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

####################################################################
####################################################################

if __name__ == "__main__":

    opts, args = command_line_arg()
    main(opts)




    
