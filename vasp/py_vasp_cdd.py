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
__version__ = "2.0"
# 1.4
# Update for target information
# Update for optparse (--input, --draw, --output)
############################################################

def check_the_input_path(target):
    if os.path.isfile(target) == False:
        print('[ERROR]  There is no input files for structure: %s' %target)
        return False
    else: 
        print('[CODE] Given File, %s, is there' %target) 
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
    print filename, npt, npp
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

1. You need to put following input in a same folder
:   (-p ) POSCAR, (-a ) CHGCAR_A, (-b ) CHGCAR_B, (-c )CHGCAR_AB
* POSCAR: need to get volume and height
* CHGCAR_AB: VASP output of CHGCAR of total system (A+B)
* CHGCAR_A: VASP output of CHGCAR of a part of system (A)
* CHGCAR_B: VASP output of CHGCAR of a part of system (B)

2. Working following    
* Read the electron density from each CHGCAR per cube 
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

    par.add_option("-c", '--chgcar_AB', 
            action='store', type="string", dest='chg_ab',
            default='./0_all/CHGCAR',
            help='Path for the CHGCAR of part AB, default: ./0_all/CHGCAR')

 
    par.add_option("-d", '--direction', 
            action='store', type="int", dest='direction',
            default=3,
            help='direction for 1D plot: detualt-3 (z-direction) \
1: x-direction\
2: y-direction\
3: z-direction')
    
    par.add_option("--pv", '--plot_volume', 
            action='store', type="int", dest='plot_v',
            default=1,
            help='Ploting unit: detualt-1 (z-direction) \
1: cdd (e/A^3)\
2: cdd (e)')
    
    par.add_option("-o", '--output', 
            action='store', type="string", dest='outname',
            default='py_1D_CDD',
            help='OUTPUT TEXT, default: py_1D_CDD.txt')
            
    par.add_option("-f", '--figure', 
            action='store_true', dest='plotting',
            default=False,
            help='Plot the 1D graph or not: Default-False')

    par.add_option("-t", '--title', 
            action='store', type="string", dest='p_title',
            default=[],
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
            help='Calculate area of chage')

    par.add_option('--dpi', 
            action='store', type="int", dest='dpi',
            default=360,
            help='resolution of the output image')

    par.add_option('--macroscopic', 
            action='store_true', dest='macroscopic',
            default=False,
            help='Calculate macroscopic average (integration 2 times)')
    

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
    
    if opts.plotting and check_the_input_path(opts.outname+'.txt'):
        import matplotlib
        import matplotlib.pyplot as plt
#set the figure format 
        width, height = opts.figsize
        xmin, xmax = opts.xlim
        dpi = opts.dpi
        if opts.macroscopic:    
            fig,(ax,ax0) = plt.subplots(2,1)
            #fig,(ax,ax0,ax1) = plt.subplots(3,1)
        #    ax  = plt.plot(1)
        #    ax0 = plt.plot(2)
        #    ax1 = plt.plot(3)
        else:  fig,(ax, ax0)=plt.subplot(1,2)
        fig.set_size_inches(width, height)
        fig.tight_layout(pad=0.50)

        def set_plot(axx):#, name): 
 #           axx.set_size_inches(width, height)
 #           axx.tight_layout(pad=0.50)
            if opts.plot_v == 2: 
                axx.set(xlabel='height(A)', ylabel=r'cdd $(e)$'    ,title='Charge Difference Density of %s' % (str2[len(str2)-1]))
            if opts.plot_v == 1:
                axx.set(xlabel='height(A)', ylabel=r'cdd $(e/A^3)$',title='Charge Difference Density of %s' % (str2[len(str2)-1]))        
            axx.grid()
            axx.legend(fontsize='small',
                      frameon=True,
                      framealpha=0.6)
            if opts.p_title:
               axx.set_title(opts.p_title)
            if xmax == 0: pass
            else: axx.set_xlim([xmin,xmax])
#            plt.savefig(name+'.png', dpi=opts.dpi)
#            plt.show()



        unitcell, compound, position = jh.r_cryst_vasp(opts.poscar_AB)
        Volume = float(unitcell[0][0]) * float(unitcell[1][1]) * float(unitcell[2][2])
         # Data for plotting
        str2 = os.getcwd().split('/')
        print str2
        if opts.plot_v == 2: 
            ax.set(xlabel='height(A)', ylabel=r'cdd $(e)$'    ,title='Charge Difference Density of %s' % (str2[len(str2)-1]))
            V = 1.
        if opts.plot_v == 1:
            ax.set(xlabel='height(A)', ylabel=r'cdd $(e/A^3)$',title='Charge Difference Density of %s' % (str2[len(str2)-1]))        
            V = Volume
 
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

           
        ax.plot(x, y,     '-',  label='charge density difference', color='black',  linewidth=1)
        ax0.plot(x, y_ab, '-', label='charge density (AB)'      , color='black',   linewidth=1)
        ax0.plot(x, y_a,  '--', label='charge density (A)'      , color='silver',    linewidth=1)
        ax0.plot(x, y_b,  '--', label='charge density (B)'      , color='grey',     linewidth=1)
        
        ax.plot([],[],linewidth=5, label='Eelectron Accumulated', color='r',alpha=0.5)           
        ax.plot([],[],linewidth=5, label='Eelectron Loss', color='g',alpha=0.5)
        ax.fill_between(x, yp, facecolor='r', alpha=0.5)
        ax.fill_between(x, ym, facecolor='g', alpha=0.5)   

    if opts.cal_area:
        x_dipole=[]
        centering=sum(x) / float(len(x))
        unitcell, compound, position_b = jh.r_cryst_vasp(opts.poscar_B);     temp=[]; 
        for a in position_b:
            if a[opts.direction] > 0 : temp.append(a[opts.direction])
            else: temp.append(a[opts.direction] + unitcell[opts.direction-1][opts.direction-1])
        xxx_dipole = [min(temp),max(temp)] ;temp1=[]
        for a in xxx_dipole:
            for num in range(len(y)-1):
                if y[num] * y[num+1] < 0.0 and abs(x[num]- float(a)) < 1. : 
                    for c in temp: 
                        if abs(x[num] - c) < 2: temp1.append(x[num])
        less_cent=[]; higher_cent=[]
        print set(temp1), min(temp), max(temp)
        for a in set(temp1):
            if   a < centering  and min(temp) > a: less_cent.append(a)
            elif a > centering  and max(temp) < a: higher_cent.append(a)
        print less_cent, higher_cent
        x_dipole.append(max(less_cent)) ; x_dipole.append(min(higher_cent))

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
            elif x[num] < max(x_dipole) and x[num] > centering or x[num] == max(x_dipole) : dist_r21.append(x[num]) ; chg_q21.append(y[num]) ; D_qor21.append(y[num]/x[num])
            else: #[num] > max(x_dipole) or x[num] == max(x_dipole):
                dist_r22.append(x[num]) ; chg_q22.append(y[num]) ; D_qor22.append(y[num]/x[num])

        print("area1-1 = %.3f; Dipole1-1 = %.3f" %(trapz(chg_q11, dist_r11), sum(D_qor11)))
        print("area1-2 = %.3f; Dipole1-2 = %.3f" %(trapz(chg_q12, dist_r12), sum(D_qor12)))
        print("area2-1 = %.3f; Dipole2-1 = %.3f" %(trapz(chg_q21, dist_r21), sum(D_qor21)))
        print("area2-2 = %.3f; Dipole2-2 = %.3f" %(trapz(chg_q22, dist_r22), sum(D_qor22)))

        if opts.plotting:
        #    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            option1 = {'horizontalalignment': 'right', 'verticalalignment': 'baseline', 'fontsize':10}#, 'bbox':props}
            option2 = {'horizontalalignment': 'left', 'verticalalignment': 'baseline', 'fontsize':10 }#, 'bbox':props}
            ax.text(x_dipole[0]*0.99, max(y)*0.95, r'%.3f (1-1)$\rightarrow$' % (trapz(chg_q11, dist_r11)), **option1 )
            ax.text(x_dipole[0]*1.01, max(y)*0.90, r'$\leftarrow$(1-2) %.3f'  % (trapz(chg_q12, dist_r12)), **option2) 
            ax.text(x_dipole[1]*0.99, max(y)*0.95, r'%.3f (2-1)$\rightarrow$' % (trapz(chg_q21, dist_r21)), **option1 ) 
            ax.text(x_dipole[1]*1.01, max(y)*0.90, r'$\leftarrow$(2-2) %.3f'  % (trapz(chg_q22, dist_r22)), **option2 )


    if opts.macroscopic:
        l1=4.2#abs(x_dipole[0] - x_dipole[1])   #length of 1st target phase
        l2=4.0#unitcell[2][2] - l1              #length of 2nd target phase#
        v_m_aver1, m_aver1 = jh.macroscopic_average(x,y,l1)
        v_m_aver2, m_aver2 = jh.macroscopic_average(x,m_aver1,l2)

        v_m_aver1_ab, m_aver1_ab = jh.macroscopic_average(x,y_ab,l1)
        v_m_aver2_ab, m_aver2_ab = jh.macroscopic_average(x,m_aver1_ab,l2)

        ax0.plot(x, m_aver1_ab,'--', label='1st integration', color='blue', linewidth=1)
        ax0.plot(x, m_aver2_ab,'-', label='2nd integration', color='red', linewidth=1)

        #ax1.plot(x, y,'--', label='cdd', color='black', linewidth=1)
        #ax1.plot(x, m_aver1,'--', label='1st integration', color='blue', linewidth=1)
        #ax1.plot(x, m_aver2,'-', label='2nd integration', color='red', linewidth=1)
        #set_plot(ax1)

    if opts.plotting:
        set_plot(ax )#,opts.outname)
        set_plot(ax0)#,opts.outname+'0')
        #set_plot(ax1)#,opts.outname+'1')
        fig.tight_layout(pad=0.50)
        fig.subplots_adjust(hspace= 0.30,\
                            left  = 0.10,\
                            bottom= 0.10,\
                            right = 0.90,\
                            top   = 0.90\
                            )

        plt.savefig(opts.outname+'.png', dpi=opts.dpi)
        plt.show()


####################################################################
####################################################################

if __name__ == "__main__":

    opts, args = command_line_arg()
    main(opts)




#python py_vasp_post_cdd.py -f --macro -s 10 20 -x 0 `head -5 0_all/CONTCAR |tail -1 | awk '{print $3}'` --pv 2 --title LOCPOT -a 1_sub/LOCPOT -b 2_on/LOCPOT -c 0_all/LOCPOT



