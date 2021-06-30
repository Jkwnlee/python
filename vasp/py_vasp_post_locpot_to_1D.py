#!/bin/python

### The python code for read potential and convert into 1D along Z-direction
### Written by Ji-Hwan Lee
### Last update: 2019. 4. 17
### Used library:os, math, matplotlib in py_vasp_post_locpot_to_1D_fun_lib
### Detail functions are in 'py_vasp_post_locpot_to_1D_fun_lib.py'
### The python code for 
### Written by Ji-Hwan Lee


import os
import math 
import sys
import JH_lib as jh
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

############################################################
__version__ = "1.8"
# 1.8 Print option for python3
# 1.7
# WF without Graph Plot
# 1.4
# SKIP WHEN SAME OUTFILE IS THERE
# Update for target information
# Update for optparse (--input, --draw, --output)
############################################################

def read_potential_line(opts):
    density=[]; i = 1; line_density=[]; NGX = 0; NGY = 0 ; NGZ = 0;
    with open(opts.inputname, 'r') as f:
        xxx=1; target = -10000
        for line in f: 
            xxx = xxx + 1
            if  line == ' \n':        target = xxx
            elif xxx == target + 1:
                NGX, NGY, NGZ =  float(line.split()[0]), float(line.split()[1]),  float(line.split()[2])
                break
    npt = NGX * NGY * NGZ
    npp = NGX * NGY
    
    print("[CODE] READING File, %20.20s , containing %10.i points  with %10.i points in a plane normal to set axis %s "\
    %(opts.inputname, npt, npp, 'z'))#opts.direction)
    with open(opts.inputname, 'r') as f:
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


def read_fermi(opts):
    with open(opts.outcar, 'r') as f:
        i = 0
        while i == 0:
            for line in f:
                if 'E-fermi' in line:
                    fermi=line.split()[2]
                    #print line, fermi
                    i = 1
    return fermi

def old_read_potential_line(vasp_locpot):
  with open(vasp_locpot, 'r') as f:
    density=[]; i = 1; n_point_total= 1; target=50; n_point_plane=0; line_density=[]
    for line in f:
      if target  ==  i:
        n_point_total = float(line.split()[0]) * float(line.split()[1]) * float(line.split()[2])
        n_point_plane = float(line.split()[0]) * float(line.split()[1]) 
      if line == ' \n':
        target = i + 1
      if  i > target:
        for a in range(len(line.split())):
          density.append( float(line.split()[a]) )
          if len(density) == n_point_plane:
            line_density.append(sum(density) / n_point_plane)
            density=[]
      i = i + 1

  return line_density



####################################################################
####################################################################

def write_potential_line(unitcell, local_potental,  output_name= 'temp'):
    f=open(output_name+'.txt','w' )
    Ef = read_fermi(opts)
    f.write('%32.16s%32.16s %+10.9s\n' %(output_name+'_'+'d',output_name+'_'+'lp', Ef))
    for a in range(len(local_potental)):
        f.write('%32.16f%32.16f \n' %(float(unitcell[2][2])*a/len(local_potental),  local_potental[a] ))
    f.close()
#    if plot: plot_1D_plot(output_name+'.txt',ax, Ef)


def  write_lines(x,xh, ys, ysh, output_name):
    if len(ysh) == len(ys) and len(x) == len(ys[0]): 
        pass
    else: print"error"
    f=open(output_name+'.txt','w' )
    f.write('%32.16s' %(output_name+'_'+xh))
    for temp in ysh:   f.write('%32.16s' %(output_name+'_'+temp))
    f.write(' \n')
    for a in range(len(x)):
        f.write('%32.16f' %(x[a]))
        for temp in ys:   f.write('%32.16f' %(temp[a]))
        f.write(' \n')
 

def plot_1D_plot(target,ax):
    f = open(target,'r')
    x=[]; y=[]; i = 1; evac = 9999
    for a in f:
        if i == 1: 
            i = i + 1
            ll= a.split()[1]
            ef= float(a.split()[2])
        else:
            x.append(float(a.split()[0]))
            y.append(float(a.split()[1]))
            if float(a.split()[0])  > opts.yyy and evac == 9999: 
                evac = a.split()[1]
    ax.plot(x,y, label=opts.label)
    if opts.yyy:
        ax.plot([min(x), max(x)], [ef,ef], '--', label='$E_{\\rm fermi}$')
        print """
        filename       = %s 
        vaccume energy = %s
        fermi energy   = %s
        work-function  = %f
        """ %(target, evac, ef, float(evac)-float(ef))
    ax.legend()
    ax.set_ylabel('Potential (eV)')
    ax.set_xlabel('Distance ($\\rm \AA$)')
    ax.axis('on')
#    ax.tick_params(top='off', bottom='on', left='on', right='off', labelleft='on', labelbottom='on')

def macro(opts):
    f = open(opts.outname+'.txt','r')
    x=[]; y=[]; i = 1
    for a in f:
        if i == 1: i = i + 1
        else:
            x.append(float(a.split()[0]))
            y.append(float(a.split()[1]))
 
    l1= opts.macro_len[0] #4.2#abs(x_dipole[0] - x_dipole[1])   #length of 1st target phase
    l2= opts.macro_len[1]# 4.0#unitcell[2][2] - l1              #length of 2nd target phase#
    v_m_aver1, m_aver1 = jh.macroscopic_average(x,y,l1)
    v_m_aver2, m_aver2 = jh.macroscopic_average(x,m_aver1,l2)
   
    return x, m_aver2



####################################################################
####################################################################

def multi_plot_1D_plot(opts):
    import matplotlib as mpl
    import numpy as np
    fig = plt.figure()
    ax  = plt.subplot(111)
    clr = mpl.cm.get_cmap('jet', len(opts.multi))
    print('[CODE] Chosen files for multi-locpot plot: \n', opts.multi)
    print('[CODE] | %20.15s | %10.8s | %20.15s | %s ' %('Vaccum-Level', 'E-Fermi', 'WF (eV)', 'filename '))
    j = 0
    for t in opts.multi: 
        f = open(t,'r')
        evac = 9999
        x=[]; y=[]; i = 1
        for a in f:
            if i == 1: 
                i = i + 1
                ll= a.split()[1]
                ef= a.split()[2]
            else:
                x.append(float(a.split()[0]))
                y.append(float(a.split()[1]))
                if float(a.split()[0])  > opts.yyy and evac == 9999: 
                    evac = a.split()[1]
        ax.plot(x,y, label=t, c=clr(j))
        ax.plot([min(x), max(x)], [ef,ef], '--', label='$E_{\\rm f}$ (%s)'%t, c=clr(j))
        print '[CODE] | %20.15s | %10.8s | %20.15s | %s ' %(evac, ef, float(ef)-float(evac), t)
#        print t, evac, ef, float(ef)-float(evac)
        j += 1
    print('[CODE] Vaccum level is defined when the position is %s Angstrom (if False is printed, please using --yyy)' %opts.yyy)
    plt.savefig(opts.outname+'.png')
    plt.yticks() 
    plt.legend()
    plt.show()
#    plt.savefig(t.replace('txt', 'png'))


####################################################################
####################################################################
def relative_permitivity(dV,z,E):
    permi=[]
    permi_z = []
    for i in range(len(dV)):
        if i +1 < len(dV):
            temp = - (dV[i])
            tempz = z[i]
            permi.append(E/temp)
            permi_z.append(tempz)
    return permi, permi_z


####################################################################
####################################################################
def check_the_input_path(target):
    if os.path.isfile(target) == False:
        print('[ERROR]  There is no input files for structure: %s' %target)
       # print('[ERROR]  Stop the calculation')
        return False
     #   sys.exit()
    else: 
        print('[CODE]  %s' %target) 
        return True

####################################################################
####################################################################

def main(opts):
    print '''\
[CODE] Chosen input file:'''

    if check_the_input_path(opts.poscar) and check_the_input_path(opts.inputname):pass
    else: 
        print('[ERROR]  Stop the calculation')
        sys.exit
    if check_the_input_path(opts.outname+'.txt'):
        print('[CODE] there is same name for output, %s' %opts.outname)
        print('[CODE] Skip the reading input file: %s' %opts.inputname)
        #print('[ERROR]  Stop the calculation')
    else:
        print('[CODE] Start to read input file: %s' %opts.inputname)
        unitcell, compound, position = jh.r_cryst_vasp(opts.poscar)
        total = read_potential_line(opts)
        write_potential_line(  unitcell     =unitcell        ,\
                            local_potental  =total           ,\
        #                    plot            =opts.plotting   ,\
                            output_name     =opts.outname )  #  ,\
    #                        ax              =axx)
       #sys.exit()

    
    if opts.plotting: 
        print('''[CODE] Drawing: ON   with  output (-o, --outname) : %s .txt''' %opts.outname   )
        #set the figure format 
        width, height = opts.figsize
        xmin, xmax = opts.xlim
        ymin, ymax = opts.ylim
        if opts.macroscopic:
            if opts.permit:
                fig,(ax0, ax1, ax2) = plt.subplots(3,1)
                axx = ax0
            else:
                fig,(ax0, ax1) = plt.subplots(2,1)
                axx = ax0
            
        else:                   
            fig, ax0 = plt.subplots(1,1)
            axx = ax0
            plot_1D_plot(opts.outname+'.txt',ax= ax0)
    else: 
        axx=''
        print('''[CODE] Drawing: OFF   (If you want: use -f with/without -o or --outname)''')
        if opts.yyy:
            f = open(opts.outname+'.txt','r')
            x=[]; y=[]; i = 1; evac = 9999
            for a in f:
                if i == 1: 
                    i = i + 1
                    ll= a.split()[1]
                    ef= float(a.split()[2])
                else:
                    x.append(float(a.split()[0]))
                    y.append(float(a.split()[1]))
                    if float(a.split()[0])  > opts.yyy and evac == 9999: 
                        evac = a.split()[1]
            print("""
            filename       = %s 
            vaccume energy = %s
            fermi energy   = %s
            work-function  = %f
            """ %(opts.outname+'.txt', evac, ef, float(evac)-float(ef)))


   #Ef = read_fermi(opts)


    if opts.macroscopic:
        print('''[CODE] Calculating Macrocsropic average (differenticate)  : ON   ''' )
        target = opts.outname+'.txt'
        x, bar_V = macro(opts)
        d_bar_V = np.zeros(bar_V.shape,np.float)
        d_bar_V[0:-1] = np.diff(bar_V)/np.diff(x) 
        d_bar_V[-1] = (bar_V[-1] - bar_V[-2])/(x[-1] - x[-2]) 

        write_lines(     x = x                            ,\
                        xh = 'dist'                       ,\
                        ys = [bar_V, d_bar_V]      ,\
                       ysh = ['b', 'db'] ,\
               output_name =  opts.outname+'_m'            )
        if opts.permit:
            p, pz =  relative_permitivity(d_bar_V,x,opts.permit)

    if opts.plotting:
        if opts.macroscopic: 
            ax0.plot(x, bar_V, '-', label=opts.label+'_macro', color= 'red', linewidth=1)
            ax1.plot(x, d_bar_V, '-', label=opts.label+'_macro_differ', color= 'red', linewidth=1)
            if opts.permit:
                ax2.plot(pz, p,         '-', label=opts.label, color= 'red', linewidth=1)
                ax2.set_xlabel('relative permittivity $(1/\varepsilon)$')
        if xmin == xmax: pass
        else: ax0.set_xlim([xmin, xmax])
        if ymin == ymax: pass
        else: ax0.set_ylim([ymin, ymax])
        fig.set_size_inches(width, height)


#        ax0.yaxis.set_ticks(np.arange(ymin, ymax, 4))
#        ax0.axes.get_yaxis().set_visible(True)
#        plt.gca().set_aspect('equal', adjustable='box')

#        fig.tight_layout(pad=0.50)
#        ax0.legend(fontsize='small',
#                      frameon=True,
#                      framealpha=0.6)
#        fig.subplots_adjust(left  = 0.10,\
#                            bottom= 0.10,\
#                            right = 0.90,\
#                            top   = 0.90 \
#                            )
        plt.grid()
#        plt.savefig(opts.outname+'.png', dpi=opts.dpi)
        plt.show()
        
 
    print('[CODE] Finish converting and writting')




####################################################################
####################################################################

   
def command_line_arg():
    usage = "usage: %prog [options] arg1"  
    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-p", '--poscar', 
            action='store', type="string", dest='poscar',
            default='./POSCAR',
            help='Path for the poscar, default: ./POSCAR')

    par.add_option("-i", '--input', 
            action='store', type="string", dest='inputname',
            default='./LOCPOT',
            help='Path for the potential file (LOCPOT,CHGCAR), default: ./LOCPOT')

    par.add_option("-o", '--output', 
            action='store', type="string", dest='outname',
            default='py_1D_pot',
            help='OUTPUT TEXT, default: py_1D_pot.txt')
             
    par.add_option("-f", '--figure', 
            action='store_true', dest='plotting',
            default=False,
            help='Plot the 1D graph or not: Default-False')
           
    par.add_option("--permi", 
            action='store_true', dest='permit',
            default=False,
            help='Calculate permitivity: Default-False')

    par.add_option('--outcar', 
        action='store', type="string", dest='outcar',
        default='./OUTCAR',
        help='OUTCAR position where the Fermi energy value is in')

    par.add_option('--label','-l', 
            action='store', type="string", dest='label',
            default='1D',
            help='label for the fitting system')

    par.add_option('--dpi', 
            action='store', type="int", dest='dpi',
            default=360,
            help='resolution of the output image')

    par.add_option("-t", '--title', 
            action='store', type="string", dest='DOStitle',
            default=False,
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

    par.add_option("-M", '--multi',
            action='append', type="string", dest='multi',
            default=[],
            help='''
Multi-Plotting graphs using 1D txt file (output of this code for each system)
1. Generate output file for each cases
2. codes -M target1_path -M target2_path -M target3_path ...''')
    
    par.add_option("--yyy",
            action='store', type="float", dest='yyy',
            default=False,
            help='To define representative position for vaccum level, give the z value in middle of vacuum')

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



if __name__ == "__main__":

    opts, args = command_line_arg()
    if len(opts.multi) == 0:   main(opts)
    else: multi_plot_1D_plot(opts)

    
