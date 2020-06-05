#!/bin/python
#Using code JH_lib
#Purpose: write a input file for vasp.dos
import JH_lib as jh
import os
import numpy as np
from optparse import OptionParser
#1. Read
import sys
############################################################
__version__ = "1.1"
############################################################

def command_line_arg():
    usage = "usage: %prog [options] arg1 arg2"  
    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-i", '--input', 
            action='store', type="string", dest='poscar',
            default='./POSCAR',
            help='location of the POSCAR')

    par.add_option("--outcar",
            action='store', type="string", dest='outcar',
            default='./OUTCAR',
            help='location of the OUTCAR (for fermi energy) ')

    par.add_option("-d", "--dist",
            action='store', type="float", dest='dist',
            default=2,
            help='Divide Structure with distance of given number to define layer (overlap is allowed by -o (0.5 A))')

    par.add_option("-s", "--shift",
            action='store', type="float", dest='shift',
            default=0.5,
            help='Shifting the DOS in Y direction (D:0.5)')

    par.add_option("-o", "--overlap",
            action='store', type="float", dest='overlap',
            default=1.0,
            help='Allowance for overlap between dos of nearest layer')

    par.add_option("--e1", '--element1',
            action='append', type="string", dest='ele1',
            default=[],
            help='specify which atoms to define place A')


    par.add_option("--e2", '--element2', 
            action='append', type="string", dest='ele2',
            default=[],
            help='specify which atoms to define interface')

    par.add_option("--sub",
            action='append', type="string", dest='subele',
            default=[],
            help='specify which atoms to define interface')

    par.add_option('--pdosoffset', 
            action='store', type="float", dest='pdosOffset',
            default=0.0,
            help='offset in pdos plot')

    par.add_option("-t", '--title', 
            action='store', type="string", dest='Dtitle',
            default=False,
            help='title of the graph')

    par.add_option('--fac', 
            action='store', type="float", dest='pdosFactor',
            default=2,
            help='scale factor of the pdos')

    par.add_option('-z', '--zero',
            action='store', type="float",
            dest='zero', default=0.0,
            help='energy reference of the band plot')

    par.add_option('--sigma',
            action='store', type="float",
            dest='sigma', default=0.02,
            help='smearing parameter, default 0.02')

    par.add_option('-n', '--nedos',
            action='store', type="int",
            dest='nedos', default=5000,
            help='number of point in DOS plot')

    par.add_option('--output',
            action='store', type="string", dest='outname',
            default='draw.sh',
            help='output shell script to draw dos "draw.sh" by default')
     
    par.add_option('--tofile',
            action='store', type="string", dest='dosname',
            default=False,
            help='make the dos result as text file')
    
    par.add_option('--skipwhite',
            action='store_false', dest='skip',
            default=True,
            help='neglect the lower vaccume region (draw only one)')

    par.add_option('--plot',
            action='store_true', dest='plotting',
            default=False,
            help='run the output file')
    
    par.add_option('--size', nargs=2,
            action='store', type="float", dest='figsize',
            default=(8, 9.0),
            help='figure size of the output plot')

    par.add_option('-x', nargs=2,
            action='store', type="float", dest='xlim',
            default=(-10, 5),
            help='x limit of the dos plot')


    return  par.parse_args( )

############################################################

def choose_cell(inputs):
  print "doping the given structure  ", inputs[1]
  filename    = inputs[1]
  target_ele  = inputs[2:]
  output_name = filename+'_targets.vasp'
  unitcell, compound, position = jh.r_cryst_vasp(filename)
  print  "Among the elements in the target file %s containing," %filename, compound[0], ", you choose:", target_ele
  new_position = []
  temp_position=[]
  for a in range(len(position)):
    for b in target_ele:
        if position[a][0] == b:
            temp_position = position[a]
            temp_position[0] = b
            print position[a]
            new_position.append(temp_position)
  if len(new_position[0]) == 7 and opts.Selec == False: 
      select =  False
  else:
      select =  opts.Selec
      
  new_compound, new_position1 = jh.component_from_positions(new_position)
  jh.w_poscar(position = new_position1,\
              compound = new_compound, \
              filename = output_name,  \
              unitcell = unitcell,     \
              Selective = select)
          

############################################################
############################################################
def temp(opts):
    i = 1  ;  spli=[]; fermi=[]
    unitcell, compound, position = jh.r_cryst_vasp(opts.poscar)
    width, height = opts.figsize
    xmin, xmax = opts.xlim
    out=open("%s" %opts.outname, "w")
    with os.popen("grep fermi %s |tail -1 | awk '{print $3}'" %opts.outcar) as a:
        temp = 1
        for line in a:
            if temp == 1: 
                fermi = line;  temp =+ 1
            else: pass
    out.write( "python /team/ptcad/jhlee/b_codework/py_vasp_post_process/py_vasp.dos.py \\\n" )
    out.write( " -y 0 %10.3f -z %8.5s --notot -q -s %i %i --fill -x %i %i \\"\
        %( (unitcell[2][2] * opts.shift +2)/opts.dist ,\
           fermi, \
           width, height,\
#           int(unitcell[2][2] / 2 / opts.dist),\
           xmin, xmax\
           )) 
    out.write("\n -p 0 --pdosoffset=%f --lc=white --fac=%.2f --pdosoffset=%f  -l '' \\" %(opts.shift, opts.pdosFactor, opts.shift))
    color = 'white'
    for a in range(int(unitcell[2][2]/opts.dist)+1): 
        temp=[]; i = 1; speci=[];zzzzz=[]
        for x in position:
            if x[3] >       a * opts.dist - opts.overlap and \
               x[3] < (a + 1) * opts.dist + opts.overlap:
                temp.append(i)
                zzzzz.append(x[3])
                speci.append(x[0])
            i += 1
        speci.sort(reverse=True)
        if sum(zzzzz) == 0 : la = 0
        else: la = sum(zzzzz)/len(zzzzz)
        if len(speci) == 0: 
            color = 'white'
            if la < unitcell[2][2]/2 and opts.skip: 
                out.write("\n -p 0 --pdosoffset=%f --lc=white --fac=%.2f --pdosoffset=%f  -l '' \\" %(opts.shift, opts.pdosFactor, opts.shift))
            if la > unitcell[2][2]/2:
                out.write("\n -p 0 --pdosoffset=%f --lc=white --fac=%.2f --pdosoffset=%f  -l '' \\" %(opts.shift, opts.pdosFactor, opts.shift))
        else:
            color = 'white'
            out.write( "\n -p '",)
            for xx in temp: out.write(" %s" %xx,)
            out.write( " '",)
            for x in opts.subele:
                if x in speci: color = 'red'
            if color != 'red':
                for x in opts.ele2:
                    if x in speci:
                        color = 'blue'
            if color != 'red' and color != 'blue':
                for x in opts.ele1:
                    if x in speci:
                        color = 'black'
#            if opts.element[0][0] in speci: 
#                if opts.element[0][1] in speci:    color ='red'
#                else:                              color ='black'
#            elif opts.element[0][1] in speci:      color ='blue'
#            else:
#                for c in opts.subelement[0]:
#                    if c in speci: 
#                        color = 'red'
            out.write( " --lc=%s  --fill  --fac=%.2f --pdosoffset=%f  -l " %(color,opts.pdosFactor, opts.shift),    )
            for yy in set(speci): out.write( "%s" %yy)
            out.write( "_%.1f \\" %la )
        if len(speci) == 0:
            print 'Dividing slab with %.2f, In the %2i (th) layer, there is no element'  %(opts.dist, a)
        else:
            print 'Dividing slab with %.2f, In the %2i (th) layer, there is element:'  %(opts.dist, a), set(speci), 'as color of %s' %color
    if opts.dosname:
        out.write( "\n  --tofile %s" %(opts.dosname))


if __name__ == '__main__':
    from time import time
    opts, args = command_line_arg()
    t0 = time()
    temp(opts) 
    t1 = time()
    print '\nCode calc completed! Time Used: %.2f [sec]\n' % (t1 - t0)
    if opts.plotting:  os.system('bash %s' %opts.outname)
  
