#!/bin/python

#Using code JH_lib

#Purpose: choose a poscar
import JH_lib as jh
import os
import numpy as np
from optparse import OptionParser
#1. Read.
import sys
############################################################
__version__ = "1.0"
############################################################

def command_line_arg():
    usage = "usage: %prog [options] arg1 arg2"  
    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-s", '--Selective', 
            action='store_true', dest='Selec',
            default=False,
            help='Turn on the Selective Dynamics to T T T(Default: False)')


    return  par.parse_args( )

############################################################

def choose_cell(inputs):
  print "doping the given structure  ", inputs[1]
  filename    = inputs[1]
  target_ele  = inputs[2:]
#  num_ta_ele  = int(inputs[3])
#  change_ele  = inputs[4]
  output_name = filename+'_targets.vasp'
  unitcell, compound, position = jh.r_cryst_vasp(filename)
  print "Among the elements in the target file %s containing," %filename, compound[0], ", you choose:", target_ele
  #print epsilon, unitcell
  new_position = []
  temp_position=[]
  for a in range(len(position)):
    for b in target_ele:
      if position[a][0] == b:
        temp_position = position[a]
        temp_position[0] = b
        new_position.append(temp_position)
#        print temp_position[0]
#  print temp_position,num_ta_ele,target_ele
#  for a in new_position:
#  for a in new_position: print a
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
def temp(filename):
    i = 1 ; dist= 2 ; overlap=6.87639/6 * 1 ; y_shift = 0.5
    unitcell, compound, position = jh.r_cryst_vasp(filename)
    spli=[]; fermi=[]
    
    with os.popen("grep fermi OUTCAR |tail -1 | awk '{print $3}'") as a:
        temp = 1
        for line in a:
            if temp == 1: 
                fermi = line;  temp =+ 1
            else: pass
    print "python /team/ptcad/jhlee/b_codework/py_vasp_post_process/py_vasp.dos.py \\"
    print "-y 0 %10.3f -z %8.5s --notot -q -s 8 %i --fill -x -10 5  \\" \
           %((unitcell[2][2] * y_shift +2)/dist , fermi, int(unitcell[2][2] / 2 / dist))
    for a in range(int(unitcell[2][2]/dist)+1): 
        temp=[]; i = 1; speci=[]
        for x in position:
            if x[3] > a * dist - overlap and  x[3] <  (a + 1) * dist + overlap:
                temp.append(i)
                speci.append(x[0])
            i += 1
        speci.sort(reverse=True)
#        print "\
#python /team/ptcad/jhlee/b_codework/py_vasp_post_process/py_vasp.dos.py \
#-y 0 1.0 -z %8.5s --notot -q -s 8 2 --fill -x -10 5  --pdosoffset=%s\
#-p " %(fermi, str(a*dist)),


        print "\\\n  -p ",
        for xx in temp: print xx,
        #print " --lc=%s --fill  --fac=2 --pdosoffset=%f  -l " %(color, y_shift),    
        print " --fill  --fac=2 --pdosoffset=%f  -l " %(y_shift),    
        for yy in set(speci):
            sys.stdout.write(yy)
        if 'Hf' in speci: color ='black'
        if 'Ti' in speci: color ='blue' 
        print " --lc=%s "  %(color),
#        print " --lc=%s --fill  --fac=2 --pdosoffset=%f  -l " %(color, y_shift),    
#        print



if __name__ == '__main__':
  from time import time
  opts, args = command_line_arg()
  t0 = time()
  temp('POSCAR') 
  t1 = time()
  print '\nCode calc completed! Time Used: %.2f [sec]\n' % (t1 - t0)
# use bash/python dos_atoms_v1.py > draw_2A.sh ; bash draw_2A.sh
  
  
