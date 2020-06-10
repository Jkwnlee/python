%matplotlib inline
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def read_doscar(folderpath):
    doscar = open(folderpath+'DOSCAR', 'r')

    i = 1 ; atoms=[]
    for line in doscar:
        try:
            if i == 6: 
                [maxdos, mindos, ndos, efermi,temp] = line.split()
                column=['level','tdos']
                tdos=np.zeros((int(ndos),2))
            elif i > 6 and i < 6 + (1 + int(ndos)) * 1: 
                tdos[i-7] = [float(xx)-float(efermi) for xx in line.split()[:2]]
            elif i == 7 + int(ndos):
                #start to read atomic density of state
                natom=1
                atom_dos=np.zeros(((int(ndos),3)))
            elif i > 6 + (1 + int(ndos)) * natom and i < 6 + (1 + int(ndos)) * (natom+1):
                atom_dos[i- 7 - (1 + int(ndos)) * (natom)] = [float(xx) for xx in line.split()[1:]]

            elif i == 6 + (1 + int(ndos)) *  (natom+1) :
                if natom != 0:
                    atoms.append(atom_dos)
                    atom_dos=np.zeros(((int(ndos),3)))
                natom =  natom + 1
        except UnboundLocalError: pass
        i = i + 1
    return tdos, atoms
