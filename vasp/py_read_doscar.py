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

def read_doscar_long(folderpath):
    import JH_lib as jh
    doscar = open(folderpath+'/DOSCAR', 'r')
    start_atom_dos=False
    len_pdos, maxdos, mindos, ndos, efermi,temp , natom= 0,0,0,0,0,0,0
    i = 1 ; atoms=[]
    for line in doscar:
#         try:
        if i == 6: 
            [maxdos, mindos, ndos, efermi,temp] = line.split()
            column=['level','tdos']
            tdos=np.zeros((int(ndos),2))
        elif i > 6 and i < 6 + (1 + int(ndos)) * 1: 
            tdos[i-7] = [float(xx)-float(efermi) for xx in line.split()[:2]]
        elif i == 7 + int(ndos):
            natom=1
        elif i > 6 + (1 + int(ndos)) * natom and i < 6 + (1 + int(ndos)) * (natom+1):
            if start_atom_dos:
#                 print(line, i- 7 - (1 + int(ndos)) * (natom), atom_dos.shape)
                atom_dos[i- 7 - (1 + int(ndos)) * (natom)] = [float(xx) for xx in line.split()[1:]]
            else: 
                start_atom_dos = True
                len_pdos = len(line.split()[1:])
                atom_dos=np.zeros(((int(ndos),len_pdos))) 
#                 print(line, i- 7 - (1 + int(ndos)) * (natom), atom_dos.shape)
        elif i == 6 + (1 + int(ndos)) *  (natom+1) :
            if natom != 0:
                atoms.append(atom_dos)
                atom_dos=np.zeros(((int(ndos),len_pdos))) 
            natom =  natom + 1
        i = i + 1

    unitcell, compound, position = jh.r_cryst_vasp(folderpath+'/POSCAR')
    df_plt= pd.DataFrame()
    df = pd.DataFrame(tdos, columns=['level','tdos'])
    i=1 ;     num = 1;     element=''
    for a in atoms:
        if position[i-1][0] == element: ele = '%s%i' %(element,i)
        else: 
            num = 1
            ele = '%s%i' %(position[i-1][0],num)
            element = position[i-1][0]
            
        if len_pdos == 9: 
            colname=['%s_s'%ele,'%s_px'%ele,    '%s_py'%ele,'%s_pz'%ele,
             '%s_dxy'%ele,'%s_dyz'%ele,'%s_dz2r2'%ele,'%s_dxz'%ele,'%s_dx2y2'%ele]
            df1 = pd.DataFrame(a, columns=[colname])
            df1['%s_p'%ele] = df1[colname[1:4]].sum(axis=1)
            df1['%s_d'%ele] = df1[colname[4:]].sum(axis=1)
            df1['%s'%ele] = df1[colname].sum(axis=1)
        elif len_pdos == 3: 
            colname=['%s_s'%ele,'%s_p'%ele, '%s_d'%ele]
            df1 = pd.DataFrame(a, columns=[colname])
            df1['%s'%ele] = df1[colname].sum(axis=1)

        df =  pd.concat([df, df1], axis=1)
        i = i+1 ;         num = num + 1
    new_colname =['level', 'tdos']
    for x in df.columns[2:]:
        new_colname.append(x[0])
    df.columns = new_colname
    return df, position

import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

%matplotlib inline

folderpath = '/team/ptcad/leejh/workspace/2_TiN/c-TiN/slab/2_TSN/100/1_top/full_coverage/'
folderpath = '/team/ptcad/jhlee/a_workspace/01_Oxide/99_screening/00_calc_pbe/db_hynix_200605/HfO2_sp137/01_scf_HSE'
df, position = read_doscar_long(folderpath)
font_size, scale = 16, 8
cols , rows = 1, 1
plt.rcParams['font.size'] = font_size
fig, axes = plt.subplots(ncols =cols, nrows=rows, figsize=(cols*scale,rows*scale*1.2))
# for x in df.columns: print(x)
try:
    plt.plot(df.level, df['Hf1_d'], label = 'Hf1$d$')
    plt.plot(df.level, df['O1_p'], label = 'O1$p$')
    plt.xlim([-10, 10])
    plt.ylim([0, 3])
    plt.legend()
#     plt.hlines(y=0, xmin = -10, xmax = 10, colors = 'k', linestyles='--')
    plt.vlines(x=0, ymin = -10, ymax = 10, colors = 'k', linestyles='--')
except KeyError:
    print(df.columns)
