#bin/python

def r_cryst_vasp(filename):
    """
    Filename: VASP, POSCAR TYPE
    Function: Read POSCAR type 
    """
    poscar=[]; unitcell=[]; compound=[]; position=[];scales=[];num_atoms=0
    with open(filename, 'r') as f:
        Selective = False
        i = 1; j = 1; k =0
        for line in f:
            if line == None: break
            if i > 2 and i <6:
                unitcell.append(line.split())
            elif i > 5 and i < 8:
                compound.append(line.split())
            elif i == 8 :
                if j == 1 :
                    if line.split()[0][0] == 'S':
                        # print 'Selective Dynamics are applied'
                        Selective= True
                        i = i - 1
                    j = 2
                if j == 2:
                    if line.split()[0][0] == 'C':
                        scales=[1.,1.,1.]
                    elif line.split()[0][0] == 'D':
                        scales=[float(unitcell[0][0]) + float(unitcell[1][0]) + float(unitcell[2][0]),\
                                float(unitcell[0][1]) + float(unitcell[1][1]) + float(unitcell[2][1]),\
                                float(unitcell[0][2]) + float(unitcell[1][2]) + float(unitcell[2][2])]          
                if num_atoms == 0:
                    for temp in compound[1]: num_atoms=num_atoms+int(temp)

            elif i > 8 :
                if i <= 8 + num_atoms:
                    x = scales[0]  * float(line.split()[0])
                    y = scales[1]  * float(line.split()[1])
                    z = scales[2]  * float(line.split()[2])
                    if k < int(compound[1][0]):
                        if Selective:
                            position.append([compound[0][0], x, y, z, line.split()[3], line.split()[4], line.split()[5]])
                        else:
                            position.append([compound[0][0], x, y, z])
                    elif int(compound[1][0]) <= k <= int(compound[1][1]):
                        if Selective:
                            position.append([compound[0][1], x, y, z, line.split()[3], line.split()[4], line.split()[5]])
                        else:
                            position.append([compound[0][1], x, y, z])
                    k= k+1

            if i == 8 + num_atoms:
                return unitcell, compound, position
            else:
                i = i + 1
    # poscar=[unitcell, compound, position]

import os
if os.path.isfile('POSCAR') == False:
    print('[ERROR]  There is no input files for structure: POSCAR')
    print('[ERROR]  Stop the calculation')
    sys.exit()
else: 
    unitcell, compound, position = r_cryst_vasp('POSCAR')

compound = {'atom':compound[0], 'n_atom':compound[1],'mass':[]}

with open('POTCAR', 'r') as f:
    for line in f:
        if    'POMASS' in line: compound['mass'].append(float(line.split()[2].replace(';','')))
            #=   16.000; ZVAL   =    6.000    mass and valenz

mm=0
vf=0
for i in range(len(compound['atom'])):
    mm = mm + int(compound['n_atom'][i]) * compound['mass'][i]
    vf = vf + int(compound['n_atom'][i])

print compound
with open('./model.ing','w') as f:
    f.write("""set notrans
mm %s
vfree %s
""" %(mm, vf))
