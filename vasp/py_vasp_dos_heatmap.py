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



def plot_heatmat_dos(folderpath, ax, points=200, len_dos = 1.5):
    import JH_lib as jh
    import pandas as pd
    
    tdos, atoms = read_doscar(folderpath)
    unitcell, compound, position = jh.r_cryst_vasp(folderpath+'/POSCAR')
    surf = 0
    for a in position: 
        if a[3] < 20 and a[3] > surf : surf = a[3]
    
    plotting  = False
    plotting2 = False; 
    plotting3 = True;  
    ndos = len(tdos)
    z_limits = [0,20]
    energy_limits=[-5,5]
    
    
    zlist = [ z_limits[1]/ points * i  for i in range(points)]
    df_plt= pd.DataFrame()
    
    noverlap = np.zeros((int(ndos)))
    for x in zlist:
        df_plt['%2.1f' %x]  = np.zeros((int(ndos)))
        
    df = pd.DataFrame(tdos, columns=['level','tdos'])
    
    
    ymax = ndos - int((df.level.max()-energy_limits[1]) / ((df.level.max()-df.level.min()) / ndos ))
    ymin = int((energy_limits[0]-df.level.min()) / ((df.level.max()-df.level.min()) / ndos ))

    print(ymax,ymin,surf)
    
    
    
    ele=1
    for a in atoms:
        df1 = pd.DataFrame(a, columns=['atom%i_s'%ele,'atom%i_p'%ele,'atom%i_d'%ele])
        df1['atom%i'%ele] = df1['atom%i_s'%ele] +df1['atom%i_p'%ele]+df1['atom%i_d'%ele]
        df =  pd.concat([df, df1], axis=1)
#         print(df1['atom%i'%ele].min(), df1['atom%i'%ele].max())

        if plotting:
            zlevel =[position[ele-1][3] for x in range(int(ndos))]
            clevel=df1['atom%i'%ele]
            plt.scatter(x=zlevel, y = df.level, s=0.1,
                        c=clevel, cmap=plt.cm.rainbow)
            plt.clim(10,30)

        if plotting2:
            for zz in zlist:
                if abs(position[ele-1][3] - zz) < 1.0:
                    zlevel = [zz for i in range(int(ndos))]
                    clevel = df1['atom%i'%ele]
                    plt.scatter(x=zlevel, y = df.level, s=0.1,
                                c=clevel, cmap=plt.cm.rainbow)
            plt.clim(10,15)
        if plotting3:
            for i in range(points):
                zz = zlist[i]
                if abs(position[ele-1][3] - zlist[i]) < len_dos:
                    weight_dist = 1 - abs(position[ele-1][3] - zlist[i])/len_dos
                    if df_plt['%2.1f' %zz].sum == 0:
#                         print('zero!')
                        df_plt['%2.1f' %zz] = df1['atom%i'%ele] * weight_dist
                    else:
                        noverlap[i] = noverlap[i] + 1
                        df_plt['%2.1f' %zz] = df1['atom%i'%ele] * weight_dist + df_plt['%2.1f' %zz]
        ele = ele+1
        
    if plotting3:
        import seaborn as sns
        df_plt.index= df.level
        if  '01_GaN_0001_' in folderpath:
            sns.heatmap(df_plt, vmin=0.10, vmax=1.5/4, ax= ax, cbar=False )
        else:
            sns.heatmap(df_plt, vmin=0.10, vmax=1.5, ax= ax, cbar=False )
        
        ax.set(ylim =[ymin, ymax],
               xlabel='Distance $(\\rm \\AA)$',
               ylabel='Energy Level (eV)',
               yticks=[-5,0,5],
#                xticks=[0,20],

              )
