#bin/python3

import pymatgen as pmg
import pandas   as pd
import seaborn  as sns
import numpy    as np
import matplotlib.pyplot as plt 
# need to check how to treat tensor of dielectric: 
# Hashin, Z., & Shtrikman, S. Physical Review, 130(1), 129–133. (1963). 
# Conductivity of Polycrystals.
# doi:10.1103/physrev.130.129 

df1 =  pd.read_json('Petousis2017_Scientific_Data.json')
# data from https://datadryad.org/stash/dataset/doi:10.5061/dryad.ph81h
# Petousis, I., Mrdjenovich, D., Ballouz, E. et al. Sci Data 4, 160134 (2017). 
# High-throughput screening of inorganic compounds for the discovery of novel dielectric and optical materials. 
# https://doi.org/10.1038/sdata.2016.134
# no data from Nature Communications volume 5, Article number: 4845 (2014) 

df2 =  pd.read_json('Jingyu2020_Scientific_Data.json')
# data from https://archive.materialscloud.org/2020.0010/v1
# Qu, J., Zagaceta, D., Zhang, W. et al. Sci Data 7, 81 (2020).
# High dielectric ternary oxides from crystal structure prediction and high-throughput screening. 
# https://doi.org/10.1038/s41597-020-0418-6

df2['band_gap'] = 0
df2['poly_total'] = 0

for i in range(df2.shape[0]):
    # User comment: 
    # (1) To plot result like Jingyu2020, getting only poly e (without adding electronics)
    # (2) The df2.e_poly[i][0] is same to (df2.e_total[i][0][0][0] + df2.e_total[i][0][1][1]+df2.e_total[i][0][2][2])/3
    # (3) Maybe.. it comes from summation of ionic and electronic term
   
    df2['poly_total'].iloc[i] = df2.e_poly[i][0] #\
                                #+ (df2.e_electronic[0][0][0][0] + df2.e_electronic[0][0][1][1]+df2.e_electronic[0][0][2][2])/3 
    df2['band_gap'].iloc[i] = df2.meta[i]['bandgap'][0][0]
    
fig, ax = plt.subplots(nrows=1, ncols=1)

plt1 = sns.scatterplot(data = df1 , x='band_gap', y='poly_total', color = 'black', size=1, 
                       ax = ax, label='Petousis2017')
plt2 = sns.scatterplot(data = df2 , x='band_gap', y='poly_total', color = 'red'  , size=1,
                       ax = ax, label='Jingyu2020')
ax.set(xscale = 'log', yscale = 'log', 
       xlim = [1e-1, 1e1] , ylim=[1e0,1e3], 
       ylabel = '$\\varepsilon_{\\rm poly}$ (arb.)', 
       xlabel = '$E_{\\rm g}$ (eV)')
ax.grid(which='major',axis='both', color='grey',linestyle='-', linewidth=0.5)
ax.grid(which='minor',axis='both', color='red',linestyle='--', linewidth=0.1)
# ax.legend(handles=[plt1, plt2])  ## Need to Update


## Text and line trend study 
x= np.array([1e-1, 1e1])
x_shift = 1.2
y_shift = 0.75
for c1 in [2,4,8,16]:
    ax.plot(x,(c1/x)**2,'g', '--', linewidth=0.5 )
    if   (c1/x[0])**2 < 1e3:
        ax.text(x[0]*x_shift, (c1/x[0])**2* y_shift, 'c=%i' %c1, color = 'g')
    else:
        ax.text(c1/(1e3)**0.5*x_shift, 1e3* y_shift, 'c=%i' %c1, color = 'g')
   

x_shift = 0.70
y_shift = 1.05
for c2 in [10,20,40,80, 160]:
    ax.plot(x,(c2/x),'r', '--', linewidth=0.5 )
    if   (c2/x[0]) < 1e3:
        ax.text(x[0]*x_shift, (c2/x[0])* y_shift, 'c=%i' %c2, color = 'r')
    else:
        ax.text(c2/1e3*x_shift, 1e3* y_shift, 'c=%i' %c2, color = 'r')


