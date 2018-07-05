#!/bin/pyton2

# Ji-Hwan Lee

# graphical tool for

import numpy as np
import os
import matplotlib.pyplot as plt

############################################################
__version__ = "1.0"
############################################################
def aiAT_plot(systems):
    import matplotlib.pyplot as plt
    E_bulk_TiN = -19.628924       # eV/TiN
    E_bulk_Ti  = -15.668424/2     # eV/Ti atom
    E_mlc_N2   = -16.650181/2     # eV/N atom
    E_form = E_bulk_TiN - E_bulk_Ti -E_mlc_N2
    num_points = 10
    add_mu_range = 0.3
    ymax =0.3


# print E_formation_TiN
    color_phase=["#396C2B", "#000080", "#00B880", "#FFD727", "#49B021", "#FF4C3D", "#00FFFF", "#00C8FF", "#00C8FF", "#0096FF", "#0064FF", "#0000E6", "#0000E6", "#0000C8", "#000096", "#000096", "#000064", "#FF9696", "#FF6464", "#FF3232", "#FFD42F", "#FF8E00", "#FF2E2F", "#FF2540", "#D50000", "#860000", "#610000", "#00FF00", "#00D200", "#00AA00", "#008200", "#366E00", "#6E6E00", "#6E4F00", "#505050"]

    mu_min, mu_max = E_form - add_mu_range , add_mu_range 
		
# Resolution: Total number of points to be considered =  num_points^2
### Setting chemical potentials
    dmu_N = np.linspace(mu_min, mu_max, num_points)
    print '%20.18s\t%10.8s\t%10.8s' %('title', 'Rich', 'Lean')
    for system in systems:
        title = system[0];     N_Ti   = system[2];    N_N   = system[3]
        area  = system[6];     energy = system[7];
        G = (energy - E_bulk_TiN * N_Ti - ( N_N - N_Ti ) * ( dmu_N + E_mlc_N2 )) / area / 2.
        plt.plot(dmu_N, G, label = title) 

        print '%20.18s\t%10.8s\t%10.8s' \
                %(title, \
                  (energy - E_bulk_TiN * N_Ti - ( N_N - N_Ti ) * ( 0 + E_mlc_N2 ) )/ area / 2., \
                  (energy - E_bulk_TiN * N_Ti - ( N_N - N_Ti ) * ( E_form + E_mlc_N2 ) ) / area / 2.   )

    plt.legend()
    plt.plot(	[0,0],[-10, 10],'--',color='black')
    plt.plot(	[E_form, E_form],[-10, 10],'--' ,color='black')
    plt.xlim(mu_min, mu_max)
    plt.ylim(0,ymax)
    plt.ylabel('Surface energy (eV/$\\rm\AA^2$)')
    plt.xlabel('Chemical potential, $\Delta\mu_{\\rm N}$, (eV)')
    plt.text(0,ymax*1.01,'Rich')
    plt.text(E_form,ymax*1.01,'Lean')
    plt.savefig('aiat.png')

############################################################



systems = [    #   f .u,  Ti,   N,Zr, O,        Area,     energy
['TiN(100)',    'N9Ti9',  9.,  9., 0, 0,  9.028582, -175.35971],
['TiN(110)',  'N18Ti18', 18., 18., 0, 0, 25.536674, -345.09568],
['TiN(111):Ti', 'N4Ti5',  5.,  4., 0, 0, 7.8189819, -84.675409],
['TiN(111):N',  'N5Ti4',  4.,  5., 0, 0, 7.8189819, -85.44504],
['TiN(211)',  'N12Ti12', 12., 12., 0, 0, 22.115394, -228.05203],
['TiN(211)-N','N10Ti12', 12., 10., 0, 0, 22.115394, -207.23992]
]



if __name__ == '__main__':
    from time import time
#    opts, args = command_line_arg()
    t0 = time()
    aiAT_plot(systems)
    t1 = time()
    print '\nCode calc completed! Time Used: %.2f [sec]\n' % (t1 - t0)
   
    os.system('/usr/bin/gnome-open aiat.png &')
