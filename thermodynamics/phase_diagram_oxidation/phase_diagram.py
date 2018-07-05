#!/bin/pyton2

# Ji-Hwan Lee

# graphical tool for ab initio plotting

import numpy as np
import os
import matplotlib.pyplot as plt

############################################################
__version__ = "2.0"
############################################################
def set_plot_2d(fig_a, label_swch=True,legend_swch =True): 
    if label_swch == True:
        # axis_font = {'fontname':'Times New Roman', 'size':'16'}
        fig_a.set_xlabel(r'$\Delta\mu_N$')#, **axis_font)
        fig_a.set_ylabel(r'$\Delta\mu_O$')#, **axis_font)
#        plt.setp(ltext_a, fontsize='small') 
    else:
        fig_a.set_xticklabels([])
        fig_a.set_yticklabels([])
    if legend_swch == True: #
        legend_a = fig_a.legend(bbox_to_anchor=(1.05,1),loc=2, ncol=1, borderaxespad=0., scatterpoints=1, markerscale=7, fancybox=True)
        frame_a = legend_a.get_frame().set_alpha(0.5)  #transparency of frame of legend
        ltext_a = legend_a.get_texts()

#    fig_a.annotate('annotate', xy=(2, 1), xytext=(3, 4), arrowprops=dict(facecolor='black', shrink=0.05))

############################################################
############################################################

def read_janaf(filename, num_points, max_temperature, n_atom):
    '''
     example: filename = 'input_JANAF.txt'
              num_points   = 500
              max_temperature = 1200
              Considering JanaF table for entropy contribution of ozone
    '''

    a =  1;		
    temperature = []; entropy = []; enthalpy =[]; Tentropy = []            ## -T \Delta S
    with open(filename, 'r') as f:
        for line in f:
            if a == 1:
                a = a + 1
            else:
                temp = line.split()
                if a == 2:
                    ground_enthalpy= float(temp[4])  
                    a = a + 1
		temperature.append(float(temp[0]))
		entropy.append(float(temp[2])  * 0.0103642688  / 1000. / n_atom)
		Tentropy.append(float(temp[2]) * float(temp[0]) * - 0.0103642688  / 1000./ n_atom ) ## J/mol to eV/atom
		enthalpy.append( ( float(temp[4])  - ground_enthalpy )  * 0.0103642688 / n_atom)   ## kJ/mol to eV/atom			
		# fit_function_dH = a * T + b * T**2 / 2- c/T + 2 * d * T ** (0.5)

			
    ## Fitting a equation of entropy as function of polynimial
    def S_shomate_eq(T,a,b,c,d,e,f,g):
        return a * ln(T / 1000.) + b * T / 1000. + c * (T / 1000.) ** 2 / 2 + d * (T / 1000.) ** 3 / 3 - e / ( T / 1000.) ** 2 / 2 + g

    def H_shomate_eq(T,a,b,c,d,e,f,h):
        return a *(T / 1000.) + b * (T / 1000.) ** 2 / 2 + c * (T / 1000.) ** 3 / 3 + d * (T / 1000.) ** 4 / 4 - e / ( T / 1000.)  + f - h

    def r2(points, fitted):
        bar  = np.sum(points)/len(points)
        ssres = np.sum( (points- fitted)**2. )
        ssreg = np.sum( (fitted- bar)**2. )
        sstot = np.sum( (points - bar)**2. )
        return 1 - ssres/sstot

    temp_range = np.linspace(min(temperature)+0.00001, max_temperature, num_points)		
    z1 = np.polyfit(temperature, entropy, 15)
    f1= np.poly1d(z1)
    new_entropy = f1(temp_range)
    r2_entropy = r2(entropy, f1(temperature))

    z2 = np.polyfit(temperature, enthalpy, 5)
    f2 = np.poly1d(z2)
    new_enthalpy = f2(temp_range)
    r2_enthalpy  = r2(enthalpy, f2(temperature))
    
    plt.plot(temperature,Tentropy, 'o',label = '-Temp x Entropy, $-T*S(T)$ (eV/atom)')
    plt.plot(temperature,enthalpy, 'o',label = 'Enthalpy, $H(T)$ (eV/atom)')
    plt.plot(temp_range, - new_entropy * temp_range, '--', label = 'Fitted $-T*S(T)$ (eV) R2=%4.4f' %r2_entropy)
    plt.plot(temp_range, new_enthalpy, '--',label = 'Fitted $H(T)$ (eV) R2=%4.4f' %r2_enthalpy)
    plt.legend()
    plt.xlim([min(temperature), max_temperature])

    plt.ylim([-3, 1])#[ymin,ymax])
    plt.xlabel('Temperature (K)')
    plt.ylabel('Energy (eV/atom)')

#    plt.axhline(0, color='black')
    plt.savefig('Fig1_Janaf.png')
    plt.close()

    return  new_enthalpy, new_entropy, temp_range

############################################################
############################################################

def aiAT_2D_plot(system, num_points):
    fig, ((fig_a))= plt.subplots(1, 1)
    fig.subplots_adjust(left=0.10, bottom = 0.15, right =0.6, top = 0.90, wspace=0.5, hspace=0.7)
   # fig.set_size_inches(3, 4) # (width, height) in inches

    matrix_G = []
    for a in system:
        Etot = a[1];	N_N = a[2];	N_O = a[3];	N_Ti = a[4];	Area = a[5]
        mu_O = dmu_O + E_mlc_O_ref 
        mu_N = dmu_N + E_mlc_N2
        G = (Etot - E_bulk_TiN * N_Ti - (N_N - N_Ti) * mu_N -  N_O * mu_O ) / Area / 2. 
        matrix_G.append(G)
    matrix_G = np.array(matrix_G)
    matrix_Gmin =  np.amin(matrix_G, axis=0)
    for num_phases in range(len(system)):
        xy_coor = []
        for i in range(num_points):
            for j in range(num_points):
                if matrix_G[num_phases][i][j] <= matrix_Gmin[i][j]:
                    x, y = dmu_N[i][j], dmu_O[i][j]
	            xy_coor.append([x, y])
        xy_coor   = np.array(xy_coor)


        if not len(xy_coor) == 0: 
            fig_a.scatter(xy_coor[:,0], xy_coor[:,1], s=4500/num_points, marker='.',\
                          label=system[num_phases][0], c=system[num_phases][6])#, c=color[num_phases], lw=0, )
        print 'Phase %2.i / %2.i : %s' %(num_phases+1,len(matrix_G), system[num_phases][0])
    xmin, xmax = mu_min1, mu_max1
    ymin, ymax = mu_min2, mu_max2
    set_plot_2d(fig_a, label_swch=True, legend_swch = True)
    fig_a.plot(	[0,0],[ymin, ymax], '--',color='black') ## N Rich condition
    fig_a.plot(	[E_form, E_form],[ymin, ymax],'--' ,color='black') ## N Lean Condition
    fig_a.set_xlim([xmin,xmax])
    fig_a.set_ylim([ymin,ymax])
    fig_a.plot(	[xmin, xmax],[0,0],'--',color='black') ## O Rich condition
    fig_a.plot(	[xmin, xmax],[E_form2,E_form2],'--' ,color='black') ## O Lean Condition

    fig_a.text(0.00   -0.11,ymax    + 0.10,'$\downarrow \mu_{\\rm N}$(Rich), $p_{\\rm N}$(high)')
    fig_a.text(E_form -0.11,ymax    + 0.10,'$\downarrow \mu_{\\rm N}$(Lean), $p_{\\rm N}$(low)')
    fig_a.text(xmax   +0.10,E_form2 + 0.00,'$\leftarrow \mu_{\\rm O}$(TiO$_{2}$), $p_{\\rm O}$(Low)')
    fig_a.text(xmax   +0.10,0       + 0.00,'$\leftarrow \mu_{\\rm O}$(Rich), $p_{\\rm O}$(High)')

    fig.savefig('Fig2_2D_aiat.png')
    plt.close()

############################################################
############################################################

def aiAT_1D_2D_plot(system, num_points, enthalpy, entropy, temp_range):
    
    fig,((fig_a1, fig_a2), (fig_b1, fig_b2))= plt.subplots(2, 2)

    fig.subplots_adjust(left=0.11, bottom = 0.10, right =0.70, top = 0.90, wspace=0.5, hspace=0.5)
    fig.set_size_inches(8,8) # (width, height) in inches

    p_min, p_max = 10**(-99), 10**(-1)
    dp_range = np.linspace(p_min, p_max, num_points)# np.linspace(p_min, p_max, num_points)
    dT, dp = np.meshgrid( temp_range, dp_range)
#    dp, m_enthalpy = np.meshgrid(dp_range,enthalpy)
#    dp, m_entropy  = np.meshgrid(dp_range, entropy)

    for mu_N in 0 , E_form:
        if mu_N == 0: 
            ffig1 = fig_a1; ffig2 = fig_a2
        elif mu_N == E_form: 
            ffig1 = fig_b1; ffig2 = fig_b2

        matrix_G = []
        for a in system:
            Etot = a[1];	N_N = a[2];	N_O = a[3];	N_Ti = a[4];	Area = a[5]; p0 = 1
            dmu_N = mu_N + E_mlc_N2 

            mu_O1 = dmu_O + E_mlc_O_ref# + ( enthalpy - T1 * entropy ) +  1/2. * kb*dT*np.log(p1/p0)            
            G1 = (Etot - E_bulk_TiN * N_Ti - (N_N - N_Ti) * dmu_N - N_O * mu_O1 ) / Area / 2. 
            ffig1.plot(dmu_O, G1, label=a[0], color=a[6]) 

            mu_O = dmu_O + E_mlc_O_ref + ( enthalpy - dT * entropy ) +  1/2. * kb * dT * np.log(dp/p0) 
            G  = (Etot - E_bulk_TiN * N_Ti - (N_N - N_Ti) * dmu_N - N_O * mu_O ) / Area / 2.  
            matrix_G.append(G)

        ffig1.set_ylim([-0.2,0.1])
        ffig1.set_xlim([mu_min2, mu_max2])

        ffig1.plot([0,0],[-100,100],'--',color='black') ## O3 condition
        ffig1.plot([E_mlc_O2-E_mlc_O3, E_mlc_O2-E_mlc_O3],[-100,100],'--',color='black') ## O2 condition
        ffig1.plot([E_form2,E_form2],[-100,100],'--' ,color='black') ## TiO2 Condition


        ffig1.set_ylabel('$\Delta G \sim G_{\\rm surf} - G_{\\rm bulk} $ (eV/${\\rm \AA}^2$)')
        ffig1.set_xlabel('$\Delta\mu_O$ (eV)')#, **axis_font)
        ffig1.text(mu_min2-0.5,.15,'[$\Delta\mu=$%2.3f (eV) condition]'%mu_N)
        ffig1.text(0.00              -0.22,0.105,'$\downarrow$O$_3$(g)')
        ffig1.text(E_mlc_O2-E_mlc_O3 -0.22,0.105,'O$_2$(g)\n$\downarrow$')
        ffig1.text(E_form2           -0.22,0.105,'TiO$_2$(bulk)\n$\downarrow$')

        
        matrix_G = np.array(matrix_G)
        matrix_Gmin =  np.amin(matrix_G, axis=0)

        for num_phases in range(len(system)):
            xy_coor = []
            for i in range(num_points):
                for j in range(num_points):
                    if matrix_G[num_phases][i][j] <= matrix_Gmin[i][j]:
#                        x, y =  dT[i][j], dp[i][j]   ##X-axis: temperature. Y-axis: pressure
                        x, y = dp[i][j], dT[i][j]    ##Y-axis: temperature. X-axis: pressure
                        xy_coor.append([x, y])
    	    xy_coor   = np.array(xy_coor)

            if not len(xy_coor) == 0:
                ffig2.scatter(xy_coor[:,0], xy_coor[:,1], s=4500/num_points, marker='.',\
                             label=system[num_phases][0], c=system[num_phases][6])
            print 'Phase %2.i / %2.i : %s' %(num_phases+1,len(matrix_G), system[num_phases][0])
        set_plot_2d(ffig2, label_swch=True, legend_swch = True)
        ffig2.set_ylim([300,800])
        ffig2.set_xlim([0,0.1])#[p_min,p_max])
        ffig2.set_xlabel('${\\rm O_3}$ pressure (atm)')#, **axis_font)
        ffig2.set_ylabel('$T\,({\\rm K})$')#, **axis_font)



    fig.savefig('Fig3_1D_aiat_2D_pT.png')
    plt.close()


############################################################
############################################################


def aiAT_boltzmann(system, num_points, enthalpy, entropy, temp_range, target_p, mu_N):
    
    fig,((fig_a, fig_b))= plt.subplots(1, 2)
    fig.subplots_adjust(left=0.15, bottom = 0.15, right =0.95, top = 0.90, wspace=0.5, hspace=0.5)
    fig.set_size_inches(10, 6) # (width, height) in inches
#    plt.axhline(0, color='black')

    dT = temp_range
    dp = target_p
    matrix_E = [];  p0 = 1
    dmu_N = mu_N + E_mlc_N2
    d_mu_O = np.linspace(mu_min2, mu_max2, num_points)
    mu_O = d_mu_O + E_mlc_O_ref + ( enthalpy - dT * entropy ) +  1/2. * kb * dT * np.log(dp/p0) 
    for a in system:
        Etot = a[1];	N_N = a[2];	N_O = a[3];	N_Ti = a[4];	Area = a[5]
        G  = (Etot - E_bulk_TiN * N_Ti - (N_N - N_Ti) * dmu_N - N_O * mu_O  ) / Area / 2. 

        E = np.exp(- G / kb / dT ) 
        matrix_E.append(E)

    sum_pop = 0
    for a in matrix_E:
	sum_pop = a + sum_pop

    i = 0 ; O_on_surf = 0
    for a in matrix_E:
        #print sum_pop[0:10]
        population = a/sum_pop * 100
        area = system[i][5] ; N_O = system[i][3];
	fig_a.plot(dT, population,label=system[i][0], c=system[i][6])
        O_on_surf = O_on_surf + population * ( N_O / area ) / 2 ## Due to symettric slab
	i = i + 1
	
    plt.legend()
    fig_a.set_xlabel(r'$T$ (K)')#, **axis_font)
    fig_a.set_xlim([0,800])
    fig_a.set_ylim([-10,110])
    fig_a.set_ylabel(r'Surface population (%)')#, **axis_font)
    legend_a = fig_a.legend(bbox_to_anchor=(.95, .95), loc=1, ncol=1, \
                            borderaxespad=0., scatterpoints=1, markerscale=7, fancybox=True)
    fig_b.plot(dT, O_on_surf, label='total O population')
    fig_b.set_ylabel(r'Relative Surface Oxidation (arb)')#, **axis_font)
    fig_b.set_xlabel(r'$T$ (K)')#, **axis_font)
    fig_b.set_xlim([0,800])
    fig_b.set_ylim([-2,22])
    fig_b.legend(bbox_to_anchor=(.95, .95), loc=1, ncol=1, \
                 borderaxespad=0., scatterpoints=1, markerscale=7, fancybox=True)
    fig.savefig('Fig4_boltzmann.png')
    plt.close()
    return O_on_surf


############################################################
############################################################


global E_bulk_TiN ; E_bulk_TiN = -19.42816              # eV/TiN
global E_bulk_Ti  ; E_bulk_Ti  = -15.674779/2           # eV/Ti atom
global E_mlc_N2   ; E_mlc_N2   = -16.633016/2           # eV/N  atom
global E_mlc_O3   ; E_mlc_O3   = -13.500444/3           # eV/O  atom
global E_mlc_O2   ; E_mlc_O2   = -8.7916228/2 - 0.9137131709999995   # eV/O  atom
global E_bulk_TiO2; E_bulk_TiO2= -.52983623E+02/2       # eV/TiO2
global E_mlc_O_ref; E_mlc_O_ref =  E_mlc_O2

global E_form     ; E_form  = E_bulk_TiN  - E_bulk_Ti - E_mlc_N2
global E_form2    ; E_form2 = (E_bulk_TiO2 - E_bulk_Ti - 2 * E_mlc_O_ref) /2
global kb         ; kb = 8.6173303 * 10**(-5)

global add_mu_range; add_mu_range = 0.5
#print E_mlc_O3 , E_mlc_O2

num_points = 500
bolzman_target_p = 5 #unit:torr
global mu_min1, mu_max1,mu_min2, mu_max2,dmu_N, dmu_O
mu_min1, mu_max1 = E_form - add_mu_range, add_mu_range
mu_min2, mu_max2 = E_form2- add_mu_range -2, add_mu_range

dmu_N, dmu_O = np.meshgrid(np.linspace(mu_min1, mu_max1, num_points), \
                           np.linspace(mu_min2, mu_max2, num_points))


    #   [      'Name',   Total E,   N_N, N_O, N_Ti, Area       , color]
				   
#symetric calculation
systems = [
	['(111):Ti',   -.84673034E+02,  4, 0, 5, 7.8189819348124, "#4d4d4d"],\
	['(111):Ti-O$_{0.8}$N$_{0.2}$', -.62316683E+03,  27, 8, 30, 7.8189819348124*5, 		"#ccccff"]#,\
  ['TiO$_{2}$', E_bulk_TiO2, 0, 2, 1, 0,"black"]
	]

#asymetric system	   
systems1 = [
	      ['(100)'  ,-76.917462, 4,   0,  4, 1*9.0245714, "#4d4d4d"],\
        ['(100):O$_{0.25}^{\\rm ad}$'  , -315.3454,  16,   1, 16, 4*9.0245714, "#ffcccc"],\
        ['(100):O$_{0.5}^{\\rm ad}$'  , -322.5520,   16,   2, 16, 4*9.0245714, "#ffb3b3"],\
        ['(100):O$_{0.75}^{\\rm ad}$'  , -328.6412,  16,   3, 16, 4*9.0245714, "#ff704d"],\
        ['(100):O$_{1.00}^{\\rm ad}$'  , -83.63019,   4,   1,  4, 1*9.0245714, "#ff3300"]
    ]




if __name__ == '__main__':
    from time import time
#    opts, args = command_line_arg()
    t0 = time()

    # Resolution: Total number of points to be considered =  num_points^2

    enthalpy, entropy, temp_range = read_janaf('O2_input_JANAF.txt', num_points, 1200, 2) #max_temperature
    aiAT_2D_plot(systems, num_points)
    aiAT_1D_2D_plot(systems, num_points, enthalpy, entropy, temp_range)
    

# mu_N = 0: N2 rich env
    O_on_surf_Nr = aiAT_boltzmann(systems, num_points, enthalpy, entropy, temp_range, bolzman_target_p/0.00131579, 0) #torr to atm
    os.system('mv Fig4_boltzmann.png Fig4_boltzmann_N_rich.png')

# mu_N = E_form: N2 lean env
    O_on_surf_Nl = aiAT_boltzmann(systems, num_points, enthalpy, entropy, temp_range, bolzman_target_p/0.00131579, E_form) #torr to atm
    os.system('mv Fig4_boltzmann.png Fig4_boltzmann_N_lean.png')


    fig,((fig_a))= plt.subplots(1, 1)
    fig.subplots_adjust(left=0.15, bottom = 0.15, right =0.95, top = 0.90, wspace=0.5, hspace=0.5)
    fig.set_size_inches(6, 6) # (width, height) in inches
#    plt.axhline(0, color='black')
    max_pop = 1 / 9.0245714 * 100
    fig_a.plot([min(temp_range),max(temp_range)],[max_pop,max_pop], '--', label='full 1ML', c='black')
    fig_a.plot(temp_range, O_on_surf_Nr,label='With N2 gas', c='blue')
    fig_a.plot(temp_range, O_on_surf_Nl,label='Without N2 gas',    c='black')
    plt.legend()
    fig_a.set_xlabel(r'$T$ (K)')#, **axis_font)
    fig_a.set_xlim([0,800])
    fig_a.set_ylim([-2,22])
    fig_a.set_ylabel(r'Relative Surface Oxidation (arb)')#, **axis_font)
    legend_a = fig_a.legend(bbox_to_anchor=(.95, .95), loc=1, ncol=1, \
                            borderaxespad=0., scatterpoints=1, markerscale=7, fancybox=True)
    fig.savefig('Fig5_pop_boltzmann.png')
    plt.close()    

    t1 = time()

    print '\nCode calc completed! Time Used: %.2f [sec]\n' % (t1 - t0)
  
    os.system('/usr/bin/gnome-open *.png &')


plt.show()
