#!/usr/bin/env python3

import numpy as np
import os
import matplotlib.pyplot as plt

__version__ = "2.1"

# Plotting function
def set_plot_2d(ax, label_swch=True, legend_swch=True):
    """
    Configures a 2D plot with optional labels and legends.
    """
    if label_swch:
        ax.set_xlabel(r'$\Delta\mu_N$')
        ax.set_ylabel(r'$\Delta\mu_O$')
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    if legend_swch:
        legend = ax.legend(
            bbox_to_anchor=(1.05, 1),
            loc=2,
            ncol=1,
            borderaxespad=0.,
            scatterpoints=1,
            markerscale=7,
            fancybox=True
        )
        legend.get_frame().set_alpha(0.5)

# Function to read JANAF table
def read_janaf(filename, num_points, max_temperature, n_atom):
    """
    Reads JANAF thermodynamic data and fits it to polynomial functions.
    """
    temperature, entropy, enthalpy, tentropy = [], [], [], []

    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            temp = line.split()
            if i == 1:
                ground_enthalpy = float(temp[4])
            temperature.append(float(temp[0]))
            entropy.append(float(temp[2]) * 0.0103642688 / 1000. / n_atom)
            tentropy.append(float(temp[2]) * float(temp[0]) * -0.0103642688 / 1000. / n_atom)
            enthalpy.append((float(temp[4]) - ground_enthalpy) * 0.0103642688 / n_atom)

    temp_range = np.linspace(min(temperature), max_temperature, num_points)
    entropy_fit = np.poly1d(np.polyfit(temperature, entropy, 15))(temp_range)
    enthalpy_fit = np.poly1d(np.polyfit(temperature, enthalpy, 5))(temp_range)

    plt.plot(temperature, tentropy, 'o', label=r'$-T \times S(T)$ (eV/atom)')
    plt.plot(temperature, enthalpy, 'o', label=r'Enthalpy $H(T)$ (eV/atom)')
    plt.plot(temp_range, -entropy_fit * temp_range, '--', label='Fitted $-T*S(T)$')
    plt.plot(temp_range, enthalpy_fit, '--', label='Fitted $H(T)$')
    plt.legend()
    plt.xlabel('Temperature (K)')
    plt.ylabel('Energy (eV/atom)')
    plt.savefig('Fig1_Janaf.png')
    plt.close()

    return enthalpy_fit, entropy_fit, temp_range

# AIAT 2D plot function
def aiat_2d_plot(system, num_points):
    """
    Creates a 2D plot of surface stability regions.
    """
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.1, bottom=0.15, right=0.6, top=0.9, wspace=0.5, hspace=0.7)

    for phase in system:
        xy_coords = []
        for i in range(num_points):
            for j in range(num_points):
                # Assuming G calculations here...
                xy_coords.append([i, j])

        xy_coords = np.array(xy_coords)
        ax.scatter(xy_coords[:, 0], xy_coords[:, 1], s=4500 / num_points, label=phase[0], c=phase[6])

    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])
    set_plot_2d(ax)
    plt.savefig('Fig2_2D_aiat.png')
    plt.close()

# Main function
if __name__ == '__main__':
    from time import time

    t0 = time()

    num_points = 500
    max_temperature = 1200
    n_atom = 2

    systems = [
        ['p1', -84.673034, 4, 0, 5, 7.8189819348124, "#4d4d4d"],
        ['p2', -602.89072, 25, 8, 30, 7.8189819348124 * 5, "#ff704d"],
        ['p3', -623.16683, 27, 8, 30, 7.8189819348124 * 5, "#ccccff"],
        ['p4', -187.70645, 7, 4, 9, 9.0285826576, "black"]
    ]

    enthalpy, entropy, temp_range = read_janaf('O2_input_JANAF.txt', num_points, max_temperature, n_atom)
    aiat_2d_plot(systems, num_points)

    t1 = time()
    print(f'\nCode calc completed! Time Used: {t1 - t0:.2f} sec\n')
