#!/bin/python
# Using code ase
# Purpose: Doping a poscar
# import os
import sys
import re
import numpy as np
import pandas as pd
from ase.io import read, write
from optparse import OptionParser
from optparse import Values as optparse_type
from ase import Atoms, neighborlist, Atom

############################################################
__version__ = 1.2
############################################################
# Define tetrahedron displacement in bcc, fcc, and hcp lattices
tetrahedron_displacement = {
    'bcc': [(0.25, 0.25, 0.25), (-0.25, -0.25, 0.25), (-0.25, 0.25, -0.25), (0.25, -0.25, -0.25)],
    'fcc': [(0.25, 0.25, 0), (0.25, 0, 0.25), (0, 0.25, 0.25), (-0.25, -0.25, 0)],
    'hcp': [(1/3, 2/3, 0.25), (2/3, 1/3, 0.25), (1/3, 2/3, -0.25), (2/3, 1/3, -0.25)]
}


octahedron_displacement = {
    'bcc': [(0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5), (-0.5, 0, 0), (0, -0.5, 0), (0, 0, -0.5)],
    'fcc': [(0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, 0, 0.5), (-0.5, 0, -0.5), (0, 0.5, 0.5), (0, -0.5, -0.5)],
    'hcp': [(2/3, 1/3, 0.5), (1/3, 2/3, -0.5), (2/3, 1/3, -0.5), (1/3, 2/3, 0.5), (1/3, 1/3, 0), (-1/3, -1/3, 0)]
}


def interstitial_doping(opts: optparse_type) -> Atoms:
    print(
        f'[CODE] You choose a mode to interstial doping mode: {opts.interstitial}')
    input_structure = read(opts.input_structure)

    n = 0  # replace with your atom index
    specific_atom_position = input_structure.get_positions()[n]
    if opts.interstitial[-2][0] == 't':
        print(
            f'[CODE] You choose an tetrahedron dict of {opts.interstitial[-1]}')
        displacement_dict = tetrahedron_displacement[opts.interstitial[-1]]
    elif opts.interstitial[-2][0] == 'o':
        print(
            f'[CODE] You choose an octahedron dict of {opts.interstitial[-1]}')
        displacement_dict = octahedron_displacement[opts.interstitial[-1]]
    specific_atom_position = input_structure.get_positions()[n]

    # Get the index of the target element in opts.interstitial
    indices = [i for i, e in enumerate(
        input_structure.get_chemical_symbols()) if e == opts.interstitial[0]]
    n = indices[int(opts.interstitial[1])-1]

    specific_atom_position = input_structure.get_positions()[n]

    new_positions = [specific_atom_position + input_structure.get_cell()
                     @ displacement for displacement in displacement_dict]

    cutoffs = neighborlist.natural_cutoffs(input_structure)
    nl = neighborlist.NeighborList(
        cutoffs, self_interaction=False, bothways=True)
    nl.update(input_structure)

    scores = []
    for position in new_positions:
        dummy = Atoms('H', positions=[position])
        temp_structure = input_structure + dummy
        cutoffs = neighborlist.natural_cutoffs(temp_structure)
        nl = neighborlist.NeighborList(
            cutoffs, self_interaction=False, bothways=True)
        nl.update(temp_structure)
        indices, offsets = nl.get_neighbors(len(temp_structure) - 1)

        # Check if there are any atoms of the same type at the new position
        same_type_atoms = np.array(temp_structure.get_chemical_symbols())[
            indices] == opts.interstitial[2]
        if any(same_type_atoms):
            score = 0
        else:
            score = np.sum(temp_structure.get_distances(
                len(temp_structure) - 1, indices))

        scores.append(score)

    selected_position = new_positions[np.argmax(scores)]
    new_Atom = Atom(opts.interstitial[2], tuple(selected_position))
    input_structure.append(new_Atom)
    input_structure = input_structure[np.argsort(
        input_structure.get_chemical_symbols())]

    write(opts.output, input_structure)
    return input_structure


def substitutional_doping(opts: optparse_type) -> Atoms:
    input_structure = read(opts.input_structure)
    positions = input_structure.get_positions()
    elements = input_structure.get_chemical_symbols()
    position = [[x[0], *x[1]] for x in list(zip(elements, positions))]

    new_position = []
    print(
        f'[CODE] You choose a mode to doping only - An-Atom \
{int(opts.select[0][1])} (st/nd/rd/th) - {opts.select[0][0]} to {opts.select[0][2]}')
    for n, i_position in enumerate(position):
        temp_position = []
        if i_position[0] == opts.select[0][0]:
            if n == int(opts.select[0][1]):
                temp_position = i_position
                temp_position[0] = opts.select[0][2]
                new_position.append(temp_position)
            else:
                new_position.append(i_position)
        else:
            new_position.append(i_position)

        new_atoms = Atoms([atom[0] for atom in new_position],
                          positions=[atom[1:] for atom in new_position])
        new_atoms.set_cell(input_structure.get_cell())
        write(opts.output, new_atoms)
    return new_atoms


def top_bot_constrained_doping(opts: optparse_type) -> Atoms:
    input_structure = read(opts.input_structure)
    positions = input_structure.get_positions()
    elements = input_structure.get_chemical_symbols()
    position = [[x[0], *x[1]] for x in list(zip(elements, positions))]

    new_position = []

    print(
        f'''[CODE] You choose a mode to doping Atoms ({opts.target_ele})
at top/bottom atoms to {opts.dopant_ele}''')
    colnames = ['e', 'x', 'y', 'z']
    pdf = pd.DataFrame(position, columns=colnames)
    zmin = pdf[pdf.e == opts.target_ele].z.min()
    zmax = pdf[pdf.e == opts.target_ele].z.max()
    pdf.loc[np.logical_and(pdf.e == opts.target_ele,
                           pdf.z == zmax), 'e'] = opts.dopant_ele
    pdf.loc[np.logical_and(pdf.e == opts.target_ele,
                           pdf.z == zmin), 'e'] = opts.dopant_ele
    new_position = pdf.values.tolist()
    new_atoms = Atoms([atom[0] for atom in new_position],
                      positions=[atom[1:] for atom in new_position])
    new_atoms.set_cell(input_structure.get_cell())
    write(opts.output, new_atoms)
    return new_atoms


def specific_region_doping(opts: optparse_type) -> Atoms:
    input_structure = read(opts.input_structure)
    positions = input_structure.get_positions()
    elements = input_structure.get_chemical_symbols()
    position = [[x[0], *x[1]] for x in list(zip(elements, positions))]

    new_position = []

    wmin, wmax = opts.window
    print(
        f'''[CODE] You choose a mode to doping Atoms ({opts.target_ele})
in specific region from {wmin} to {wmax} 
along z-axis to {opts.dopant_ele}''')
    for n, i_position in enumerate(position):
        temp_position = []
        if i_position[0] == opts.target_ele and \
                i_position[3] < wmax and \
                i_position[3] > wmin:
            temp_position = i_position
            temp_position[0] = opts.dopant_ele
            new_position.append(temp_position)
        else:
            new_position.append(i_position)

    new_atoms = Atoms([atom[0] for atom in new_position],
                      positions=[atom[1:] for atom in new_position])
    new_atoms.set_cell(input_structure.get_cell())
    write(opts.output, new_atoms)
    return new_atoms


def doping_cell(opts):
    if opts.output:
        pass
    else:
        opts.output = opts.input_structure+'_new.vasp'

    print("[CODE] doping the given structure  ", opts.input_structure)
    if len(opts.select) > 0:
        new_atoms = substitutional_doping(opts)
    elif len(opts.interstitial) > 0:
        cs_list = list(tetrahedron_displacement.keys()) + \
            list(octahedron_displacement.keys())
        cs_list = list(set(cs_list))
        opts.interstitial = list(opts.interstitial[0])
        if opts.interstitial[-1] in cs_list:
            new_atoms = interstitial_doping(opts)
        else:
            print(
                f"[CODE] the target crystal system ({opts.interstitial[-1]}) is not in {cs_list}")
            new_atoms = Atoms()
    elif opts.topbottom:
        new_atoms = top_bot_constrained_doping(opts)
    elif opts.target_ele and opts.dopant_ele:
        new_atoms = specific_region_doping(opts)
    else:
        new_atoms = Atoms()
        print('[CODE] YOU NEED TO CHECK help (-h) option')
    return new_atoms

############################################################
############################################################


def command_line_arg():
    usage = """
[CODE] run this commend with target input atomc structure (ex. POSCAR file) and doping which you want:
[CODE] for example: python %prog -i POSCAR -s Hf 3 Si
[CODE]              means: Change 3rd Hf atom to Si (the number of element is able to be found by VESTA)
    usage: %prog [options] arg1 arg2"""

    par = OptionParser(usage=usage, version=__version__)

    par.add_option("-i", '--input_structure',
                   action='store', type="string", dest='input_structure',
                   default='POSCAR',
                   help='location of the POSCAR')

    par.add_option("-o", '--output',
                   action='store', type="string", dest='output',
                   default=False,
                   help='Filename of the Doped_Structure')

    par.add_option("-a", '--target',
                   action='store', type="string", dest='target_ele',
                   default=False,
                   help='specify which atoms to substitute')

    par.add_option("-d", '--dopant',
                   action='store', type="string", dest='dopant_ele',
                   default=False,
                   help='specify which atoms to dope')

    par.add_option("--tb", '--topbottom',
                   action='store_true', dest='topbottom',
                   default=False,
                   help='specify which atoms to dope')

    par.add_option("-w", '--window', nargs=2,
                   action='store', type="float", dest='window',
                   default=(0, 0),
                   help='Minimul and maximum z-value where the target element in')

    par.add_option('-s', '--substitution', nargs=3,
                   action='append',  dest='select',
                   default=[],
                   help='''choose specific atom to others:
--select (-s) (target-element) (target-element number) (element-to-dope)
            ''')

    par.add_option('--interstitial', nargs=5,
                   action='append',  dest='interstitial',
                   default=[],
                   help='''choose specific atom to others:
--interstitial (-s) [target-element] [target-element number] [element-to-dope] [position: t(tetrahedron), o(octahedrond)] [crystal_system: bcc fcc hcp]
            ''')
    return par.parse_args()

############################################################


if __name__ == '__main__':
    from time import time
    opts, args = command_line_arg()
    # print(type(opts) == optparse_type)

    t0 = time()
    doping_cell(opts)

    t1 = time()
#    print('\nCode calc completed! Time Used: %.2f [sec]\n' % (t1 - t0))
