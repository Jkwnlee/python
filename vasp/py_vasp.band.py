#!/usr/bin/env python

import os
import re
import numpy as np
from optparse import OptionParser

############################################################
__version__ = "1.0"
############################################################

def WeightFromPro(infile='PROCAR', whichAtom=None, spd=None, lsorbit=False):
    """
    Contribution of selected atoms to the each KS orbital
    """

    assert os.path.isfile(infile), '%s cannot be found!' % infile
    FileContents = [line for line in open(infile) if line.strip()]

    # when the band number is too large, there will be no space between ";" and
    # the actual band number. A bug found by Homlee Guo.
    # Here, #kpts, #bands and #ions are all integers
    nkpts, nbands, nions = [int(xx) for xx in re.sub('[^0-9]', ' ', FileContents[1]).split()]

    if spd:
        Weights = np.asarray([line.split()[1:-1] for line in FileContents
                              if not re.search('[a-zA-Z]', line)], dtype=float)
        Weights = np.sum(Weights[:,spd], axis=1)
    else:
        Weights = np.asarray([line.split()[-1] for line in FileContents
                              if not re.search('[a-zA-Z]', line)], dtype=float)
    
    nspin = Weights.shape[0] / (nkpts * nbands * nions)
    nspin /= 4 if lsorbit else 1

    if lsorbit:
        Weights.resize(nspin, nkpts, nbands, 4, nions)
        Weights = Weights[:,:,:,0,:]
    else:
        Weights.resize(nspin, nkpts, nbands, nions)
    
    if whichAtom is None:
        return np.sum(Weights, axis=-1)
    else:
        whichAtom = [xx - 1 for xx in whichAtom]
        return np.sum(Weights[:,:,:,whichAtom], axis=-1)



############################################################
def get_bandInfo(inFile = 'OUTCAR'):
    """
    extract band energies from OUTCAR
    """

    outcar = [line for line in open(inFile) if line.strip()]

    for ii, line in enumerate(outcar):
        if 'NKPTS =' in line:
            nkpts = int(line.split()[3])
            nband = int(line.split()[-1])

        if 'ISPIN  =' in line:
            ispin = int(line.split()[2])

        if "k-points in reciprocal lattice and weights" in line:
            Lvkpts = ii + 1

        if 'reciprocal lattice vectors' in line:
            ibasis = ii + 1

        if 'E-fermi' in line:
            Efermi = float(line.split()[2])
            LineEfermi = ii + 1
            # break

    # basis vector of reciprocal lattice
    B = np.array([line.split()[3:] for line in outcar[ibasis:ibasis+3]], 
                 dtype=float)
    # k-points vectors and weights
    tmp = np.array([line.split() for line in outcar[Lvkpts:Lvkpts+nkpts]],
                   dtype=float)
    vkpts = tmp[:,:3]
    wkpts = tmp[:,-1]

    # for ispin = 2, there are two extra lines "spin component..."
    N = (nband + 2) * nkpts * ispin + (ispin - 1) * 2
    bands = []
    # vkpts = []
    for line in outcar[LineEfermi:LineEfermi + N]:
        if 'spin component' in line or 'band No.' in line:
            continue
        if 'k-point' in line:
            # vkpts += [line.split()[3:]]
            continue
        bands.append(float(line.split()[1]))

    bands = np.array(bands, dtype=float).reshape((ispin, nkpts, nband))

    if os.path.isfile('KPOINTS'):
        kp = open('KPOINTS').readlines()

    if os.path.isfile('KPOINTS') and kp[2][0].upper() == 'L':
        Nk_in_seg = int(kp[1].strip())
        Nseg = nkpts // Nk_in_seg
        vkpt_diff = np.zeros_like(vkpts, dtype=float)
        
        for ii in range(Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            vkpt_diff[start:end, :] = vkpts[start:end,:] - vkpts[start,:]

        kpt_path = np.linalg.norm(np.dot(vkpt_diff, B), axis=1)
        # kpt_path = np.sqrt(np.sum(np.dot(vkpt_diff, B)**2, axis=1))
        for ii in range(1, Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            kpt_path[start:end] += kpt_path[start-1]

        # kpt_path /= kpt_path[-1]
        kpt_bounds =  np.concatenate((kpt_path[0::Nk_in_seg], [kpt_path[-1],]))
    else:
        # get band path
        vkpt_diff = np.diff(vkpts, axis=0)
        kpt_path = np.zeros(nkpts, dtype=float)
        kpt_path[1:] = np.cumsum(np.linalg.norm(np.dot(vkpt_diff, B), axis=1))
        # kpt_path /= kpt_path[-1]

        # get boundaries of band path
        xx = np.diff(kpt_path)
        kpt_bounds = np.concatenate(([0.0,], kpt_path[np.isclose(xx, 0.0)], [kpt_path[-1],]))

    return kpt_path, bands, Efermi, kpt_bounds

############################################################
def bandplot(kpath, bands, efermi, kpt_bounds, opts, whts=None):
    '''
    Use matplotlib to plot band structure
    '''

    width, height = opts.figsize
    ymin, ymax = opts.ylim
    dpi = opts.dpi

    fig = plt.figure()
    fig.set_size_inches(width, height)
    ax = plt.subplot(111)

    nspin, nkpts, nbands = bands.shape

    clrs = ['r', 'b']

    if opts.occLC and (whts is not None):
        from matplotlib.collections import LineCollection
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        LW = opts.occLC_lw
        DELTA = 0.3
        EnergyWeight = whts[0]
        norm = mpl.colors.Normalize(vmin=EnergyWeight.min(),
                                    vmax=EnergyWeight.max())
        # norm = mpl.colors.Normalize(0, 1)
        # create a ScalarMappable and initialize a data structure
        s_m = mpl.cm.ScalarMappable(cmap=opts.occLC_cmap, norm=norm)
        s_m.set_array([EnergyWeight])

        for Ispin in range(nspin):
            for jj in range(nbands):
                x = kpath
                y = bands[Ispin,:,jj]
                z = EnergyWeight[Ispin,:,jj]

                ax.plot(x, y,
                        lw=LW + 2 * DELTA,
                        color='gray', zorder=1)

                points = np.array([x, y]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = LineCollection(segments,
                                    # cmap=opts.occLC_cmap, # alpha=0.7,
                                    colors=[s_m.to_rgba(ww) for ww in (z[1:] + z[:-1])/2.]
                                    # norm=plt.Normalize(0, 1)
                                    )
                # lc.set_array((z[1:] + z[:-1]) / 2)
                lc.set_linewidth(LW)
                ax.add_collection(lc)

        divider = make_axes_locatable(ax)
        ax_cbar = divider.append_axes(opts.occLC_cbar_pos.lower(),
                size=opts.occLC_cbar_size, pad=opts.occLC_cbar_pad)

        if opts.occLC_cbar_pos.lower() == 'top' or opts.occLC_cbar_pos.lower() == 'bottom':
            ori = 'horizontal'
        else:
            ori = 'vertical'
        cbar = plt.colorbar(s_m, cax=ax_cbar,
                # ticks=[0.0, 1.0],
                orientation=ori)
        # cbar.ax.set_xticklabels([])
        # cbar.ax.set_yticklabels([])

    else:
        for Ispin in range(nspin):
            for Iband in range(nbands):
                lc = None if Iband == 0 else line.get_color()
                line, = ax.plot(kpath, bands[Ispin, :, Iband], lw=opts.linewidth, zorder=0,
                                alpha=0.8,
                                color=lc,
                                )
                if whts is not None:
                    for ii in range(len(opts.occ)):
                        ax.scatter(kpath, bands[Ispin, :, Iband],
                                color=opts.occMarkerColor[ii],
                                s=whts[ii][Ispin,:,Iband] * opts.occMarkerSize[ii],
                                marker=opts.occMarker[ii], zorder=1, lw=0.0,
                                alpha=0.5)

    for bd in kpt_bounds:
        ax.axvline(x=bd, ls='-', color='k', lw=0.5, alpha=0.5)

    ax.set_ylabel('Energy [eV]', # fontsize='small',
            labelpad=5)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(kpath.min(), kpath.max())

    ax.set_xticks(kpt_bounds)
    if opts.kpts:
        kname = [x.upper() for x in opts.kpts]
        for ii in range(len(kname)):
            if kname[ii] == 'G':
                kname[ii] = r'$\mathrm{\mathsf{\Gamma}}$'
            else:
                kname[ii] = r'$\mathrm{\mathsf{%s}}$' % kname[ii]
        ax.set_xticklabels(kname)
    else:
        ax.set_xticklabels([])

    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    plt.tight_layout(pad=0.20)
    plt.savefig(opts.bandimage, dpi=opts.dpi)

############################################################
def saveband_dat(kpath, bands):
    '''
    save band info to txt files
    '''
    prefix = 'pyband'
    spinSuffix = ['up', 'do']
    nspin, nkpts, nbands = bands.shape

    if nspin == 1:
        with open(prefix + '.dat', 'w') as out:
            for ii in range(nkpts):
                line = '%8.4f ' % kpath[ii]
                for jj in range(nbands):
                    line += '%10.4f' % bands[0,ii,jj]
                line += '\n'
                out.write(line)
    else:
        for Ispin in range(nspin):
            filename = prefix + '_' + spinSuffix[Ispin] + '.dat'
            with open(filename, 'w') as out:
                for ii in range(nkpts):
                    line = '%10.4f ' % kpath[ii]
                    for jj in range(nbands):
                        line += '%8.4f' % bands[0,ii,jj]
                    line += '\n'
                    out.write(line)

############################################################
def command_line_arg():
    usage = "usage: %prog [options] arg1 arg2"  
    par = OptionParser(usage=usage, version= __version__)

    par.add_option('-f', '--file',
            action='store', type="string",
            dest='filename', default='OUTCAR',
            help='location of OUTCAR')

    par.add_option('--procar', 
            action='store', type="string", dest='procar',
            default='PROCAR',
            help='location of the PROCAR')

    par.add_option('-z', '--zero',
            action='store', type="float",
            dest='efermi', default=None,
            help='energy reference of the band plot')

    par.add_option('-o', '--output',
            action='store', type="string", dest='bandimage',
            default='band.png',
            help='output image name, "band.png" by default')

    par.add_option('-k', '--kpoints',
            action='store', type="string", dest='kpts',
            default=None,
            help='kpoint path')

    par.add_option('-s', '--size', nargs=2,
            action='store', type="float", dest='figsize',
            default=(3.0, 4.0),
            help='figure size of the output plot')

    par.add_option('-y', nargs=2,
            action='store', type="float", dest='ylim',
            default=(-3, 3),
            help='energy range of the band plot')

    par.add_option('--lw',
            action='store', type="float", dest='linewidth',
            default=1.0,
            help='linewidth of the band plot')

    par.add_option('--dpi', 
            action='store', type="int", dest='dpi',
            default=360,
            help='resolution of the output image')

    par.add_option('--occ', 
            action='append', type="string", dest='occ',
            default=[],
            help='orbital contribution of each KS state')

    par.add_option('--occL', 
            action='store_true', dest='occLC',
            default=False,
            help='use Linecollection or Scatter to show the orbital contribution')

    par.add_option('--occLC_cmap', 
            action='store', type='string', dest='occLC_cmap',
            default='jet',
            help='colormap of the line collection')

    par.add_option('--occLC_lw', 
            action='store', type='float', dest='occLC_lw',
            default=2.0,
            help='linewidth of the line collection')

    par.add_option('--occLC_cbar_pos', 
            action='store', type='string', dest='occLC_cbar_pos',
            default='top',
            help='position of the colorbar')

    par.add_option('--occLC_cbar_size', 
            action='store', type='string', dest='occLC_cbar_size',
            default='3%',
            help='size of the colorbar, relative to the axis')

    par.add_option('--occLC_cbar_pad', 
            action='store', type='float', dest='occLC_cbar_pad',
            default=0.02,
            help='pad between colorbar and axis')

    par.add_option('--occM', 
            action='append', type="string", dest='occMarker',
            default=[],
            help='the marker used in the plot')

    par.add_option('--occMs', 
            action='append', type="int", dest='occMarkerSize',
            default=[],
            help='the size of the marker')

    par.add_option('--occMc', 
            action='append', type="string", dest='occMarkerColor',
            default=[],
            help='the color of the marker')

    par.add_option('--spd', 
            action='store', type="string", dest='spdProjections',
            default=None,
            help='Spd-projected wavefunction character of each KS orbital.')

    par.add_option('--lsorbit',
            action='store_true', dest='lsorbit',
            help='Spin orbit coupling on, special treament of PROCAR')

    par.add_option('-q', '--quiet', 
            action='store_true', dest='quiet',
            help='not show the resulting image')


    return  par.parse_args( )

############################################################
if __name__ == '__main__':
    opts, args = command_line_arg()

    if opts.occ:

        Nocc = len(opts.occ)
        occM  = ['o' for ii in range(Nocc)]
        occMc = ['r' for ii in range(Nocc)]
        occMs = [20  for ii in range(Nocc)]
        for ii in range(min(len(opts.occMarker), Nocc)):
            occM[ii] = opts.occMarker[ii]
        for ii in range(min(len(opts.occMarkerSize), Nocc)):
            occMs[ii] = opts.occMarkerSize[ii]
        for ii in range(min(len(opts.occMarkerColor), Nocc)):
            occMc[ii] = opts.occMarkerColor[ii]
        opts.occMarker = occM
        opts.occMarkerColor = occMc
        opts.occMarkerSize = occMs

        whts = []
        for occ in opts.occ:
            alist = np.array(occ.split(), dtype=int)
            nlist = [x for x in alist if not x == -1]
            cmark, = np.where(alist == -1)
            for ii in cmark:
                nlist += range(alist[ii + 1], alist[ii + 2] + 1)
            occAtom = set(nlist)

            # occAtom = [int(x) for x in opts.occ.split()]
            # whts = WeightFromPro(opts.procar, whichAtom = occAtom)
            if opts.spdProjections and (Nocc == 1):
                angularM = [int(x) for x in opts.spdProjections.split()]
                # print angularM
                whts.append(WeightFromPro(opts.procar, whichAtom=occAtom,
                    spd=angularM, lsorbit=opts.lsorbit))
            else:
                whts.append(WeightFromPro(opts.procar, whichAtom=occAtom,
                    lsorbit=opts.lsorbit))
    else:
        # if opts.spdProjections:
        #     angularM = [int(x) for x in opts.spdProjections.split()]
        #     whts = WeightFromPro(opts.procar, whichAtom=None, spd=angularM)
        # else:
        #     whts = None
        whts = None

    kpath, bands, efermi, kpt_bounds = get_bandInfo(opts.filename)
    if opts.efermi is None:
        bands -= efermi
    else:
        bands -= opts.efermi

    import matplotlib as mpl
    from matplotlib.ticker import AutoMinorLocator
    # Use non-interactive backend in case there is no display
    mpl.use('agg')
    import matplotlib.pyplot as plt
    mpl.rcParams['axes.unicode_minus'] = False

    bandplot(kpath, bands, efermi, kpt_bounds, opts, whts)
    saveband_dat(kpath, bands)


 
    os.system('/usr/bin/gnome-open band.png &')

#    if not opts.quiet:
#        from subprocess import call
#        call(['feh', '-xdF', opts.bandimage])
