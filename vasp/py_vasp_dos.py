#!/usr/bin/env python

import os
import re
import numpy as np
from optparse import OptionParser

import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
mpl.use('agg')
import matplotlib.pyplot as plt

import matplotlib.colors as mcolors
from matplotlib.patches import Polygon

mpl.rcParams['axes.unicode_minus'] = False

############################################################
__version__ = "1.2"
############################################################

def WeightFromPro(infile='PROCAR', lsorbit=False):
    """
    Contribution of selected atoms to the each KS orbital
    """

    assert os.path.isfile(infile), '%s cannot be found!' % infile
    FileContents = [line for line in open(infile) if line.strip()]

    # when the band number is too large, there will be no space between ";" and
    # the actual band number. A bug found by Homlee Guo.
    # Here, #kpts, #bands and #ions are all integers
    nkpts, nbands, nions = [int(xx) for xx in re.sub('[^0-9]', ' ', FileContents[1]).split()]

    # Weights = np.asarray([line.split()[-1] for line in FileContents
    #                       if not re.search('[a-zA-Z]', line)], dtype=float)
    Weights = np.asarray([line.split()[1:-1] for line in FileContents
                          if not re.search('[a-zA-Z]', line)], dtype=float)

    kpt_weight = np.asarray([line.split()[-1] for line in FileContents if 'weight' in line], dtype=float)
    
    energies = np.asarray([line.split()[-4] for line in FileContents
                            if 'occ.' in line], dtype=float)

    nlmax = Weights.shape[-1]
    nspin = Weights.shape[0] / (nkpts * nbands * nions)
    nspin /= 4 if lsorbit else 1

    if lsorbit:
        Weights.resize(nspin, nkpts, nbands, 4, nions, nlmax)
        Weights = Weights[:,:,:,0,:,:]
    else:
        Weights.resize(nspin, nkpts, nbands, nions, nlmax)

    kpt_weight.resize(nspin, nkpts)
    energies.resize(nspin, nkpts, nbands)
    
    return energies, kpt_weight, Weights

############################################################
def gradient_fill(x, y, fill_color=None, ax=None, direction=1, **kwargs):
    """
    Plot a line with a linear alpha gradient filled beneath it.

    Parameters
    ----------
    x, y : array-like
        The data values of the line.
    fill_color : a matplotlib color specifier (string, tuple) or None
        The color for the fill. If None, the color of the line will be used.
    ax : a matplotlib Axes instance
        The axes to plot on. If None, the current pyplot axes will be used.
    Additional arguments are passed on to matplotlib's ``plot`` function.

    Returns
    -------
    line : a Line2D instance
        The line plotted.
    im : an AxesImage instance
        The transparent gradient clipped to just the area beneath the curve.
    """

    line, = ax.plot(x, y, **kwargs)
    if fill_color is None:
        fill_color = line.get_color()

    # print fill_color
    zorder = line.get_zorder()
    alpha = line.get_alpha()
    alpha = 1.0 if alpha is None else alpha

    z = np.empty((100, 1, 4), dtype=float)
    rgb = mcolors.colorConverter.to_rgb(fill_color)
    z[:,:,:3] = rgb
    if direction == 1:
        z[:,:,-1] = np.linspace(0, alpha, 100)[:,None]
    else:
        z[:,:,-1] = np.linspace(alpha, 0, 100)[:,None]

    xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(), y.max()
    im = ax.imshow(z, aspect='auto', extent=[xmin, xmax, ymin, ymax],
                   origin='lower', zorder=zorder)

    xy = np.column_stack([x, y])
    if direction == 1:
        xy = np.vstack([[xmin, ymin], xy, [xmax, ymin], [xmin, ymin]])
    else:
        xy = np.vstack([[xmin, ymax], xy, [xmax, ymax], [xmin, ymax]])
    clip_path = Polygon(xy, lw=0.0, facecolor='none', edgecolor='none', closed=True)
    ax.add_patch(clip_path)
    im.set_clip_path(clip_path)

    ax.autoscale(True)

    return line, im
############################################################
def lorentz_smearing(x, x0, sigma=0.03):
    '''
    Lorentz smearing of a Delta function.
    '''
    return sigma / np.pi / ((x-x0)**2 + sigma**2)

def gaussian_smearing(x, x0, sigma=0.05):
    '''
    Gaussian smearing of a Delta function.
    '''
    smear = np.zeros(x.size)
    condition = np.abs(x - x0) < (sigma * 5.) 
    smear[condition] = 1. / (np.sqrt(2*np.pi) * sigma) * np.exp(-(x[condition] - x0)**2 / (2*sigma**2))

    return smear

def gaussian_smearing_org(x, x0, sigma=0.05):
    '''
    Gaussian smearing of a Delta function.
    '''

    return 1. / (np.sqrt(2*np.pi) * sigma) * np.exp(-(x - x0)**2 / (2*sigma**2))

def generateDos(opts):
    '''
    generate dos
    '''

    ens, kptw, whts = WeightFromPro(opts.procar, opts.lsorbit)
    nspin, nkpts, nbands, nions, nlmax = whts.shape

    emin =  ens.min()
    emax =  ens.max()
    eran = emax - emin
    emin = emin - eran * opts.extra
    emax = emax + eran * opts.extra

    xen = np.linspace(emin, emax, opts.nedos)
    # tDOS = np.empty((opts.nedos, nspin))
    pDOS = []

    tdos_smear = np.empty((nspin, nkpts, nbands, opts.nedos))
    for IS in range(nspin):
        sign = 1 if IS == 0 else -1
        for Ik in range(nkpts):
            for Ib in range(nbands):
                x0 = ens[IS, Ik, Ib]
                # tdos_smear[IS, Ik, Ib] = sign * gaussian_smearing(xen, x0, opts.sigma) * kptw[IS, Ik]
                tdos_smear[IS, Ik, Ib] = sign * gaussian_smearing_org(xen, x0, opts.sigma) * kptw[IS, Ik]
                # tdos_smear[IS, Ik, Ib] = sign * lorentz_smearing(xen, x0, opts.sigma) * kptw[IS, Ik]
    tDOS = np.sum(tdos_smear, axis=(1, 2)).T

    # tdos_smear = np.empty((nspin, nkpts, nbands, opts.nedos))
    # for IS in range(nspin):
    #     sign = 1 if IS == 0 else -1
    #     for Ik in range(nkpts):
    #         tdos_smear[IS, Ik, :] = sign * gaussian_smearing_org(xen[np.newaxis, ...],
    #                                         ens[IS,Ik,:,np.newaxis], opts.sigma) * kptw[IS, Ik]
    # tDOS = np.sum(tdos_smear, axis=(1, 2)).T

    # tdos_smear = gaussian_smearing_org(xen[np.newaxis,...], ens[...,np.newaxis],
    #         opts.sigma) * kptw[..., np.newaxis, np.newaxis]
    # tDOS = np.sum(tdos_smear, axis=(1, 2)).T
    # if nspin == 2: tDOS[:,1] *= -1

    if len(opts.pdosAtom) != 0:

        factor = np.ones(len(opts.pdosAtom))
        for ii in range(min(len(opts.pdosAtom), len(opts.pdosFactor))):
            factor[ii] = opts.pdosFactor[ii]

        for ia, atoms in enumerate(opts.pdosAtom):
            # p = np.zeros((opts.nedos, nspin))

            alist = np.array(atoms.split(), dtype=int)
            nlist = [x for x in alist if not x == -1]
            cmark, = np.where(alist == -1)
            for ii in cmark:
                nlist += range(alist[ii + 1], alist[ii + 2] + 1)
            nlist = [x - 1 for x in set(nlist)]

            if ia <= len(opts.spdProjections) - 1:
                spdList = [int(x) for x in opts.spdProjections[ia].split()]
                pwhts = np.sum(whts[..., spdList], axis=-1)
            else:
                pwhts = np.sum(whts, axis=-1)

            pwhts = np.sum(pwhts[:,:,:, nlist], axis=-1)

            p = np.sum(pwhts[..., np.newaxis] * tdos_smear, axis=(1, 2)).T

            for IS in range(nspin):
                sign = 1 if IS == 0 else -1

                p[:,IS] = opts.pdosOffset * ia * sign + p[:,IS] * factor[ia]
                if ia == 0:
                    tDOS[:,IS] += sign * opts.pdosOffset * len(opts.pdosAtom)

            pDOS += [p]


    return xen, tDOS, pDOS

############################################################
def readDOSFromFile(opts):
    '''
    Read DOS info from file.

    the format of the DOS file:
    first line: ISPIN, NEDOS
    second line: labels of the dos
    next lines:
        if ISPIN = 1:
            Energy pDOS1 PDOS2 ... TotalDOS
        else:
            Energy pDOS1_up PDOS2_up ... TotalDOS_up pDOS1_down PDOS2_down ... TotalDOS_down
    '''

    inp = open(opts.dosFromFile).readlines()

    # the dos basic info
    nspin, nedos = [int(x) for x in inp[0].split()[1:]]
    labels = inp[1].split()[1:]
    # data
    DOS = np.array([line.split() for line in inp[2:] if line.strip()],
            dtype=float)
    NoPdos = (DOS.shape[1] - 1) / nspin - 1

    tDOS = np.empty((nedos, nspin))
    pDOS = []
    xen = DOS[:,0]
    for ii in range(nspin):
        tDOS[:,ii] = DOS[:, (ii + 1) * (NoPdos + 1)]
    for pp in range(NoPdos):
        tmp = []
        for ii in range(nspin):
            tmp += [DOS[:, (pp + 1) + ii * (NoPdos + 1)]]
        pDOS += [np.array(tmp).T]

    opts.nedos = nedos
    opts.pdosAtom = ['' for x in range(NoPdos)]
    opts.pdosLabel = labels

    return xen, tDOS, pDOS

def saveDOSToFile(opts, xen, tDOS, pDOS):
    '''
    save DOS info to file.

    the format of the DOS file:
    first line: ISPIN, NEDOS
    second line: labels of the dos
    next lines:
        if ISPIN = 1:
            Energy pDOS1 PDOS2 ... TotalDOS
        else:
            Energy pDOS1_up PDOS2_up ... TotalDOS_up pDOS1_down PDOS2_down ... TotalDOS_down
    '''

    nspin = tDOS.shape[1]
    nedos = tDOS.shape[0]
    NoPdos = len(pDOS)

    out = open(opts.dosToFile, 'w')
    out.write('# %5d %8d\n' % (nspin, nedos))
    labels = '# ' + ' '.join(opts.pdosLabel) + '\n'
    out.write(labels)

    for nn in range(nedos):
        line = '%8.4f ' % xen[nn]
        for ii in range(nspin):
            for p in pDOS:
                line += '%8.4f ' % p[nn,ii]
            line += '%8.4f ' % tDOS[nn,ii]
        line += '\n'
        out.write(line)

############################################################
def dosplot(xen, tdos, pdos, opts):
    '''
    Use matplotlib to plot band structure
    '''

    width, height = opts.figsize
    xmin, xmax = opts.xlim
    dpi = opts.dpi

    plt.style.use(opts.mpl_style)
    # DO NOT use unicode minus regardless of the style
    mpl.rcParams['axes.unicode_minus'] = False

    fig = plt.figure()
    fig.set_size_inches(width, height)
    ax = plt.subplot(111)

    LINES = []
    nspin = tdos.shape[1]

    plabels = []
    LWs = []
    LCs = []
    if opts.pdosAtom:
        plabels = ['p_%d' % ii for ii in range(len(opts.pdosAtom))]
        LWs = [0.5 for ii in range(len(opts.pdosAtom))]
        LCs = [None for ii in range(len(opts.pdosAtom))]

        for ii in range(min(len(opts.pdosAtom), len(opts.pdosLabel))):
            plabels[ii] = opts.pdosLabel[ii]
        for ii in range(min(len(opts.pdosAtom), len(opts.linewidth))):
            LWs[ii] = opts.linewidth[ii]
        for ii in range(min(len(opts.pdosAtom), len(opts.linecolors))):
            LCs[ii] = opts.linecolors[ii]

    plabels += ['total']

    xen -= opts.zero

    for ip, p in enumerate(pdos):
        for ii in range(nspin):
            fill_direction = 1 if ii == 0 else -1
            lc = LCs[ip] if ii == 0 else line.get_color()

            if opts.fill:
                line, im = gradient_fill(xen, p[:,ii], ax=ax, lw=LWs[ip],
                        color=lc,
                        direction=fill_direction)
            else:
                line, = ax.plot(xen, p[:,ii], lw=LWs[ip], alpha=0.6,
                        color=lc)
            if ii == 0: LINES += [line]

    if opts.showtotal:
        for ii in range(nspin):
            fill_direction = 1 if ii == 0 else -1
            lc = 'k' if ii == 0 else line.get_color()

            if opts.fill:
                line, im = gradient_fill(xen, tdos[:,ii], ax=ax, 
                        color=lc,
                        lw=0.5,
                        # zorder=-1,
                        direction=fill_direction,
                        )
            else:
                line, = ax.plot(xen, tdos[:,ii], color=lc,
                        lw=0.5, alpha=0.6)
            if ii == 0: LINES += [line]

    ax.set_xlabel('Energy [eV]', # fontsize='small',
            labelpad=5)
    ax.set_ylabel('DOS [arb. unit]', # fontsize='small',
            labelpad=10)
#    ax.tick_params(which='both', labelsize='small')

    ax.tick_params(which='x', labelsize='small')

    ax.set_xlim(xmin, xmax)
    if opts.ylim is not None:
        ymin, ymax = opts.ylim
        ax.set_ylim(ymin, ymax)
    else:
        ymin, ymax =  y.min(), y.max()

    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.plot([0,0], [ymin,ymax], dashes=[6,2],color='black',lw=0.5)#, alpha=0.6)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    #ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    #ax.set_ytick([])

    opts.pdosLabel = plabels
    ax.legend(LINES, plabels,
              loc=opts.legendloc, 
              fontsize='small',
              frameon=True,
              framealpha=0.6)
    if opts.DOStitle:
       plt.title(opts.DOStitle)
    plt.tight_layout(pad=0.50)
    plt.savefig(opts.dosimage, dpi=opts.dpi)

############################################################
def command_line_arg():
    usage = "usage: %prog [options] arg1 arg2"  
    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-i", '--input', 
            action='store', type="string", dest='procar',
            default='PROCAR',
            help='location of the PROCAR')

    par.add_option("-p", '--pdos', 
            action='append', type="string", dest='pdosAtom',
            default=[],
            help='specify which atoms to plot the pdos')

    par.add_option('--pdosoffset', 
            action='store', type="float", dest='pdosOffset',
            default=0.0,
            help='offset in pdos plot')

    par.add_option("-l", '--label', 
            action='append', type="string", dest='pdosLabel',
            default=[],
            help='label of the pdos')

    par.add_option("-t", '--title', 
            action='store', type="string", dest='DOStitle',
            default=False,
            help='title of the graph')

    par.add_option('--lloc', 
            action='store', type="string", dest='legendloc',
            default='upper right',
            help='legend location of dos plot')

    par.add_option('--fac', 
            action='append', type="float", dest='pdosFactor',
            default=[],
            help='scale factor of the pdos')

    par.add_option('-z', '--zero',
            action='store', type="float",
            dest='zero', default=0.0,
            help='energy reference of the band plot')

    par.add_option('--sigma',
            action='store', type="float",
            dest='sigma', default=0.02,
            help='smearing parameter, default 0.02')

    par.add_option('-n', '--nedos',
            action='store', type="int",
            dest='nedos', default=5000,
            help='number of point in DOS plot')

    par.add_option('-o', '--output',
            action='store', type="string", dest='dosimage',
            default='dos_h.png',
            help='output image name, "dos_h.png" by default')

    par.add_option('-s', '--size', nargs=2,
            action='store', type="float", dest='figsize',
            default=(4.8, 3.0),
            help='figure size of the output plot')

    par.add_option('-x', nargs=2,
            action='store', type="float", dest='xlim',
            default=(-10, 10),
            help='x limit of the dos plot')

    par.add_option('-y', nargs=2,
            action='store', type="float", dest='ylim',
            default=(0, 10),
            help='energy range of the band plot')

    par.add_option('-e',
            action='store', type="float", dest='extra',
            default=0.05,
            help='extra energy range of the band plot')

    par.add_option('--lw',
            action='append', type="float", dest='linewidth',
            default=[],
            help='linewidth of the band plot')

    par.add_option('--lc',
            action='append', type="string", dest='linecolors',
            default=[],
            help='linecolors of the band plot')

    par.add_option('--fill',
            action='store_true', dest='fill',
            default=True,
            help='fill under the DOS')

    par.add_option('--nofill',
            action='store_false', dest='fill',
            help='no fill under the DOS')

    par.add_option('--dpi', 
            action='store', type="int", dest='dpi',
            default=360,
            help='resolution of the output image')

    par.add_option('--tot', 
            action='store_true', dest='showtotal',
            default=True,
            help='show total dos')

    par.add_option('--notot', 
            action='store_false', dest='showtotal',
            help='not show total dos')

    par.add_option('--style', 
            action='store', type='string', dest='mpl_style',
            default='default',
            help='plot style of matplotlib. See "plt.style.available" for list of available styles.')

    par.add_option('--fromfile', 
            action='store', type='string', dest='dosFromFile',
            default=None,
            help='plot the dos contained in the file')

    par.add_option('--tofile', 
            action='store', type='string', dest='dosToFile',
            default=None,
            help='save DOS to file.')

    par.add_option('--spd', 
            action='append', type="string", dest='spdProjections',
            default=[],
            help='''
Spd-projected wavefunction character of each KS orbital.
    s orbital: 0\n
    py, pz, px orbital: 1 2 3\n
    dxy, dyz, dz2, dxz, dx2 orbital: 4 5 6 7 8 ''')

    par.add_option('--lsorbit',
            action='store_true', dest='lsorbit',
            help='Spin orbit coupling on, special treament of PROCAR')

    par.add_option('-q', '--quiet', 
            action='store_true', dest='quiet',
            help='not show the resulting image')

    return  par.parse_args( )

############################################################
def grep_e_fermi():
    with open('OUTCAR') as f:
        lines = f.readlines()
        for line in lines:
            if 'E-fermi' in line: e_fermi = line.split()[2]
    return e_fermi

	

if __name__ == '__main__':
    from time import time
    opts, args = command_line_arg()

    if opts.dosFromFile:
        xen, tdos, pdos = readDOSFromFile(opts)
    else:
        t0 = time()
        xen, tdos, pdos = generateDos(opts)
        t1 = time()
        print '\nDOS calc completed! Time Used: %.2f [sec]\n' % (t1 - t0)
    
    t0 = time()
    dosplot(xen, tdos, pdos, opts)
    t1 = time()

    print '\nDOS plot completed! Time Used: %.2f [sec]\n' % (t1 - t0)
    e_fermi = float(grep_e_fermi())
    print 'If you want to shift the picture with respect to fermi energy, %.6f, vasp.dos_h -z %.6f' \
           %(e_fermi, e_fermi)
    os.system('/usr/bin/gnome-open dos_h.png &')
    
    # save dos to file
    if opts.dosToFile:
        saveDOSToFile(opts, xen, tdos, pdos)

#    if not opts.quiet:
#        from subprocess import call
#        call(['feh', '-xdF', opts.dosimage])
