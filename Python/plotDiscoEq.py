import os
import sys
import argparse as ag
import h5py as h5
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discoUtil as du

def primPlot(fig, ax, rjph, piph, r, q, label, vmin=None, vmax=None, noGhost=False, 
                colorbar=True, xlabel=None, ylabel=None, log=False, rmax=None):
        
    if vmin is None:
        vmin = q.min()
    if vmax is None:
        vmax = q.max()

    Rs = np.unique(r)
    if noGhost:
        Rs = Rs[:-2]

    if rmax is None or rmax <=0.0:
        lim_float = True
    else:
        lim_float = False

    if log:
        norm = mpl.colors.LogNorm(vmin, vmax)
    else:
        norm = mpl.colors.Normalize(vmin, vmax)

    for i, R in enumerate(Rs):
        ind = r==R
        phif = np.empty(len(piph[ind])+1)
        phif[1:] = piph[ind]
        phif[0] = piph[ind][-1]
        rf = rjph[i:i+2]

        x = rf[:,None] * np.cos(phif)[None,:]
        y = rf[:,None] * np.sin(phif)[None,:]

        C = ax.pcolormesh(x, y, q[None,ind], 
                cmap=plt.cm.inferno, vmin=vmin, vmax=vmax, norm=norm)

        if lim_float and rf.max() > rmax:
            rmax = rf.max()

        
    ax.set_aspect('equal')
    ax.set_xlim(-rmax, rmax)
    ax.set_ylim(-rmax, rmax)
    
    if colorbar:
        cb = fig.colorbar(C)
        cb.set_label(label, fontsize=24)

    if xlabel == None:
        xlabel = r'$x$'
    if ylabel == None:
        ylabel = r'$y$'

    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)


def plotCheckpoint(file, vars=None, logvars=None, noGhost=False, om=None,
                    bounds=None, rmax=None):
    
    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim, dat = du.loadCheckpoint(file)
    rjph = dat[0]
    zkph = dat[1]
    primPhi0 = dat[2]
    piph = dat[3]

    if om is not None:
        phi1 = phi - om*t
        piph1 = piph - om*t
    else:
        phi1 = phi
        piph1 = piph

    nq = prim.shape[1]

    varnames, vartex, num_c, num_n = du.getVarNames(file)

    Zs = np.unique(z)
    z_eq = Zs[len(Zs)/2]
    eq_ind = (z==z_eq)
    title = "DISCO t = {0:.1f}".format(t)
    name = file.split('/')[-1].split('.')[0].split('_')[-1]

    if vars is None:
        vars = range(nq)

    if logvars is None:
        logvars = []

    for q in range(nq):
        if q in vars or q in logvars:
            print("   Plotting...")

            if bounds is not None:
                vmin, vmax = bounds[q]
            elif q >= num_c:
                vmin, vmax = 0.0, 1.0
            else:
                vmin, vmax = None, None

            if q in vars:

                fig, ax = plt.subplots(1,1, figsize=(12,9))

                primPlot(fig, ax, rjph, piph1, r, prim[:,q], vartex[q], 
                            vmin=vmin, vmax=vmax, rmax=rmax)
                fig.suptitle(title, fontsize=24)
                plotname = "plot_eq_{0:s}_lin_{1:s}.png".format(name, varnames[q])
                
                print("   Saving {0:s}...".format(plotname))
                fig.savefig(plotname)
                plt.close(fig)

            if q in logvars:
                fig, ax = plt.subplots(1,1, figsize=(12,9))

                primPlot(fig, ax, rjph, piph1, r, prim[:,q], vartex[q], 
                            vmin=vmin, vmax=vmax, rmax=rmax, log=True)
                fig.suptitle(title, fontsize=24)
                plotname = "plot_eq_{0:s}_log_{1:s}.png".format(name, varnames[q])

                print("   Saving {0:s}...".format(plotname))
                fig.savefig(plotname)
                plt.close(fig)

def calcBounds(files):

    f = files[0]
    t, r, phi, z, prim, dat = du.loadCheckpoint(f)

    num_q = prim.shape[1]
    bounds = [[prim[:,q].min(), prim[:,q].max()] for q in range(num_q)]
    bounds = np.array(bounds)

    for f in files[1:]:
        t, r, phi, z, prim, dat = du.loadCheckpoint(f)
        for q in range(num_q):
            bounds[q,0] = min(bounds[q,0], prim[:,q].min())
            bounds[q,1] = max(bounds[q,1], prim[:,q].max())
    
    return bounds

def writeBoundsFile(filename, names, bounds):

    lines = ["{0:s} {1:.12g} {2:.12g}\n".format(names[q], bounds[q,0], bounds[q,1])
                for q in range(len(names))]

    f = open(filename, "w")
    for line in lines:
        f.write(line)

    f.close()

def readBoundsFile(filename, num_q):

    bounds = np.loadtxt(filename, usecols=[1,2])

    bounds = bounds[:num_q]

    return bounds
    
def getBounds(use_bounds, names, files):

    print(use_bounds)
    print(names)
    print(files)

    num_q = len(names)

    if use_bounds is not None:
        if use_bounds is True:
            bounds = calcBounds(files)
        else:
            if os.path.isfile(use_bounds):
                bounds = readBoundsFile(use_bounds, num_q)
            else:
                bounds = calcBounds(files)
                writeBoundsFile(use_bounds, names, bounds)
    else:
        bounds = None

    return bounds

if __name__ == "__main__":

    parser = ag.ArgumentParser(description="Create 2D plots of Disco variables.")
    parser.add_argument('checkpoints', nargs='+', 
                            help="Checkpoint (.h5) files to plot.")
    parser.add_argument('-v', '--vars', nargs='+', type=int,
                            help="Variables to plot.")
    parser.add_argument('-l', '--logvars', nargs='+', type=int,
                            help="Variables to plot logscale.")
    parser.add_argument('-b', '--bounds', nargs='?', const=True,
                            help="Use global max/min for bounds, optionally set by BOUNDS file.")
    parser.add_argument('-r', '--rmax', type=float, 
                            help="Set plot limits to RMAX.")
    parser.add_argument('-o', '--omega', type=float, 
                            help="Rotate frame at rate OMEGA.")
    parser.add_argument('--noghost', action='store_true', 
                            help="Do not plot ghost zones.")

    args = parser.parse_args()

    vars = args.vars
    logvars = args.logvars
    om = args.omega
    rmax = args.rmax
    use_bounds = args.bounds
    noghost = args.noghost

    files = args.checkpoints

    names, texnames, num_c, num_n = du.getVarNames(files[0])

    bounds = getBounds(use_bounds, names, files)

    for f in files:
        plotCheckpoint(f, vars=vars, logvars=logvars, bounds=bounds, om=om, 
                        rmax=rmax, noGhost=noghost)

