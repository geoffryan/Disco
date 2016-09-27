import sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

xscale = "linear"
yscale = "log"
GAM = 1.66666666667
RMAX = 2.0

def loadCheckpoint(filename):

    f = h5.File(filename, "r")

    piph = f['Data']['Cells'][:,-1][...]
    prim = f['Data']['Cells'][:,:-1][...]
    index = f['Grid']['Index'][...]
    idPhi0 = f['Grid']['Id_phi0'][...]
    nphi = f['Grid']['Np'][...]
    t = f['Grid']['T'][0]
    riph = f['Grid']['r_jph'][...]
    ziph = f['Grid']['z_kph'][...]

    r = np.zeros(piph.shape)
    z = np.zeros(piph.shape)
    R = 0.5*(riph[1:] + riph[:-1])
    Z = 0.5*(ziph[1:] + ziph[:-1])
    for i in range(index.shape[0]):
        for k in range(index.shape[1]):
            ind0 = index[i,k]
            ind1 = ind0 + nphi[i,k]
            r[ind0:ind1] = R[i]
            z[ind0:ind1] = Z[k]

    if RMAX > 0:
        ind = r < RMAX
        r = r[ind]
        prim = prim[ind,:]

    return t, r, z, prim

def plotCheckpoint(file):
    
    print("Loading {0:s}...".format(file))

    t, r, z, prim = loadCheckpoint(file)

    rho = prim[:,0]
    P = prim[:,1]
    ur = prim[:,2]
    up = prim[:,3]
    uz = prim[:,4]
    Br = prim[:,5]
    Bp = prim[:,6]
    Bz = prim[:,7]

    try:
        Phi0 = prim[:,9]
        Phi1 = prim[:,10]
        Phi2 = prim[:,11]
    except IndexError:
        Phi0 = np.zeros(r.shape)
        Phi1 = np.zeros(r.shape)
        Phi2 = np.zeros(r.shape)

    print(prim.shape)

    B2 = Br*Br + Bp*Bp + Bz*Bz
    Bt = np.sqrt(Br*Br+Bp*Bp);
    #u2 = ur*ur + r*r*up*up + uz*uz
    u2 = ur*ur + up*up/(r*r) + uz*uz
    ut = np.sqrt(ur*ur + up*up/(r*r));
    w2 = 1.0 + u2
    #uB = ur*Br + r*up*Bp + uz*Bz
    uB = ur*Br + up*Bp/r + uz*Bz
    b2 = (B2 + uB*uB) / w2
    rhoh = rho + GAM/(GAM-1.0)*P
    s = np.log(P * np.power(rho, -GAM)) / (GAM-1.0)

    cs = np.sqrt(GAM*P/rhoh)
    cA = np.sqrt(b2/(rhoh+b2))
    Ma = np.sqrt(u2) / (cs/np.sqrt(1-cs*cs))

    print("   Plotting...")
    nq = prim.shape[1]



    fig, ax = plt.subplots(4,4,figsize=(16,12))
    plotAx(ax[0,0], z, rho, xscale, yscale, r"$z$", r"$\rho$", 'k+')
    plotAx(ax[0,1], z, P, xscale, yscale, r"$z$", r"$P$", 'k+')
    plotAx(ax[0,2], z, P/rho, xscale, yscale, r"$z$", r"$P/\rho$", 'k+')
    plotAx(ax[0,3], z, cs, xscale, yscale, r"$z$", r"$c_s$", 'k+')
    #if nq > 8:
    #    plotAx(ax[0,2], z, prim[:,8], xscale, "linear", r"$z$", r"$q$", 'k+')
    plotAx(ax[1,0], z, uz, xscale, "linear", r"$z$", r"$u_z$", 'k+')
    plotAx(ax[1,1], z, ut, xscale, yscale, r"$z$", r"$u_T$",'k+')
    #plotAx(ax[1,2], z, uz, xscale, "linear", r"$z$", r"$u_z$", 'k+')
    plotAx(ax[1,2], z, s, xscale, "linear", r"$z$", r"$s$", 'k+')
    plotAx(ax[1,3], z, Ma, xscale, yscale, r"$z$", r"$\mathcal{M}$", 'k+')
    plotAx(ax[2,0], z, Bz, xscale, "linear", r"$z$", r"$B^z$", 'k+')
    plotAx(ax[2,1], z, Bt, xscale, yscale, r"$z$", r"$B_T$",'k+')
    #plotAx(ax[2,2], z, Bz, xscale, "linear", r"$z$", r"$B^z$", 'k+')
    plotAx(ax[2,2], z, b2/P, xscale, yscale, r"$z$", r"$b^2/P$", 'k+')
    plotAx(ax[2,3], z, cA, xscale, yscale, r"$z$", r"$c_A$", 'k+')
    plotAx(ax[3,0], z, Phi0, xscale, "linear", r"$z$", r"$\Phi_0$", 'k+')
    plotAx(ax[3,1], z, Phi1, xscale, "linear", r"$z$", r"$\Phi_1$", 'k+')
    plotAx(ax[3,2], z, Phi2, xscale, "linear", r"$z$", r"$\Phi_2$", 'k+')

    title = "DISCO t = {0:.3g}".format(t)

    #fig.suptitle(title, fontsize=18)

    plt.tight_layout()

    name = file.split('/')[-1].split('.')[0].split('_')[-1]
    plotname = "plot_r_{0:s}.png".format(name)
    
    print("   Saving {0:s}...".format(plotname))
    fig.savefig(plotname)

    plt.close(fig)

def plotAx(ax, x, y, xscale, yscale, xlabel, ylabel, *args, **kwargs):
    ax.plot(x, y, *args, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    if (y>0).any():
        ax.set_yscale(yscale)
    else:
        ax.set_yscale("linear")

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Makes plots of Disco prims as a function of r.")
        print("usage: python plotDiscoR.py <checkpoint.h5 ...>")
        sys.exit()

    files = sys.argv[1:]
    for f in files:
        plotCheckpoint(f)
