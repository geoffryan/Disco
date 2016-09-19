import sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

xscale = "log"
yscale = "log"
GAM = 1.66666666667
RMIN = 0.0
RMAX = np.inf
M = 1.0

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
    phi = np.zeros(piph.shape)
    R = 0.5*(riph[1:] + riph[:-1])
    Z = 0.5*(ziph[1:] + ziph[:-1])
    for i in range(index.shape[0]):
        for k in range(index.shape[1]):
            ind0 = index[i,k]
            ind1 = ind0 + nphi[i,k]
            r[ind0:ind1] = R[i]
            z[ind0:ind1] = Z[k]
            piph_strip = piph[ind0:ind1]
            pimh = np.roll(piph_strip, 1)
            pimh[pimh>piph_strip] -= 2*np.pi
            phi[ind0:ind1] = 0.5*(pimh+piph_strip)

    if RMAX > 0:
        R2 = r*r + z*z
        ind = (R2 < RMAX*RMAX) * (R2 > RMIN*RMIN)
        r = r[ind]
        phi = phi[ind]
        z = z[ind]
        prim = prim[ind,:]

    return t, r, phi, z, prim

def plotCheckpoint(file):
    
    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim = loadCheckpoint(file)

    rho = prim[:,0]
    P = prim[:,1]
    lr = prim[:,2]
    lp = prim[:,3]
    lz = prim[:,4]
    Br = prim[:,5]
    Bp = prim[:,6]/r
    Bz = prim[:,7]

    th = np.arctan2(r, z)
    R = np.sqrt(r*r+z*z)
    sinth = r/R
    costh = z/R

    try:
        Phi0 = prim[:,9]
        Phi1 = prim[:,10]
        Phi2 = prim[:,11]
    except IndexError:
        Phi0 = np.zeros(r.shape)
        Phi1 = np.zeros(r.shape)
        Phi2 = np.zeros(r.shape)

    print(prim.shape)

    al = 1.0 / np.sqrt(1.0 + 2*M/R)
    ber = 2*M/R * sinth / (1.0 + 2*M/R)
    bep = np.zeros(r.shape)
    bez = 2*M/R * costh / (1.0 + 2*M/R)
    gamrr = 1.0 + 2*M/R * sinth*sinth
    gamrp = np.zeros(r.shape)
    gamrz = 2.0*M/R * sinth*costh
    gampp = r*r
    gampz = np.zeros(r.shape)
    gamzz = 1.0 + 2*M/R * costh*costh
    igamrr = gamzz / (1.0 + 2*M/R)
    igamrp = np.zeros(r.shape)
    igamrz = -gamrz / (1.0 + 2*M/R)
    igampp = 1.0/(r*r)
    igampz = np.zeros(r.shape)
    igamzz = gamrr / (1.0 + 2*M/R)


    B2 = gamrr*Br*Br + 2*gamrp*Br*Bp + 2*gamrz*Br*Bz \
            + gampp*Bp*Bp + gampz*Bp*Bz + gamzz*Bz*Bz
    u2 = igamrr*lr*lr + 2*igamrp*lr*lp + 2*igamrz*lr*lz \
            + igampp*lp*lp + igampz*lp*lz + igamzz*lz*lz
    w2 = 1.0 + u2
    w = np.sqrt(w2)
    uB = lr*Br + lp*Bp + lz*Bz
    b2 = (B2 + uB*uB) / w2

    u0 = w / al
    l0 = -al*w + ber*lr + bep*lp + bez*lz
    ur = igamrr*lr + igamrp*lp + igamrz*lz - ber*u0
    up = igamrp*lr + igampp*lp + igampz*lz - bep*u0
    uz = igamrz*lr + igampz*lp + igamzz*lz - bez*u0
    uR = sinth*ur + costh*uz

    rhoh = rho + GAM/(GAM-1.0)*P
    s = np.log(P * np.power(rho, -GAM)) / (GAM-1.0)

    cs = np.sqrt(GAM*P/rhoh)
    cA = np.sqrt(b2/(rhoh+b2))
    Ma = np.sqrt(u2) / (cs/np.sqrt(1-cs*cs))
    Mdot = -R*R*rho*uR

    print("   Plotting...")
    nq = prim.shape[1]

    fig, ax = plt.subplots(4,4,figsize=(16,12))
    plotAx(ax[0,0], R, rho, xscale, yscale, r"$R$", r"$\rho$", 'k+')
    plotAx(ax[0,1], R, P, xscale, yscale, r"$R$", r"$P$", 'k+')
    plotAx(ax[0,2], R, P/rho, xscale, yscale, r"$R$", r"$P/\rho$", 'k+')
    plotAx(ax[0,3], R, cs, xscale, yscale, r"$R$", r"$c_s$", 'k+')
    plotAx(ax[1,0], R, uR, xscale, "linear", R"$R$", R"$u^R$", 'k+')
    plotAx(ax[1,1], R, up, xscale, "linear", R"$R$", R"$u^\phi$",'k+')
    plotAx(ax[1,2], R, s, xscale, yscale, R"$R$", R"$s$", 'k+')
    plotAx(ax[1,3], R, Ma, xscale, yscale, R"$R$", R"$\mathcal{M}$", 'k+')
    plotAx(ax[2,0], R, Br, xscale, "linear", R"$R$", R"$B^R$", 'k+')
    plotAx(ax[2,1], R, Bp, xscale, "linear", R"$R$", R"$B^\phi$",'k+')
    plotAx(ax[2,2], R, b2/P, xscale, yscale, R"$R$", R"$b^2/P$", 'k+')
    plotAx(ax[2,3], R, cA, xscale, yscale, R"$R$", R"$c_A$", 'k+')
    plotAx(ax[3,0], r, Mdot, xscale, "linear", r"$R$", r"$\dot{M} = R^2\rho u^R$", 'k+')

    title = "DISCO t = {0:.3g}".format(t)

    #fig.suptitle(title, fontsize=18)

    plt.tight_layout()

    name = file.split('.')[0].split('_')[-1]
    plotname = "plot_bondi_{0:s}.png".format(name)
    
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