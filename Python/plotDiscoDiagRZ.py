import sys
import math
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import discoUtil as du

xscale = "log"
yscale = "log"
GAM = 1.66666666667
RMAX = 0.8
GR = False

def plotCheckpoint(file):
    
    print("Loading {0:s}...".format(file))

    t, r, z, diag, rjph, zkph = du.loadDiagRZ(file)

    nq = diag.shape[-1]

    name = ".".join(file.split("_")[-1].split(".")[:-1])

    print("   Plotting")

    labels = [r"$\rho$", r"$P$", r"$v^r$", r"$v^\phi$", r"$v^z$",
                r"$B^r$", r"$B^\phi$", r"$B^z$", r"$D$", r"$E$",
                r"$\rho v^r$", r"$\rho v^\phi$", r"$\rho v^z$",
                r"$\rho v^r v^\phi$", r"$B^r B^\phi$", r"$B^2$", r"$c_s$",
                r"$c_A$"]


    for q in range(nq):
        f = diag[:,:,q]

        fig, ax = plt.subplots(1,1)
        C = ax.pcolormesh(rjph, zkph, f, cmap=mpl.cm.inferno)
        cb = fig.colorbar(C)
        ax.set_xlabel(r'$r$')
        ax.set_ylabel(r'$z$')
        cb.set_label(labels[q])
        ax.set_aspect('equal')
        fig.savefig("plot_diag_rz_{0:s}_{1:02d}.png".format(name,q))
        plt.close(fig)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Makes plots of Disco prims as a function of r.")
        print("usage: python plotDiscoR.py <checkpoint.h5 ...>")
        sys.exit()

    files = sys.argv[1:]
    dat = np.empty((len(files), 4))
    t = np.empty(len(files))
    for i,f in enumerate(files):
        ret = plotCheckpoint(f)
        #t[i] = ret[0]
        #dat[i,:] = np.array(ret[1])
"""
    Rrp = dat[:,0]
    Mrp = dat[:,1]
    Pg = dat[:,2]
    Pb = dat[:,3]

    fig, ax = plt.subplots(2,3)
    ax[0,0].plot(t, Rrp, label=r"$R_{r\phi}$")
    ax[0,0].plot(t, Mrp, label=r"$M_{r\phi}$")
    ax[0,0].legend()
    ax[0,1].plot(t, Pg, label=r"$P_g$")
    ax[0,1].plot(t, Pb, label=r"$P_b$")
    ax[0,1].legend()
    ax[1,0].plot(t, (Rrp+Mrp)/Pg, label=r"$\alpha$")
    ax[1,0].plot(t, Mrp/Pg, label=r"$\alpha_M$")
    ax[1,0].legend()
    ax[1,1].plot(t, Pg/Pb, label=r"$\beta$")
    ax[1,1].legend()
    ax[1,2].plot(t, 180.0/np.pi * 0.5*np.arcsin(Mrp / Pb), label=r"$\theta_B$")
    ax[1,2].legend()
    fig.tight_layout()

    name = "mri"
    fig.savefig(name+".pdf")
"""

