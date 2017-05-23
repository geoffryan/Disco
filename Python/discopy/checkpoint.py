import sys
import math
import h5py as h5
import numpy as np

def loadPars(f):
    pars = dict()

    parG = f['Pars']

    for key in parG:
        pars[key] = parG[key][0]

    return pars

def loadOpts(f):
    opts = dict()

    optG = f['Opts']

    for key in optG:
        opts[key] = optG[key][0]

    return opts

def loadGrid(f):
    g = f['Grid']
    id_phi0 = g['Id_phi0'][...]
    index = g['Index'][...]
    nphi = g['Np'][...]
    t = g['T'][0]
    r_jph = g['r_jph'][...]
    z_kph = g['z_kph'][...]

    grid = (t, id_phi0, index, nphi, r_jph, z_kph)

    return grid

def loadData(f):

    g = f['Data']
    ct_mode = f['Opts']['CT_MODE'][0]
    if ct_mode == 1:
        nfaces = 3
        bflux = g['Cells'][:,-1-nfaces:-1]
    elif ct_mode == 2:
        nfaces = 5
        bflux = g['Cells'][:,-1-nfaces:-1]
    else:
        nfaces = 0
        bflux = None

    prim = g['Cells'][:,:-1-nfaces]
    piph = g['Cells'][:,-1]
    planets = g['Planets'][...]
    diagRZ = g['Poloidal_Diagnostics'][...]
    diagR = g['Radial_Diagnostics'][...]
    diagZ = g['Vertical_Diagnostics'][...]

    data = (prim, bflux, piph, planets, diagRZ, diagR, diagZ)

    return data

def loadCheckpointAll(filename):

    f = h5.File(filename, "r")

    GIT_VERSION = f['GIT_VERSION'][0]
    pars = loadPars(f)
    opts = loadOpts(f)
    grid = loadGrid(f)
    data = loadData(f)

    f.close()

    checkpoint = (GIT_VERSION, pars, opts, grid, data)

    return checkpoint

def loadCheckpointData(filename):

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
    primPhi0 = np.zeros((index.shape[0], index.shape[1], prim.shape[1]))
    for k in range(index.shape[0]):
        for j in range(index.shape[1]):
            ind0 = index[k,j]
            ind1 = ind0 + nphi[k,j]
            r[ind0:ind1] = R[j]
            z[ind0:ind1] = Z[k]
            piph_strip = piph[ind0:ind1]
            pimh = np.roll(piph_strip, 1)
            pimh[pimh>piph_strip] -= 2*np.pi
            phi[ind0:ind1] = 0.5*(pimh+piph_strip)
            primPhi0[k,j,:] = prim[idPhi0[k,j],:]

    return t, r, phi, z, prim, (riph, ziph, primPhi0, piph)

def loadDiagRZ(filename):

    f = h5.File(filename, "r")

    t = f['Grid']['T'][0]
    rjph = f['Grid']['r_jph'][...]
    zkph = f['Grid']['z_kph'][...]
    diag = f['Data']['Poloidal_Diagnostics'][...]

    f.close()

    Nr = rjph.shape[0]-1
    Nz = zkph.shape[0]-1
    Nq = diag.shape[-1]
    #diag = np.resize(diag, (Nz,Nr,Nq))

    R = 0.5*(rjph[1:]+rjph[:-1])
    Z = 0.5*(zkph[1:]+zkph[:-1])

    r = np.empty((Nz,Nr))
    z = np.empty((Nz,Nr))

    r[:,:] = R[None,:]
    z[:,:] = Z[:,None]

    return t, r, z, diag, rjph, zkph

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: $ python checkpoint.py <checkpoint.h5>")
        print("Loads checkpoint data and prints to screen.")
        sys.exit()

    checkpoint = loadCheckpointAll(sys.argv[1])

    print(checkpoint[0])
    print(checkpoint[1])
    print(checkpoint[2])
    print(checkpoint[3])

