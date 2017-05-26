import math
import sys
import numpy as np
import checkpoint as cp

class DiscoDomain:

    # A class to contain a Disco Domain

    git = None
    _pars = None
    _opts = None

    t = 0.0

    rjph = None
    zkph = None
    piph = None
    _index = None
    
    prim = None
    bflux = None
    planets = None
    
    diagRZ = None
    diagR = None
    diagZ = None

    def __init__(self, checkpointFile=None, parFile=None):
        if checkpointFile is not None:
            self._loadFromCheckpoint(checkpointFile)
        elif parFile is not None:
            self._loadFromParfile(parFile)

    def __str__(self):

        nr = self.numR
        nz = self.numZ

        optkeys = self._opts.keys()
        optkeys.sort()

        optStr = ", ".join(["{0:s}: {1:s}".format(
                        key, str(self._opts[key])) for key in optkeys])
        parStr = "NumR: {0:d}, NumZ: {1:d}".format(nr, nz)

        outStr = "DiscoDomain{{t: {0:.12g}, {1:s}, {2:s}}}".format(
                        self.t, parStr, optStr)

        return outStr

    def writeCheckpoint(self, filename):

        grid = self._packGridTuple()
        data = self._packDataTuple()
        checkpoint = (self.git, self._pars, self._opts, grid, data)
        cp.saveCheckpoint(checkpoint, filename)

    @property
    def numQ(self):
        return self.prim[0][0].shape[1]

    @property
    def numR(self):
        return self.rjph.shape[0]-1

    @property
    def numZ(self):
        return self.zkph.shape[0]-1

    @property
    def N(self):
        return self.numPhi.sum()

    @property
    def numPhi(self):
        return np.array([[self.piph[k][j].shape[0] for j in range(self.numR)]
                            for k in range(self.numZ)], dtype=np.int)

    @property
    def numBf(self):
        if self.bflux is None:
            return 0
        return self.bflux[0][0].shape[1]

    @property
    def R(self):
        return np.atleast_1d(0.5*(self.rjph[:-1] + self.rjph[1:]))

    @property
    def Z(self):
        return np.atleast_1d(0.5*(self.zkph[:-1] + self.zkph[1:]))

    @property
    def Rmax(self):
        return self.rjph.max()
    
    @property
    def Rmin(self):
        return self.rjph.min()
    
    @property
    def Zmax(self):
        return self.zkph.max()
    
    @property
    def Zmin(self):
        return self.zkph.min()

    @property
    def Phimax(self):
        return self._pars['Phi_Max']

    def _loadFromCheckpoint(self, filename):

        checkpoint = cp.loadCheckpointAll(filename)

        git = checkpoint[0]
        pars = checkpoint[1]
        opts = checkpoint[2]
        grid = checkpoint[3]
        data = checkpoint[4]

        rjph = grid[4]
        zkph = grid[5]
        index = grid[2]
        nphi = grid[3]

        prim, bflux, piph, planets, diagRZ, diagR, diagZ = \
                                        self._unpackDataTuple(grid, data)

        self.git = git
        self._pars = pars
        self._opts = opts

        self.t = grid[0]
        self.rjph = rjph
        self.zkph = zkph
        self._index = index

        self.piph = piph
        self.prim = prim
        self.bflux = bflux

        self.planets = planets
        self.diagRZ = diagRZ
        self.diagR = diagR
        self.diagZ = diagZ

    def _loadFromParfile(self, filename):

        return None

    def _unpackDataTuple(self, grid, data):

        rjph = grid[4]
        zkph = grid[5]
        index = grid[2]
        nphi = grid[3]

        nr = rjph.shape[0]-1
        nz = zkph.shape[0]-1

        flatpiph = data[2]
        flatprim = data[0]
        flatbflux = data[1]
        
        piph = []
        prim = []
        bflux = None
        if flatbflux is not None:
            bflux = []


        for k in range(nz):
            primsheet = []
            piphsheet = []
            fluxsheet = []
            for j in range(nr):
                n1 = index[k,j]
                n2 = index[k,j] + nphi[k,j]
                primsheet.append(flatprim[n1:n2,:])
                if bflux is not None:
                    fluxsheet.append(flatbflux[n1:n2,:])
                piphsheet.append(flatpiph[n1:n2])
            piph.append(piphsheet)
            if bflux is not None:
                bflux.append(fluxsheet)
            prim.append(primsheet)

        planets = data[3]
        diagRZ = data[4]
        diagR = data[5]
        diagZ = data[6]

        return (prim, bflux, piph, planets, diagRZ, diagR, diagZ)

    def _packGridTuple(self):

        idphi0 = self._getIndPhi(0.0)

        grid = (self.t, idphi0, self._index, self.numPhi, self.rjph, self.zkph)
        
        return grid

    def _packDataTuple(self):

        nr = self.numR
        nz = self.numZ
        N = self.N
        nq = self.numQ

        piph = np.empty((N,))
        prim = np.empty((N,nq))
        bflux = None
        if self.bflux is not None:
            nbf = self.numBf
            bflux = np.empty((N,nbf))

        nphi = self.numPhi

        for k in range(nz):
            for j in range(nr):
                n1 = self._index[k,j]
                n2 = self._index[k,j] + nphi[k,j]
                piph[n1:n2] = self.piph[k][j][:]
                prim[n1:n2,:] = self.prim[k][j][:,:]
                if bflux is not None:
                    bflux[n1:n2,:] = self.bflux[k][j][:,:]

        data = (prim, bflux, piph, self.planets, self.diagRZ, self.diagR,
                self.diagZ)

        return data

    def _getIndPhi(self, phi):
        nr = self.numR
        nz = self.numZ
        nphi = self.numPhi

        ind = np.empty((nz, nr), dtype=np.int)
        index = self._index

        phimax = self.Phimax

        for k in range(nz):
            for j in range(nr):
                ann = self.piph[k][j] - phi
                while ann.max() > 0.5*phimax:
                    ann[ann > 0.5*phimax] -= phimax
                while ann.min() < -0.5*phimax:
                    ann[ann < -0.5*phimax] += phimax
                ind[k,j] = index[k,j] + np.where(ann == ann[ann>0].min())[0]

        return ind


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: $ python domain.py <checkpoint.h5>")
        print("Builds DiscoDomain object from given checkpoint")

    dom = DiscoDomain(sys.argv[1])
    print(dom)

    print(dom.R)
    print(dom.Z)
    print(dom.numPhi)
    print(dom._pars)

    dom.writeCheckpoint("discopy_dom_test.h5")
