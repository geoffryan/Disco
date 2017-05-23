import math
import sys
import numpy as np
import checkpoint as cp

class DiscoDomain:

    # A class to contain a Disco Domain

    def __init__(self, checkpointFile=None, parFile=None):
        if checkpointFile is not None:
            self._loadFromCheckpoint(checkpointFile)
        elif parFile is not None:
            self._loadFromParfile(parFile)

    def __str__(self):

        nr = len(self.rjph)-1
        nz = len(self.zkph)-1

        optkeys = self._opts.keys()
        optkeys.sort()

        optStr = ", ".join(["{0:s}: {1:s}".format(
                        key, str(self._opts[key])) for key in optkeys])
        parStr = "NumR: {0:d}, NumZ: {1:d}".format(nr, nz)

        outStr = "DiscoDomain{{t: {0:.8g}, {1:s}, {2:s}}}".format(
                        self.t, parStr, optStr)

        return outStr

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

        self.git = git
        self._pars = pars
        self._opts = opts

        self.t = grid[0]
        self.rjph = rjph
        self.zkph = zkph

        self.piph = piph
        self.prim = prim
        self.bflux = bflux

        self.planets = planets
        self.diagRZ = diagRZ
        self.diagR = diagR
        self.diagZ = diagZ

    def _loadFromParfile(self, filename):

        return None


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: $ python domain.py <checkpoint.h5>")
        print("Builds DiscoDomain object from given checkpoint")

    dom = DiscoDomain(sys.argv[1])
    print(dom)
