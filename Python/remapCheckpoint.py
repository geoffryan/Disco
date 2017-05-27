import sys
import discopy as dp

def remap(d1, d2):

    z = d2.Z
    r = d2.R
    nq = d1.numQ

    prim = []
    for k in range(d2.numZ):
        sheet = []
        for j in range(d2.numR):
            sheet.append(np.empty(d2.numPhi[k,j], nq))
        prim.append(sheet)

    for k in range(d2.numZ):
        for j in range(d2.numR):
            phi = d2.Phi(k, j)
            for i in range(d2.numPhi[k,j]):
                prim1 = d1.getPrimAt(r[j], phi[i], z[k])
                prim[k][j][i,:] = prim1[:]

    d2.prim = prim


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("usage: $ python remapCheckpoint.py <checkpoint.h5> <in.par> <filename>")
        print("Remaps checkpoint to grid given by in.par. Saves to file <filename>")
        sys.exit()

    chkin = sys.argv[1]
    parfile = sys.argv[2]
    chkout = sys.argv[3]

    d1 = dp.DiscoDomain(chkin)
    d2 = dp.DiscoDomain(parFile=parfile)
    remap(d1, d2)
    d2.writeCheckpoint(chkout)
