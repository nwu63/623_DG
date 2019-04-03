import numpy as np

#-----------------------------------------------------------
def readMesh(fname):
    f = open(fname+'.gri', 'r')
    Nn, Ne, dim = [int(s) for s in f.readline().split()]
    # read vertices
    node = np.array([[float(s) for s in f.readline().split()] for n in range(Nn)])
    # read boundaries
    NB = int(f.readline())
    bdy = []
    for i in range(NB):
        s = f.readline().split()
        Nb = int(s[0])
        Bi = np.array([[int(s) for s in f.readline().split()] for n in range(Nb)])
        bdy.append({'nbface':Nb,'nnode':s[1],'name':s[2],'NB':Bi})
    # read elements
    Ne0 = 0; E2N = []
    while (Ne0 < Ne):
        s = f.readline().split(); ne = int(s[0])
        Ei = np.array([[int(s) for s in f.readline().split()] for n in range(ne)])
        E2N = Ei if (Ne0==0) else np.concatenate((E2N,Ei), axis=0)
        Ne0 += ne
    f.close()

    return node, E2N, bdy

def readCurvedMesh(fname):
    f = open(fname+'.gri', 'r')
    Nn, Ne, dim = [int(s) for s in f.readline().split()]
    # read vertices
    node = np.array([[float(s) for s in f.readline().split()] for n in range(Nn)])
    # read boundaries
    NB = int(f.readline())
    bdy = []
    for i in range(NB):
        s = f.readline().split()
        Nb = int(s[0])
        Bi = np.array([[int(s) for s in f.readline().split()] for n in range(Nb)])
        bdy.append({'nbface':Nb,'nnode':s[1],'name':s[2],'NB':Bi})
    # read elements
    E2N = []
    for i in range(2):
        s = f.readline().split(); ne = int(s[0])
        Ei = np.array([[int(s) for s in f.readline().split()] for n in range(ne)])
        E2N.append(Ei)
    q = f.readline().split();nq = int(q[0])
    qlist = np.zeros(nq,dtype=int)
    for i in range(nq):
        s = f.readline()
        qlist[i] = int(s)
    f.close()

    return node, E2N, bdy, qlist

def readMeshMatrices(fname):
    f = open(fname+'.dat', 'r')
    # read I2E
    l1 = [s for s in f.readline().split()]
    nI2E = int(l1[1])
    I2E = np.array([[int(s) for s in f.readline().split()] for n in range(nI2E)])

    # read B2E
    l2 = [s for s in f.readline().split()]
    nB2E = int(l2[1])
    B2E = np.array([[int(s) for s in f.readline().split()] for n in range(nB2E)])

    # read In
    l3 = [s for s in f.readline().split()]
    nIn = int(l3[1])
    In = np.array([[float(s) for s in f.readline().split()] for n in range(nIn)])

    # read I2E
    l4 = [s for s in f.readline().split()]
    nBn = int(l4[1])
    Bn = np.array([[float(s) for s in f.readline().split()] for n in range(nBn)])

    # read Area
    l5 = [s for s in f.readline().split()]
    nA = int(l5[1])
    Area = np.array([[float(s) for s in f.readline().split()] for n in range(nA)])

    f.close()
    return I2E, B2E, In, Bn, Area

def writeSolution(filename,**kwargs):
    np.savez(filename,**kwargs)


def readSolution(filename):
    npzfile = np.load(filename+'.npz')
    data = dict(npzfile)
    return data