from __future__ import division
from dg_solver import dg,roeflux,eulerflux,basis2d,gauss2d_pre,integrate
import numpy as np
import matplotlib.pyplot as plt
from fileIO import readMesh, readCurvedMesh, readMeshMatrices, writeSolution, readSolution
from processMesh import signedArea, curveMesh
from constants import GAS_CONSTANT, GAMMA, getIC, getBC
import argparse
from plotting import plotMach, plotCp
import time



def initSolution(p,restart=False,filename=None):
    if restart != 0:
        d = readSolution(filename)
        q = d['q']
    else:
        qIC = getIC()
        q = np.tile(qIC,(nelem,int((p+1)*(p+2)/2),1))
    return q
def testFlux():
    qIC = getIC()
    np.random.seed(0)
    qL = qIC
    qR = qIC * (np.random.rand()*0.1+1)
    nrm = np.array([1,0])
    nrm = nrm/np.linalg.norm(nrm)
    print(roeflux(qL,qL,nrm,GAMMA)[0] - eulerflux(qL,nrm,GAMMA)[0])
    print(roeflux(qL,qR,nrm,GAMMA)[0] + roeflux(qR,qL,-nrm,GAMMA)[0])
    print(roeflux(qL,qR,nrm,GAMMA)[0] - eulerflux(qL,nrm,GAMMA)[0])

def testmat():
    Mref = getmassmatrix(0)
    for ielem in range(nelem):
        _,_,J = getjacobian(node[E2N[ielem,:]-1,:])
        M = Mref * J
        print(np.sum(M)/area[ielem] - 1)
    


def test_basis():
    xy = np.random.rand(10,2)
    p = 4
    phi = basis2d(xy, p)
    gphi = gbasis2d(xy, p)
    print(np.sum(phi,axis=1))
    print(np.sum(gphi,axis=1))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, default=0)
    parser.add_argument("--q", type=int, default=-1)
    parser.add_argument("--mesh", type=str, default='test')
    parser.add_argument("--task", choices=['run', 'post'], default='run')
    args = parser.parse_args()
    meshFile = args.mesh
    p = args.p
    if args.q == -1:
        if args.p >= 2:
            geom = 2
        else:
            geom = args.p + 1
    else:
        geom = args.q
    
    saveFile = '../solution/'+meshFile+'_'+str(args.p)+'_sol'
    # restartFile = '../solution_0/'+meshFile+'_1_sol'
    restartFile = saveFile

    I2E, B2E, In, Bn, area = readMeshMatrices('../../grid/'+meshFile+'_mat')
    if geom == 1:
        node, E2N, bdy = readMesh('../../grid/'+meshFile+'_'+str(geom))
        E2N = [E2N, np.zeros((1,3))]
        qlist = np.array([])
    elif geom >= 1:
        node, E2N, bdy,qlist = readCurvedMesh('../../grid/'+meshFile+'_'+str(geom))
        
    nelem = E2N[0].shape[0]
    nnode = node.shape[0]
    niface = I2E.shape[0]
    nbface = B2E.shape[0]

    rBC = getBC()
    q = initSolution(p,restart=False)
    CFL = 1/(1+p)
    convtol = 1e-7
    miniter = 1e3
    maxiter = 1e6
    t = time.time()
    if args.task == 'run':
        print('        |-------------------------------------------------|')
        print('        |        Running DG Solver on mesh',args.mesh,'         |')
        print('        |        p = ',p,'                                  |')
        print('        |        q = ',geom,'                                  |')
        print('        |-------------------------------------------------|\n')
        q,resids,resnorm = dg(q,p,geom,node,qlist,E2N[0],E2N[1],I2E,B2E,In,Bn,rBC,GAMMA,GAS_CONSTANT,CFL,convtol,miniter,maxiter)
        t2 = time.time() - t
    else:
        d = readSolution(restartFile)
        q = d['q']
        resids = d['resids']
        resnorm = d['resnorm']
        t2 = d['time']
    cl,cd,Es = integrate(q,p,geom,node,qlist,E2N[0],E2N[1],B2E,Bn,rBC,GAMMA,GAS_CONSTANT)
    print(cl,cd,Es)
    if args.task == 'run':
        writeSolution(saveFile,q=q,p=p,geom=geom,resids=resids,resnorm=resnorm,time=t2,cl=cl,cd=cd,Es=Es)
    else:
        qlist -= 1
        E2N[0] -= 1
        E2N[1] -= 1
        B2E[:,0:2] -= 1
        plotMach(q,node,qlist,E2N[0],E2N[1],GAMMA,p,geom,4)
        plotCp(q,node,qlist,E2N[0],E2N[1],B2E,rBC,GAMMA,p,geom,12)
        plt.show()