from __future__ import division
from fortran import dg, getjacobian,basis2d,getrefmassmatrix,roeflux,eulerflux,basis2d,gbasis2d,integrate,getmassinv
import numpy as np
import matplotlib.pyplot as plt
from fileIO import readMesh, readMeshMatrices, writeSolution, readSolution
from processMesh import signedArea, computeCentroidVec
from constants import GAS_CONSTANT, GAMMA, getIC, getBC
import argparse
from plotting import plotSolution, plotCp

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
    parser.add_argument("--mesh", type=str, default='test')
    args = parser.parse_args()
    meshFile = args.mesh
    p = args.p
    # saveFile = '../solution_0/'+meshFile+'_'+str(order)+'_sol'
    # restartFile = '../solution_0/'+meshFile+'_1_sol'
    # restartFile = saveFile
    restart = 0 # 0=no, 1=yes but continue, 2=only plot
    writeSol = False


    node, E2N, bdy = readMesh('../../grid/'+meshFile)
    I2E, B2E, In, Bn, area = readMeshMatrices('../../grid/'+meshFile+'_mat')
    
    # I2E = np.array([[0,0,0,0]])
    # In = np.array([[0,0,0]])

    nelem = E2N.shape[0]
    nnode = node.shape[0]
    niface = I2E.shape[0]
    nbface = B2E.shape[0]

    # Minv = getmassinv(args.p,np.ones(2))
    # print(Mref)
    # print(Minv[0,:,:])
    # # print(np.linalg.inv(Mref) - Minv[0,:,:])
    # print(np.linalg.det(Mref))
    # exit()
    

    rBC = getBC()
    q = initSolution(p,restart=restart)
    CFL = 1/(1+p)
    convtol = 1e-7
    miniter = 1e3
    maxiter = 1e5

    q,resids,maxres,detJ = dg(q,p,node,E2N,I2E,B2E,In,Bn,rBC,GAMMA,GAS_CONSTANT,CFL,convtol,miniter,maxiter,nnode,nelem,niface,nbface)

    cl,cd,Es,cp,mach = integrate(q,p,B2E,Bn,rBC,detJ,GAMMA,GAS_CONSTANT,nelem,nbface)
    print(np.max(np.abs(resids)))
    # print(resids[:,:,0])
    # print(resids[:,:,1])
    # print(resids[:,:,2])
    # print(resids[:,:,3])
    print(cl,cd,Es)

    E2N -= 1
    B2E[0:1] -= 1
    plotSolution(node,E2N,mach)
    plotCp(node,E2N,B2E,cp)
    plt.show()

    # print(np.squeeze(resids))
    # print(np.max(np.abs(resids)))
    # Mref = getrefmassmatrix(0)
    # for ielem in range(nelem):
    #     _,_,J = getjacobian(node[E2N[ielem,:]-1,:])
    #     M = Mref * J
    #     print(np.sum(M)/area[ielem] - 1)