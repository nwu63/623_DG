import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation,UniformTriRefiner
from fileIO import readSolution
import numpy as np
from dg_solver import getx,getm
from matplotlib import rc
import sys
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
figsize = (6,3)
def plotSolution(q,node,qlist,E2N,E2N2,gamma,p,geom,level):
    fig = plt.figure(figsize=figsize)
    std_tri = refineTriangle(level)
    
    nelem = E2N.shape[0]
    xieta = np.stack([std_tri.x,std_tri.y],axis=1)
    n_xieta = xieta.shape[0]
    # xy = np.zeros((nelem,n_xieta,2))
    xy = []
    mach = []
    # mach = np.zeros((nelem,n_xieta))
    new_E2N = []
    for ielem in range(nelem):
        if np.any(qlist==ielem):
            idx = np.squeeze(np.argwhere(qlist == ielem))
            xy_elem = getx(node[E2N2[idx,:],:],geom,xieta)
        else:
            xy_elem = getx(node[E2N[ielem,:],:],1,xieta)
        xy.append(xy_elem)
        new_E2N.append(std_tri.triangles+n_xieta*ielem)
        mach.append(getm(q[ielem,:,:],p,xieta,gamma))
    xy = np.vstack(xy)
    mach = np.hstack(mach)
    new_E2N = np.vstack(new_E2N)
    triang = Triangulation(xy[:,0], xy[:,1], triangles=new_E2N)
    n_levels = 20
    plt.tricontourf(triang,mach,100)
    plt.tricontour(triang,mach,10,colors='k')
    # plt.triplot(triang, lw=0.5, color='white')

    plt.axis('equal')
    plt.title('Mach Number')
    plt.xlim(-1.5,1.5)
    plt.ylim(0,0.8)
    # plt.tight_layout()
    # plt.savefig('modern_art.pdf')

def refineTriangle(level):
    x = np.array([0,1,0])
    y = np.array([0,0,1])
    E2N = np.array([[0,1,2]])
    triang = Triangulation(x,y, triangles=E2N)
    refiner = UniformTriRefiner(triang)
    triang2 = refiner.refine_triangulation(subdiv=level)
    return triang2

def plotCp(node,E2N,B2E,cp):
    plt.figure(figsize=figsize)
    cp = np.squeeze(cp[np.where(cp[:,0]>0),:])
    elem = cp[:,0].astype(int)
    x = np.zeros_like(elem,dtype=float)
    nbface = B2E.shape[0]
    for iface in range(nbface):
        ielem = B2E[iface,0]
        face = B2E[iface,1]
        btype = B2E[iface,2]
        if btype == 4:
            idx = np.squeeze(np.argwhere(elem == ielem))
            n1 = (face-1)%3
            n2 = (face-2)%3
            xvals = node[E2N[elem[idx],[n1,n2]],0]
            x[idx] = np.mean(xvals)
    cp[:,0] = x
    cp = cp[cp[:,0].argsort()]
    plt.plot(cp[:,0],cp[:,1])
    plt.gca().invert_yaxis()
    plt.xlabel('$x$')
    plt.ylabel('$c_p$')
    # plt.tight_layout()
    # plt.savefig('../plots/cp_bump2_1.pdf')

def plotConvergence():
    meshes = ['bump'+str(x) for x in range(4)]
    orders = [1,2]
    cl = np.zeros((len(orders),len(meshes)))
    cd = np.zeros_like(cl)
    Es = np.zeros_like(cl)
    cl_err = np.zeros((len(orders),len(meshes)))
    cd_err = np.zeros_like(cl)
    Es_err = np.zeros_like(cl)
    ndof = np.zeros(len(meshes))
    cl_exact = 1.537095
    cd_exact = 2.94278E-6
    Es_exact = 0.0

    for i,order in enumerate(orders):
        for j,mesh in enumerate(meshes):
            saveFile = '../solution_0/'+mesh+'_'+str(order)+'_sol'
            sol = readSolution(saveFile)
            cl[i][j] = sol['cl']
            cl_err[i][j] = np.abs(sol['cl'] - cl_exact)
            cd[i][j] = sol['cd']
            cd_err[i][j] = np.abs(sol['cd'] - cd_exact)
            Es[i][j] = sol['Es']
            Es_err[i][j] = np.abs(sol['Es'] - Es_exact)
            ndof[j] = sol['q'].shape[0]

    p_CL = np.polyfit(np.log(np.sqrt(ndof)),np.log(cl_err).T,1)
    p_CD = np.polyfit(np.log(np.sqrt(ndof)),np.log(cd_err).T,1)
    p_Es = np.polyfit(np.log(np.sqrt(ndof)),np.log(Es_err).T,1)
    print(p_CL[0,:])
    print(p_CD[0,:])
    print(p_Es[0,:])

    plt.figure(figsize=figsize)
    plt.plot(np.sqrt(ndof),cl_err.T,'-o')
    plt.ylabel('$c_l$')
    plt.xlabel('$\sqrt{\mathrm{dof}}$')
    plt.gca().set_xscale("log")
    plt.gca().set_yscale("log")
    # plt.tight_layout()
    # plt.savefig('../plots/cl.pdf')

    plt.figure(figsize=figsize)
    plt.plot(np.sqrt(ndof),cd_err.T,'-o')
    plt.ylabel('$c_d$')
    plt.gca().set_xscale("log")
    plt.xlabel('$\sqrt{\mathrm{dof}}$')
    plt.gca().set_yscale("log")
    # plt.tight_layout()
    # plt.savefig('../plots/cd.pdf')

    plt.figure(figsize=figsize)
    plt.plot(np.sqrt(ndof),Es_err.T,'-o')
    plt.ylabel('$E_s$')
    plt.xlabel('$\sqrt{\mathrm{dof}}$')
    plt.gca().set_xscale("log")
    plt.gca().set_yscale("log")
    # plt.tight_layout()
    # plt.savefig('../plots/Es.pdf')

def plotResidual():
    saveFile = '../solution/bump0_2_sol'
    sol = readSolution(saveFile)
    plt.figure(figsize=figsize)
    plt.semilogy(sol['resnorm'])
    plt.xlabel('Iterations')
    plt.ylabel('Residual Norm')
    # plt.tight_layout()
    # plt.savefig('../plots/resnorm_2.pdf')



if __name__ == '__main__':
    # plotConvergence()
    # plotResidual()
    # plt.show()
    refineTriangle(1)