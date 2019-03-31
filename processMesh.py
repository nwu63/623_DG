import numpy as np
import scipy.sparse

def signedArea(node):
    mat = np.hstack((node,np.ones((3,1))))
    area = np.linalg.det(mat)/2
    return area
def testMesh(node,E2N,I2E,B2E,In,Bn):
    nelem = E2N.shape[0]
    niface = I2E.shape[0]
    nbface = B2E.shape[0]
    checkval = np.zeros((nelem,2))
    for iface in range(niface):
        elemL = I2E[iface,0]
        faceL = I2E[iface,1]
        elemR = I2E[iface,2]
        n1 = (faceL-1)%3
        n2 = (faceL-2)%3
        length = np.linalg.norm(node[E2N[elemL,n1],:] - node[E2N[elemL,n2],:])
        normal = In[iface,:2]
        if np.abs(length - In[iface,2]) > 1e-6:
            print('I2E Fail!')
        checkval[elemL,:] = checkval[elemL,:] + length*normal
        checkval[elemR,:] = checkval[elemR,:] - length*normal
    
    for iface in range(nbface):
        elem = B2E[iface,0]
        face = B2E[iface,1]
        n1 = (face-1)%3
        n2 = (face-2)%3
        length = np.linalg.norm(node[E2N[elem,n1],:] - node[E2N[elem,n2],:])
        normal = Bn[iface,:2]
        if np.abs(length - Bn[iface,2]) > 1e-6:
            print('B2E Fail!')
        checkval[elem,:] += length*normal
    print(np.max(np.linalg.norm(checkval,axis=1)))

def computeQuality(node,E2N,Area):
    nelem = E2N.shape[0]
    qual = np.zeros(nelem)
    for ielem in range(nelem):
        p1 = node[E2N[ielem,0],:]
        p2 = node[E2N[ielem,1],:]
        p3 = node[E2N[ielem,2],:]
        L1 = np.linalg.norm(p3-p2)
        L2 = np.linalg.norm(p3-p1)
        L3 = np.linalg.norm(p2-p1)
        qual[ielem] = 4*np.sqrt(3)*Area[ielem]/(L1**2 + L2**2 + L3**2)
    return qual


def computeCentroidVec(node,E2N):
    nelem = E2N.shape[0]
    centroidvec = np.zeros((nelem,3,2))
    for ielem in range(nelem):
        centroid = np.mean(node[E2N[ielem,:],:],axis=0)
        # print(centroid);exit()
        for iedge in range(3):
            n1 = (iedge+1)%3
            n2 = (iedge+2)%3
            p1 = node[E2N[ielem,n1],:]
            p2 = node[E2N[ielem,n2],:]
            midpoint = (p1+p2)/2
            centroidvec[ielem,iedge] = midpoint - centroid
    return centroidvec




def curveMesh(node,E2N,B2E,q):
    E2N -= 1 # convert to 0-based
    B2E[:,[0,1]] -= 1 
    nbface = B2E.shape[0]
    qlist = [] #high q list of elems
    nnode = node.shape[0]
    node_num = nnode + 1
    E2N_new = np.array([])
    for ib in range(nbface):
        if B2E[ib,2] == 4: # bottom wall
            elem = B2E[ib,0]
            face = B2E[ib,1]
            qlist.append(elem)
            elem_E2N = list(E2N[elem,:]) # copy elem row in E2N
            for iface in range(3):
                xy = getNode(node,E2N,elem,face) # node on face and elem
                new_xy = np.zeros((q+1,2))
                new_xy[:,0] = np.linspace(xy[0,0],xy[1,0],q+1) # new nodes on face + old nodes
                new_xy[:,1] = np.linspace(xy[0,1],xy[1,1],q+1)
                new_nodes = new_xy[2:-1,:]; # only new nodes
                if iface == face: # project to bump
                    new_xy[:,1] = bump(new_xy[:,0])
                node = np.vstack([node,new_xy]) # add node entry
                idx = (iface+1)%3
                if idx == 3:
                    idx = 0
                insert = idx + iface*(q-1)
                elem_E2N = [elem_E2N[0:insert],np.arange(node_num,node_num+q-2),elem_E2N[insert+1:]]
                node_num = node_num + q - 1
            if E2N_new.size == 0:
                E2N_new = elem_E2N
                print(elem_E2N)
            else:
                E2N_new = np.vstack([E2N_new,elem_E2N])
    E2N += 1
    E2N_new += 1
    qlist = np.array(qlist)
    return node,E2N,E2N_new,qlist

def getNode(node,E2N,elem,face):
    print(face)
    idx = E2N[elem,:]
    xy = np.zeros((2,2))
    j = 0
    idx2 = np.zeros(2)
    print(node.shape)
    for i in range(3):
        if i != face:
            xy[j,:] = node[idx[i],:]
            idx2[j] = i
            j = j + 1
    return xy


    
def bump(x):
    return 0.0625*np.exp(-25*x**2)