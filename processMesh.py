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