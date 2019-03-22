import numpy as np

GAMMA = 1.4
GAS_CONSTANT = 1.0
M_INF = 0.5
P_INF = 1.0
RHO_INF = 1.0
ALPHA = 0
def getIC():
    c = np.sqrt(GAMMA*P_INF/RHO_INF)
    u = M_INF*c*np.cos(np.deg2rad(ALPHA))
    v = M_INF*c*np.sin(np.deg2rad(ALPHA))
    q = np.zeros(4)
    q[0] = RHO_INF
    q[1] = RHO_INF*u
    q[2] = RHO_INF*v
    q[3] = P_INF/(GAMMA-1) + 0.5*RHO_INF*(u**2+v**2)
    return q

def getBC():
    Tt = 1 + (GAMMA-1)/2*M_INF**2
    pt = Tt**(GAMMA/(GAMMA-1))
    rBC = np.zeros(5)
    rBC[0] = P_INF
    rBC[1] = RHO_INF
    rBC[2] = Tt
    rBC[3] = pt
    rBC[4] = np.deg2rad(ALPHA) # we do the conversion here since FORTRAN does not handle degrees
    return rBC