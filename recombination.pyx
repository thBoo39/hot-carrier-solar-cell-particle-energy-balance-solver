import numpy as np
from scipy import constants as sc
import math


cdef double sc_k = sc.k
cdef double _c1 = 2*np.pi/(sc.h**3*sc.c**2)  # Jrec
# sc_k = sc.k
# _c1 = 2*np.pi/(sc.h**3*sc.c**2)


cpdef double frec(double E, double mu, double absb_T):
    """
    recombination flux at energy E
    note that symmetric two band model is used -> 2*mu
    generalized Planck law
    """
    cdef double denom, numer
    denom = _c1*E**2
    numer = (math.exp((E-2*mu)/(sc_k*absb_T))-1)
    # print 'denom, numer', denom, numer, E, mu, self.absb.T
    return denom/numer

def urec(E, mu, absb_T):
    """energy flux at specific energy E
    """
    # generalized Planck law
    return E*frec(E, mu, absb_T)
