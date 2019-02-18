from __future__ import division
import numpy as np
from scipy import constants as sc
from scipy import integrate
import nu
from matplotlib import pyplot as plt
import math
# ------------------------------------------------------------
# default unit
# cm
# ------------------------------------------------------------
cdef double sc_k = sc.k


class bulk(object):
    """docstring for mat"""
    def __init__(self, T=300):
        super(bulk, self).__init__()
        # electron effective mass
        self.m_e = .1*sc.m_e
        # bandgap
        self.Eg = 0.65*sc.e
        # carrier temperature in the quantum well
        self.T = T
        # conduction band edge
        self.Ec = self.Eg/2
        # constant for 2d density of states
        self._cdos = 8*np.sqrt(2)*self.m_e**(1.5)/sc.h**3

    # carrier density at given Fermi level
    def f_d(self, E, muc):
        cdef double Ec = self.Ec

        return math.sqrt(E-self.Ec)/(math.exp((E-muc)/(sc_k*self.T))+1)

    def density(self, Efn):
        a = self.Ec
        b = 20*sc.k*self.T+a
        ret = integrate.quad(self.f_d, a, b, args=(Efn, ))
        # print ret
        return self._cdos*ret[0]


def main():
    return

if __name__ == '__main__':
    main()
