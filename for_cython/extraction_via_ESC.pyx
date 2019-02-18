"""
JextTherm(muc, mue)
    Single barrier ESC -> thermionic current
UextTherm(muc, mue)

j_therm_e_FD(Ez, mu, Efl, Efr, Tl, Tr)
    thermionic currrent
J_therm_e_FD(mu, Efl, Efr, Tl, Tr)
s_therm_e_FD(Ez, mu, Efl, Efr, Tl, Tr)
S_therm_e_FD(mu, Efl, Efr, Tl, Tr)
J_therm_e(VB, Tl, Tr, Efl, Efr)
    thermionic emission current
    Thermionic emission electron current and Tunneling
    Conduction band edge 0 at higher of the two
S_therm_e(VB, Tl, Tr, Efl, Efr)
    energy flux
    Thermionic emission electron current
    Conduction band edge 0 at higher of the two
"""
from scipy import constants as sc
import numpy as np
from scipy import integrate
from math import exp


cdef double sc_k = sc.k


def J_therm_e_FD(double mu, double Efl, double Efr, double Tl, double Tr):
    cdef double a
    cdef double b
    a = mu
    b = mu+15*sc_k*Tl
    ret = integrate.quad(j_therm_e_FD, a, b,
                         args=(mu, Efl, Efr, Tl, Tr))
    return ret[0]


cdef double j_therm_e_FD(double Ez, double mu, double Efl, double Efr, double Tl, double Tr):
    """thermionic currrent with Fermi Dirac distribution
    """
    cdef double Alr
    cdef double Arl
    Alr = sc_k*Tl*np.log(exp(-(Ez-Efl)/(sc_k*Tl))+1)
    Arl = sc_k*Tr*np.log(exp(-(Ez-Efr)/(sc_k*Tr))+1)
    # print Ez/nu.eV, Efl/nu.eV, Efr/nu.eV, Tl, Alr, Arl, Alr-Arl
    return (Alr-Arl)


cdef double s_therm_e_FD(double Ez, double mu, double Efl, double Efr, double Tl, double Tr):
    return (Ez*j_therm_e_FD(Ez, mu, Efl, Efr, Tl, Tr))


def S_therm_e_FD(double mu, double Efl, double Efr, double Tl, double Tr):
    cdef double a
    cdef double b
    a = mu
    b = mu+15*sc_k*Tl
    ret = integrate.quad(s_therm_e_FD, a, b,
                         args=(mu, Efl, Efr, Tl, Tr))
    return ret[0]


def J_therm_e(double VB, double Tl, double Tr, double Efl, double Efr):
    """thermionic emission current
    Thermionic emission electron current
    and Tunneling
    Conduction band edge 0 at higher of the two
    """
    cdef double kT1
    cdef double kTr
    cdef ret
    kTl = sc_k*Tl
    kTr = sc_k*Tr
    ret = (kTl**2*exp(-(VB-Efl)/kTl) -
           kTr**2*exp(-(VB-Efr)/kTr))
    return ret


def S_therm_e(double VB, double Tl, double Tr, double Efl, double Efr):
    """energy flux
    Thermionic emission electron current
    Conduction band edge 0 at higher of the two
    """
    kTl = sc_k*Tl
    kTr = sc_k*Tr
    # Jlr = A*kT*exp(-(sc.e*Vl))
    ret = (kTl**3*(VB/kTl+1)*exp(-(VB-Efl)/kTl) -
           kTr**3*(VB/kTr+1)*exp(-(VB-Efr)/kTr))
    return ret


# -----------------------------------------------------
cdef double j_tunnel_e_FD(Ez, mu, sig, Efl, Efr, Tl, Tr):
    Alr = sc_k*Tl*np.log(exp(-(Ez-Efl)/(sc_k*Tl))+1)
    Arl = sc_k*Tr*np.log(exp(-(Ez-Efr)/(sc_k*Tr))+1)
    return (Alr-Arl)


cdef double J_tunnel_e_FD(mu, sig, Efl, Efr, Tl, Tr):
    a = mu-(sig/2)
    b = mu+(sig/2)
    ret = integrate.quad(j_tunnel_e_FD, a, b,
                         args=(mu, sig, Efl, Efr, Tl, Tr))
    return ret[0]


cdef double s_tunnel_e_FD(Ez, mu, sig, Efl, Efr, Tl, Tr):
    return (Ez*j_tunnel_e_FD(Ez, mu, sig, Efl, Efr, Tl, Tr))


cdef double S_tunnel_e_FD(mu, sig, Efl, Efr, Tl, Tr):
    a = mu-(sig/2)
    b = mu+(sig/2)
    ret = integrate.quad(s_tunnel_e_FD, a, b,
                         args=(mu, sig, Efl, Efr, Tl, Tr))
    return ret[0]


def J_e(VB, ew, Tl, Tr, Efl, Efr):
    """current flux by tunneling or thermionic emission throu ESC
    VB extraction energy height
    Tl, Tr carrier temperature left and right
    Efl, Efr fermi energy left and right
    """
    # ret = self.J_therm_e(VB, Tl, Tr, Efl, Efr)
    # ret = self.J_tunnel_e0(VB, ew, Efl, Efr, Tl, Tr)
    ret = J_tunnel_e_FD(VB, ew, Efl, Efr, Tl, Tr)
    # ret = ret1+ret2
    # current is negative for electron current
    return -ret


def S_e(VB, ew, Tl, Tr, Efl, Efr):
    """energy flux
    """
    # ret = self.S_therm_e(VB, Tl, Tr, Efl, Efr)
    # ret = self.S_tunnel_e0(VB, ew, Efl, Efr, Tl, Tr)
    ret = S_tunnel_e_FD(VB, ew, Efl, Efr, Tl, Tr)
    return ret
