"""
J_therm_e_FD(mu, Efl, Efr, Tl, Tr)
    interface to hcsc.py
    Current flux through ESC via thermionic emission with Fermi Dirac
S_therm_e_FD(mu, Efl, Efr, Tl, Tr)
    interface to hcsc.py
    Energy flux through ESC via thermionic emission with Fermi Dirac
J_e(VB, ew, Tl, Tr, Efl, Efr)
    current flux via tunneling throu ESC
    VB extraction energy height
    Tl, Tr carrier temperature left and right
    Efl, Efr fermi energy left and right
S_e(VB, ew, Tl, Tr, Efl, Efr)
    energy flux via tunneling
j_therm_e_FD(Ez, mu, Efl, Efr, Tl, Tr)
    thermionic currrent
s_therm_e_FD(Ez, mu, Efl, Efr, Tl, Tr)
J_therm_e(VB, Tl, Tr, Efl, Efr)
    thermionic emission current
    Thermionic emission electron current and Tunneling
    Conduction band edge 0 at higher of the two
S_therm_e(VB, Tl, Tr, Efl, Efr)
    energy flux
    Thermionic emission electron current
    Conduction band edge 0 at higher of the two

j_tunnel_e_FD(Ez, mu, sig, Efl, Efr, Tl, Tr)
J_tunnel_e_FD(mu, sig, Efl, Efr, Tl, Tr)
s_tunnel_e_FD(Ez, mu, sig, Efl, Efr, Tl, Tr)
S_tunnel_e_FD(mu, sig, Efl, Efr, Tl, Tr)

J_tunnel_e0(self, mu, sig, Efl, Efr, Tl, Tr)
    tunneling current using step shape transmission coefficient
S_tunnel_e0(self, mu, sig, Efl, Efr, Tl, Tr)
    energy flux
    Conduction band edge 0 at higher of the two
J_tunnel_e(self, mu, sig, Efl, Efr, Tl, Tr)

TC(self, E, E0, Ew)
    transmission coefficient

RNf1(self, muc, Tc, mue)
    Ross and Nozik
    define two equations to solve
    solve for muc
    since use optimize muc in unit eV
RNf2(self, Tc, muc, mue)
RNf3(self, Tc, muc, mue)
    solve for Tc in particle conservation equation
    similar to f1 but solve for Tc

"""
from scipy import constants as sc
import numpy as np
from scipy import integrate


def j_therm_e_FD(Ez, mu, Efl, Efr, Tl, Tr):
    """thermionic currrent with Fermi Dirac distribution
    """
    Alr = sc.k*Tl*np.log(np.exp(-(Ez-Efl)/(sc.k*Tl))+1)
    Arl = sc.k*Tr*np.log(np.exp(-(Ez-Efr)/(sc.k*Tr))+1)
    # print Ez/nu.eV, Efl/nu.eV, Efr/nu.eV, Tl, Alr, Arl, Alr-Arl
    return (Alr-Arl)


def J_therm_e_FD(mu, Efl, Efr, Tl, Tr):
    a = mu
    b = mu+15*sc.k*Tl
    ret = integrate.quad(j_therm_e_FD, a, b,
                         args=(mu, Efl, Efr, Tl, Tr))
    return ret[0]


def s_therm_e_FD(Ez, mu, Efl, Efr, Tl, Tr):
    return (Ez*j_therm_e_FD(Ez, mu, Efl, Efr, Tl, Tr))


def S_therm_e_FD(mu, Efl, Efr, Tl, Tr):
    a = mu
    b = mu+15*sc.k*Tl
    ret = integrate.quad(s_therm_e_FD, a, b,
                         args=(mu, Efl, Efr, Tl, Tr))
    return ret[0]


def J_therm_e(VB, Tl, Tr, Efl, Efr):
    """thermionic emission current
    Thermionic emission electron current
    and Tunneling
    Conduction band edge 0 at higher of the two
    """
    kTl = sc.k*Tl
    kTr = sc.k*Tr
    ret = (kTl**2*np.exp(-(VB-Efl)/kTl) -
           kTr**2*np.exp(-(VB-Efr)/kTr))
    return ret


def S_therm_e(VB, Tl, Tr, Efl, Efr):
    """energy flux
    Thermionic emission electron current
    Conduction band edge 0 at higher of the two
    """
    kTl = sc.k*Tl
    kTr = sc.k*Tr
    # Jlr = A*kT*np.exp(-(sc.e*Vl))
    ret = (kTl**3*(VB/kTl+1)*np.exp(-(VB-Efl)/kTl) -
           kTr**3*(VB/kTr+1)*np.exp(-(VB-Efr)/kTr))
    return ret


# -----------------------------------------------------
def j_tunnel_e_FD(Ez, mu, sig, Efl, Efr, Tl, Tr):
    Alr = sc.k*Tl*np.log(np.exp(-(Ez-Efl)/(sc.k*Tl))+1)
    Arl = sc.k*Tr*np.log(np.exp(-(Ez-Efr)/(sc.k*Tr))+1)
    return (Alr-Arl)


def J_tunnel_e_FD(mu, sig, Efl, Efr, Tl, Tr):
    a = mu-(sig/2)
    b = mu+(sig/2)
    ret = integrate.quad(j_tunnel_e_FD, a, b,
                         args=(mu, sig, Efl, Efr, Tl, Tr))
    return ret[0]


def s_tunnel_e_FD(Ez, mu, sig, Efl, Efr, Tl, Tr):
    return (Ez*j_tunnel_e_FD(Ez, mu, sig, Efl, Efr, Tl, Tr))


def S_tunnel_e_FD(mu, sig, Efl, Efr, Tl, Tr):
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

def J_tunnel_e0(self, mu, sig, Efl, Efr, Tl, Tr):
    """tunneling current using step shape transmission coefficient
    """
    a = mu-sig/2
    b = mu+sig/2
    Blr = (np.exp(-a/(sc.k*Tl))-np.exp(-b/(sc.k*Tl)))
    Brl = (np.exp(-a/(sc.k*Tr))-np.exp(-b/(sc.k*Tr)))
    Jlr = Blr*(sc.k*Tl)**2*np.exp(Efl/(sc.k*Tl))
    Jrl = -Brl*(sc.k*Tr)**2*np.exp(Efr/(sc.k*Tr))
    ret = self._cJ*(Jlr+Jrl)
    # print a, b, Blr, sc.k*Tl, Tl
    return ret


def S_tunnel_e0(self, mu, sig, Efl, Efr, Tl, Tr):
    """energy flux
    Conduction band edge 0 at higher of the two
    """
    a = mu-sig/2
    b = mu+sig/2
    kTl = sc.k*Tl
    kTr = sc.k*Tr
    Blr = (a/kTl+1)*np.exp(-a/kTl)-(b/kTl+1)*np.exp(-b/kTl)
    Brl = (a/kTr+1)*np.exp(-a/kTr)-(b/kTr+1)*np.exp(-b/kTr)
    Slr = kTl**3*Blr*np.exp(Efl/kTl)
    Srl = -kTr**3*Brl*np.exp(Efr/kTr)
    # Slr = kTl**3*Blr
    # Srl = -kTr**3*Brl
    ret = self._cS*(Slr+Srl)
    return ret

def J_tunnel_e(self, mu, sig, Efl, Efr, Tl, Tr):
    El = mu-sig/2
    Eh = mu+sig/2
    a0 = sc.e*sc.m_e/(2*np.pi**2*sc.hbar**3)
    a1 = .5*np.sqrt(np.pi)*(np.sqrt(2*(sig)**2))

    def a2(T):
        return np.exp((2*sig**2-4*mu*sc.k*T)/(4*(sc.k*T)**2))

    def a3(Ez, T):
        return ((-2*mu*sc.k*T+2*sc.k*T*Ez+2*sig**2) /
                (2*sc.k*T*np.sqrt(2*sig**2)))

    Jlr = a2(Tl)*erf(a3(Eh, Tl))-a2(Tl)*erf(a3(El, Tl))
    Jrl = -(a2(Tr)*erf(a3(Eh, Tr))-a2(Tr)*erf(a3(El, Tr)))
    J_tunnel = (a0*a1*sc.k*Tl*(Jlr*np.exp(Efl/(sc.k*Tl)) +
                Jrl*np.exp(Efr/(sc.k*Tr))))
    return J_tunnel

# ---------------------------------------------------------
def JextRN(muc, mue):
    """dummy functions for RN model
    """
    return 0


def UextRN(muc, mue):
    """dummy function
    """
    return 0

# ---------------------------------------------------------
def TC(self, E, E0, Ew):
    """transmission coefficient
    """
    return Ew**2/4/(Ew**2/4+(E-E0)**2)

# ---------------------------------------------------------
def RNf1(self, muc, Tc, mue):
    """
    Ross and Nozik
    define two equations to solve
    solve for muc
    since use optimize muc in unit eV
    """
    self.absb.T = Tc
    try:
        ret1 = 2*self.resc.E*self.Jout(muc*nu.eV)/nu.eV
        ret2 = self.Uabs-self.Urec(muc*nu.eV)
        ret = ret1-ret2
    except FloatingPointError as e:
        print(e)
        print(traceback.format_exc())
        print('f1', muc/sc.e, self.absb.T, mue/sc.e)
        ret = -1
    # print 'Jext', ret2
    # print 'f1', muc, Tc, self.Jabs, -self.Jrec(muc*nu.eV), -ret2, ret
    return ret


def RNf2(self, Tc, muc, mue):
    self.absb.T = Tc
    try:
        Trt = self.rcnt.T
        # ret1 = self.Uabs-self.Urec(muc)
        # ret2 = self.Jout(muc)
        # ret3 = ret1/(ret2/nu.eV)
        ret3 = 2*self.resc.E
        ret4 = (ret3*(1-Trt/Tc)+2*muc*Trt/Tc)
        ret = -(ret4-2*mue)
        # print ret1, ret2, ret3/nu.eV, muc/nu.eV, Tc, ret
    except FloatingPointError as e:
        print(e)
        print(traceback.format_exc())
        print('f2', muc/sc.e, self.absb.T, mue/sc.e)
        ret = -1
    # print 'f2', self.Uabs, self.Urec(muc), self.Uext(muc, mue), ret
    # print 'Uext', self.Uext(muc, mue)
    return ret


def RNf3(self, Tc, muc, mue):
    """
    solve for Tc in particle conservation equation
    similar to f1 but solve for Tc
    """
    self.absb.T = Tc
    try:
        ret1 = 2*self.resc.E*self.Jout(muc)/nu.eV
        ret2 = self.Uabs-self.Urec(muc)
        ret = -(ret1-ret2)
    except FloatingPointError as e:
        print(e)
        print(traceback.format_exc())
        print('f3', muc/sc.e, self.absb.T, mue/sc.e)
        ret = -1
    return ret
