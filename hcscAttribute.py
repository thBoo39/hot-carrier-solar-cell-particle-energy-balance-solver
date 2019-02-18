"""Defining hot carrier solar cell properties

Most attributes should be self explanetory.
If choose thermionic emission,
    self.Jext = self.JextTherm
    self.Uext = self.UextTherm
If choose tunneling
    self.Jext = self.JextESC
    self.Uext = self.UextESC

display_attributes(self) method in hcscAttribute class maybe useful when needs showing current condition.

This class will be inherited in hcsc class.
"""
import numpy as np
from scipy import constants as sc
# custom modules
import nu
import bulk
import contacts
import extraction_via_ESC as viaESC


class hcscAttribute:
    def __init__(self):
        """
        Attributes
        ----------
        Trm : float
            solar cell temperature
        tau_th : float
            carrier thermalization time
        tau_pp : float
            optical phonon decay time
        NM : float
            number of phonon modes
        Tph : float
            optical phonon temperature
        EP : bool
             if true, solve electron phonon balance equations
        d : float
            thickness of the cell if bulk structure is assumed
        Ephn : float
            dominant phonon energy (only used in polar optical (POP) scattering)
        hw0 : float
            dominant longitudinal optical (LO) phonon energy
        absb : class
            absorber class object defined in a separate file
        lcnt : class
            left electrical contact defined in a separte file
        rcnt : class
            right electrical contacts defined in a separete file
        lesc : class
            left energy selective contact defined in a separate file
        resc :class
            right energy selective contact defined in a separate file
        Jext : function
            current through energy selective contact
        Uext : function
            energy through energy selective contact
        Uth : function
            carrier energy thermalization rate
        _c1 : float
            constant variables
        _cJ : float
            constant variables
        _cS : float
            constant variables

        Methods
        -------
        display_attributes(self):
            print attributes
        JextTherm(muc, mue)
            current flux extraction via thermionic emission
        UextTherm(muc, mue)
            energy flux extraction via thermionic emission
        JextESC(self, muc, mue)
            extracted current through tunneling
            assume carrier selective contact
            electron current flows from absb to right contact (absb->resc->rcnt)
        UextESC(self, muc, mue)
            extracted energy flux
        JextRN(self, muc, mue)
            unused
        UextRN(self, muc, mue)
            unused
        Uthstandard(self, muc):
            standard thermalization
            consider dimension
            try 3D
            parameter is tau_th thermalization time

        """
        # Current and energy flux
        # Choose either via thermionic emission or tunneling
        # via tunneling
        self.Jext = self.JextESC
        self.Uext = self.UextESC
        # via thermionic emission
        self.Jext = self.JextTherm
        self.Uext = self.UextTherm
        # thermalization loss in bulk
        self.Uth = self.Uthstandard
        self.Trm = 300.0
        self.tau_th = 1.0*nu.ps
        self.tau_pp = 8.0*nu.ps
        self.NM = 1.0e17  # good guess number
        self.Tph = 300.0
        self.EP = False
        self.d = 100*nu.nm
        self.Ephn = 28.8*nu.meV
        self.hw0 = self.Ephn
        # For future purpose, materials in each components are defined
        # in separate files in bulk and contacts
        self.absb = bulk.bulk()
        # electrical contact
        self.lcnt = contacts.contact()
        self.rcnt = contacts.contact()
        # energy selective contacts left and right
        self.lesc = contacts.energy_selective_contact()
        self.resc = contacts.energy_selective_contact()
        # These are constants
        self._c1 = 2*np.pi/(sc.h**3*sc.c**2)  # Jrec
        self._cJ = sc.e*self.absb.m_e/(2*np.pi**2*sc.hbar**3)  # Je
        self._cS = self.absb.m_e/(2*np.pi**2*sc.hbar**3)  # Se
        return

    def display_attributes(self):
        print("Band gap (eV):", self.absb.Eg/nu.eV)
        print("ESC(right) extraction energy (eV):", self.resc.E/nu.eV)
        print("ESC(right) extraction energy width (eV):", self.resc.Ew/nu.eV)
        print("Thermalization time (ps):", self.tau_th/nu.ps)
        print("Electron phonon balance model (bool):", self.EP)
        print("Optical phonon decay time (ps):", self.tau_pp/nu.ps)
        print("Effective mass (electron mass):", self.absb.m_e/sc.m_e)
        return

    def JextTherm(self, muc, mue):
        """Single barrier ESC -> thermionic current
        """
        return -self._cJ*viaESC.J_therm_e_FD(self.resc.E,
                                             muc,
                                             mue,
                                             self.absb.T,
                                             self.rcnt.T)

    def UextTherm(self, muc, mue):
        """Energy flux extracted via thermionic emission
        """
        return self._cS*viaESC.S_therm_e_FD(self.resc.E,
                                            muc,
                                            mue,
                                            self.absb.T,
                                            self.rcnt.T)

    def JextESC(self, muc, mue):
        """current flux by tunneling or thermionic emission throu ESC
        assume carrier selective contact
        electron current flows from absb to right contact (absb->resc->rcnt)
        """
        return self._cJ*viaESC.J_e(self.resc.E,
                                   self.resc.Ew,
                                   self.absb.T,
                                   self.rcnt.T,
                                   muc,
                                   mue)

    def UextESC(self, muc, mue):
        """energy flux by tunneling or thermionic emission throu ESC
        """
        return self._cS*viaESC.S_e(self.resc.E,
                                   self.resc.Ew,
                                   self.absb.T,
                                   self.rcnt.T,
                                   muc,
                                   mue)

    def Uthstandard(self, muc):
        """
        bulk thermalization
        consider dimension
        try 3D
        parameter is tau_th thermalization time
        """
        ret2 = self.absb.density(muc)
        ret1 = 3*sc.k*(self.absb.T-self.Trm)*ret2*self.d
        ret = ret1/self.tau_th
        return ret


def main():
    hcsc = hcscAttribute()
    hcsc.display_attributes()
    return

if __name__ == '__main__':
    main()
