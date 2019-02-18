""" Hot Carrier Solar Cell

Platform: Python 3.7
"""
import numpy as np
from scipy import constants as sc
from scipy import interpolate
from scipy import integrate
from scipy import optimize
from scipy.special import erf
from sys import exit
import traceback
import time
import os
import warnings
from matplotlib import pyplot as plt
# ------------------------------------------------------------
import myutil  # basic utility
import nu  # numerical unit
# defining photon flux, concentratoin factor etc
import photon_in
# current and energy flux through ESC
import extraction_via_ESC as viaESC
# radiative recombination by generalized planck law
import recombination as rec
# defining material properties, thermalization time, ESC energy level etc
from hcscAttribute import hcscAttribute
# ------------------------------------------------------------


class hcsc(hcscAttribute):
    """ A class used to calculate hot carrier solar cell characteristics

    Attributes
    ----------
    f1 : function
        used finding carrier chemical potential in absorber and carrier temperature. swapable function
    f2 : function
    f3 : function
    success : bool
        variables used to verify convergence

    Methods
    -------
    shine(self, ph_in)
        do this first so that photons are absorbed in the cell
        set light source defined by photon_in class
        Calculate photon flux, energy flux, etc
    photon_flux(self, ph_in)
        photon flux reaching the surface at the absorber
        Quantum efficiency = 0 for E < Eg, 1 for E > Eg
    energy_flux(self, ph_in)
        energy flux
        QE = 0 for E < Eg, 1 for E > Eg
    Frec(self, mu)
        recombination flux at given chemical potential at a contact
        measured from the center of the bandgap. Twice of mu divided by q (charge) equals to applied bias
        integrating over from Eg to infty (Eg+15xkTc)
    Jrec(self, mu)
        recombination current at absorber mu
    Urec(self, mu)
    Jout(self, mu)
        current density output at given mu
        to be consistent, electron current is negative
    Jsc(self)
        short circuit current
    Uout(self, mu)
        extracted energy flux at bias mu, Uabs-Uem(mu)
    QWDBRTSf1(self, muc, Tc, mue)
        These functions are used to find muc and Tc
        define two equations to solve
        solve for muc
        since use optimize muc in unit eV
    QWDBRTSf2(self, Tc, muc, mue)
    QWDBRTSf3(self, Tc, muc, mue)
        solve for Tc in particle conservation equation
        similar to f1 but solve for Tc
    N0(self, _Tph)
    solveTph(self, muc, Tc)
    UthPOP(self, muc, Tph=300)
        energy loss rate by optical phonon
    Uph(self, Tph=300)
        optical phonon decay into acoustic modes
    Nph(self, Tph=300)
        phonon population given by Bose Einstein distribution
    solve_mu_Tc(self, mue)
        solve absorber carrier chemical potential mu and carrier temperature T at given mue
        solve particle and energy equations simulataneously
        Once mu and T are determined, Jout -> Pout -> Efficiency
        Jout = Jtunnel(mu, T) = q(Fabs-Fem)
        Uout = Utunnel = Uabs-Eem
        set reference potential 0 at mid bandgap
        symmetric effective mass two-band model
    Jouthc(self, mue)
        Current output of the hot carrier
        return
    Vochc(self, fmaxP=False, dsp_msg=True)
        solve Voc
        mue should be in unit eV so that brentq gives accurate value
    Pouthc(self, mue)
        return poewr output at given mue
        note mue is measured from Eg/2
        symmetric two band model
        hence mue difference at contacts is 2*mue
    maxPouthc(self)
        solve max power at given Eg, ESC E and Ew
    """

    def __init__(self):
        hcscAttribute.__init__(self)
        self.f1 = self.QWDBRTSf1  # these are used in muc and Tc
        self.f2 = self.QWDBRTSf2
        self.f3 = self.QWDBRTSf3
        # variables used to verify convergence
        self.success = True
        return

    def shine(self, ph_in):
        """Shine light source defined in photon_in class
        Calculate photon flux, energy flux, etc
        """
        self.Fabs = self.photon_flux(ph_in)  # #/sm^2
        self.Jabs = nu.q*self.Fabs  # A/sm^2
        self.Uabs = self.energy_flux(ph_in)  # W/m^2
        return

    def photon_flux(self, ph_in):
        """photon flux into the absorber
        Quantum efficiency (QE) = 0 for E < Eg, 1 for E > Eg

        Parameters
        ----------
        ph_in : class Photon_in defined in photon_in.py
            defining spectrum (AM15), concentratoin factor and it does
            integration
        """
        return ph_in.photon_flux(self.absb.Eg)

    def energy_flux(self, ph_in):
        """energy flux into the absorber
        QE = 0 for E < Eg, 1 for E > Eg
        """
        return ph_in.photon_intensity(self.absb.Eg)

    def Frec(self, mu):
        """recombination flux at given chemical potential
        integrating over dhbarw from Eg to infty
        """
        infty = self.absb.Eg+15*sc.k*self.absb.T
        ret = integrate.quad(rec.frec, self.absb.Eg,
                                infty, args=(mu, self.absb.T))
        # print 'Frec', ret, mu/nu.eV, self.absb.T
        return ret[0]

    def Jrec(self, mu):
        """recombination current at absorber

        Args:
            mu: chemical potential
        """
        return nu.q*self.Frec(mu)

    def Urec(self, mu):
        """Energy flux by radiative recombination at chemical potential mu at absorber
        """
        infty = self.absb.Eg+15*sc.k*self.absb.T
        ret = integrate.quad(rec.urec, self.absb.Eg,
                                infty, args=(mu, self.absb.T))
        return ret[0]

    def Jout(self, mu):
        """current density output at given mu
        to be consistent, electron current is negative
        """
        ret = self.Fabs-self.Frec(mu)
        return -nu.eV*ret

    def Jsc(self):
        """short circuit current
        """
        return self.Jout(0)

    def Uout(self, mu):
        """extracted energy
        """
        ret = self.Uabs-self.Urec(mu)
        return ret

    def QWDBRTSf1(self, muc, Tc, mue):
        """
        define two equations to solve
        solve for muc
        since use optimize muc in unit eV
        """
        self.absb.T = Tc
        try:
            ret1 = self.Jout(muc*nu.eV)
            ret2 = self.Jext(muc*nu.eV, mue)
            ret = -(ret1-ret2)
        except FloatingPointError as e:
            print(e)
            ret = -1
        return ret

    def QWDBRTSf2(self, Tc, muc, mue):
        self.absb.T = Tc
        try:
            ret1 = self.Uabs-self.Urec(muc)-2*self.Uext(muc, mue)
            ret2 = self.Uth(muc)
            ret = ret1-ret2
        except FloatingPointError as e:
            print(e)
            ret = -1
        return ret

    def QWDBRTSf3(self, Tc, muc, mue):
        """
        solve for Tc in particle conservation equation
        similar to f1 but solve for Tc
        """
        self.absb.T = Tc
        try:
            ret1 = self.Jout(muc)
            ret2 = self.Jext(muc, mue)
            ret = -(ret1-ret2)
        except FloatingPointError as e:
            print(e)
            # print traceback.format_exc()
            # print 'f3', muc/sc.e, self.absb.T, mue/sc.e
            ret = -1
        return ret

    def N0(self, _Tph):
        return 1/(np.exp(self.hw0/(sc.k*_Tph))-1)

    def solveTph(self, muc, Tc):
        self.Neq = self.N0(self.Trm)
        self.absb.T = Tc

        def f(Tph, muc):
            ret1 = self.UthPOP(muc, Tph)
            ret2 = self.Uph(Tph)
            return ret1-ret2
        ret, r = optimize.brentq(f, 300, 1400, args=(muc,), full_output=True)
        if r.converged is not True:
            return self.Trm, False
        return ret, True

    def UthPOP(self, muc, Tph=300):
        """
        energy loss by optical phonon
        """
        ret1 = self.absb.density(muc)
        N0 = self.Nph(Tph)
        ret2 = self.Ephn*((N0+1)*np.exp(-self.Ephn/(sc.k*self.absb.T))-N0)
        # factor of 2 for electrons and holes
        ret = 2*ret2*ret1*self.d/self.tau_th
        return ret

    def Uph(self, Tph=300):
        """
        optical phonon decay into acoustic modes
        """
        N0 = self.Nph(Tph)
        Neq = self.Nph(self.Trm)
        return self.hw0*self.NM/self.tau_pp*(N0-Neq)

    def Nph(self, Tph=300):
        """
        phonon population given by Bose Einstein distribution
        """
        return 1/(np.exp(self.hw0/(sc.k*Tph))-1)

    def solve_mu_Tc(self, mue):
        """
        solve absorber mu (chemical potential at absorber) and T (carrier temperature at absorber) at given mu_e (chemical potential at contact)

        solve current and energy equations simulataneously
        3 unknown variables mu_c, Tc, and mu_e to solve
        since mu_e is assumed, mu_c and Tc can be determined
        which leads to Jout -> Pout -> Efficiency
        Jout = Jtunnel(mu, T) = q(Fabs-Fem)
        Uout = Utunnel = Uabs-Eem
        set potential 0 at Eg/2
        symmetric effective mass two-band model
        solve for Tc
        """
        f1 = self.f1
        f2 = self.f2
        f3 = self.f3

        def narrowmu(mu1, mu2, Tc):
            # the curve made by f1 function should be positive to negative
            mu = np.linspace(mu1, mu2)
            ret = np.zeros(mu.size)
            for cnt, _ in enumerate(mu):
                try:
                    ret[cnt] = f1(_/nu.eV, Tc)
                    if (ret[cnt] < 0):
                        break
                except FloatingPointError as e:
                    print(e)
                    print(traceback.format_exc())
                    print('narrowmu', muc/sc.e, self.absb.T, mue/sc.e)
                    return 0, 0, False
            ret1 = mu[(ret > 0)]
            ret2 = mu[(ret < 0)]
            return ret1[-1], ret2[0], True

        def solvemuc(Tc):
            # to keep exp from blows up
            # (E-2mu)/(kT) has to be less than 300 or so
            minmu = self.absb.Eg/2-100*sc.k*Tc
            # test_plot_mu(minmu, 5*nu.eV, Tc)
            muc1, muc2, success = narrowmu(minmu, 5*nu.eV, Tc)
            if success is False:
                return -1, False
            # print 'guess', muc1/nu.eV, muc2/nu.eV
            ret, r = optimize.brentq(f1, muc1/nu.eV, muc2/nu.eV, args=(Tc, ),
                                     full_output=True)
            # print 'result muc:{:g} (eV)'.format(muc/nu.eV)
            if r.converged is False:
                print('Convergence Failure!!')
                return -1, False
            return ret*nu.eV, True

        def narrowT(T1, T2, muc, func=f2, fcnt=0):
            # the curve made by f2 function should be positive to negative
            # the curve made by f3 function should be positive to negative
            T = np.linspace(T1, T2)
            ret = np.ones(T.size)*-1
            for cnt, _ in enumerate(T):
                try:
                    ret[cnt] = func(_, muc, mue)
                    if (ret[cnt] < 0):
                        break
                except FloatingPointError as e:
                    print(e)
                    print(traceback.format_exc())
                    print('func in narrowT', muc/sc.e, self.absb.T, mue/sc.e)
                    return 0, 0, False
            # print T1, T2, muc/nu.eV, cnt, func
            # plt.plot(T[:cnt+1], ret[:cnt+1])
            # # # plt.ylim(-50, max(ret[:cnt]))
            # plt.show()
            # print "ret", ret[:cnt+1]
            # print "T", T[ret > 0]
            if cnt > 0:
                # ret1 = T[(ret[:cnt+1] > 0)]
                # ret2 = T[(ret[:cnt+1] < 0)]
                ret1 = T[(ret > 0)]
                ret2 = T[(ret < 0)]
                # print 'narrowT', muc/nu.eV, ret[:cnt], ret1, ret2, self.rcnt.T
                # rcnt.T ... right contact temperature
                if (ret1[-1] < self.rcnt.T) or (ret2[0] < self.rcnt.T):
                    return 0, 0, False
                # if there are only less than 3 points but more than 2 points
                # are taken, go for a detail to have a better solution
                elif cnt < 5:
                    # 1st time try
                    if fcnt < 1:
                        retT1, retT2, success = narrowT(ret1[-1], ret2[0],
                                                        muc, func, fcnt+1)
            else:
                # no solution exists
                # print 'no solution in narrowT'
                return 0, 0, True
            return ret1[-1], ret2[0], True

        def solveT(muc, minT=300, maxT=7200, func=f2):
            T1, T2, success = narrowT(minT, maxT, muc, func)
            # no solution, numerical error is detected
            if success is False:
                return 0, False
            if (T1 == 0) or (T2 == 0):
                return 0, True
            # print 'guess', T1, T2
            ret, r = optimize.brentq(func, T1, T2, args=(muc, mue),
                                     full_output=True)
            # print 'result Tc:{0:4.15f}(K), r:{1:}'.format(ret, r.converged)
            if r.converged is False:
                print('Convergence Failure!!')
                return 0, False
            return ret, True

        def solvef2T(muc, minT=300, maxT=7200):
            return solveT(muc, minT, maxT, f2)

        def solvef3T(muc, minT=300, maxT=7200):
            return solveT(muc, minT, maxT, f3)

        def test_plot_T(T1, T2, muc, func=f2):
            T = np.linspace(T1, T2)
            ret = np.zeros(T.size)
            for cnt, _ in enumerate(T):
                ret[cnt] = func(_, muc, mue)
                if ret[cnt] < 0:
                    break
            # print ret
            print(ret)
            plt.plot(T[:cnt+1], ret[:cnt+1])
            # plt.ylim(-50, max(ret[:cnt]))
            plt.xlabel(r'Carrier temperature $T_c$ (K)')
            plt.ylabel(r'Value of f$_p$')
            plt.show()
            return

        def test_plot_mu(mu1, mu2, Tc):
            mu = np.linspace(mu1, mu2, 100)
            ret = np.zeros(mu.size)
            for cnt, _ in enumerate(mu):
                ret[cnt] = f1(_/nu.eV, Tc)
                if ret[cnt] > 0:
                    cnt += 1
                    break
            plt.plot(mu[:cnt]/nu.eV, ret[:cnt])
            plt.ylim(min(ret[:cnt]), 50)
            plt.show()
            return

        def test_plot(minmu, maxmu, minT=300, maxT=7200, fcnt=0):
            """
            the function giving the best guess to solve two equations
            if g_msg is True:
                print 'test_plot', minmu/nu.eV, maxmu/nu.eV, minT, maxT
            maxmuc has to be less than 2*maxmu < Eg
            to avoid divergent behavior in integral in f2 and f3
            1/exp((E-2*muc)/(kT)-1) from Eg to Eg + 15kT
            """
            if maxmu > self.absb.Eg*.49:
                maxmu = self.absb.Eg*.49
            # 100 and 50 gave similar result
            if g_msg is True:
                print('minmu, maxmu', minmu/nu.eV, maxmu/nu.eV)
            mu = np.linspace(minmu, maxmu, 25*(fcnt+1))
            Tf3 = np.zeros(mu.size)
            Tf2 = np.zeros(mu.size)
            # d = np.zeros(mu.size)
            sign = 1
            success = True
            maxrepeat = 3
            for cnt, _mu in enumerate(mu):
                # test_plot_T(minT, 3000, _mu, f2)
                # test_plot_T(minT, 3000, _mu, f3)
                # exit()
                Tf2[cnt], success1 = solvef2T(_mu, minT, maxT)
                Tf3[cnt], success2 = solvef3T(_mu, minT, maxT)
                if g_msg is True:
                    print('test_plot:result', Tf2[cnt], Tf3[cnt], _mu/nu.eV, cnt)
                # any error detected
                if success1 is False or success2 is False:
                    success = False
                    print('error')
                    break
                if Tf2[cnt] > 0 and Tf3[cnt] > 0:
                    # there is a solution if sign changes
                    sign = (Tf2[0]-Tf3[0])*(Tf2[cnt]-Tf3[cnt])
                    # print sign, Tf2[0]-Tf3[0], Tf2[cnt]-Tf3[cnt]
                else:
                    sign = 1
                # 0 return -> maybe solution close to 300 K
                if (sign < 0) or (Tf2[cnt] == 0) or (Tf3[cnt] == 0):
                    # try if solver can find a solution with T at cnt-1
                    # if sign changes, there is a solution
                    # if not, maybe need more detail
                    if g_msg is True:
                        print(('test_plot:Tf2 or Tf3 = 0 or sign < 0',
                               Tf2[cnt], Tf3[cnt], sign))
                        print(Tf2[cnt-1], Tf3[cnt-1])
                    if sign < 0:
                        if fcnt > 0:
                            fcnt += maxrepeat
                            break
                    if cnt < 1:
                        # print 'no solution'
                        break
                    # Tf2[cnt] = Tf2[cnt-1]
                    # Tf3[cnt] = Tf3[cnt-1]
                    if fcnt > maxrepeat:
                        # give up narrowing down
                        # return closest values so far
                        # which is cnt-1
                        break
                    minT = 300
                    maxT1 = max([Tf2[cnt], Tf3[cnt]])
                    maxT2 = max([Tf2[cnt-1], Tf3[cnt-1]])
                    maxT = max([maxT1, maxT2])*1.01
                    # test_plot_T(300, maxT, mu[cnt-1], f2)
                    # test_plot_T(300, maxT, mu[cnt-1], f3)
                    muc, Tc, success, sign = test_plot(mu[cnt-1],
                                                       mu[cnt], minT,
                                                       maxT, fcnt+1)
                    break
            # __mu = 0.0105653101669*nu.eV
            # test_plot_T(300, 7200, __mu, f2)
            # test_plot_T(300, 7200, __mu, f3)
            # double check. if no solution, return -1 K
            if success is False:
                mu[cnt], -1
            if fcnt > maxrepeat:
                if g_msg is True:
                    print(('test_plot:considered values',
                           Tf2[cnt], Tf3[cnt], Tf2[cnt-1], Tf3[cnt-1]))
                    print(mu[cnt]/nu.eV, mu[cnt-1]/nu.eV)
                muc = (mu[cnt]+mu[cnt-1])/2
                Tc = ((Tf2[cnt]+Tf2[cnt-1])/2+(Tf3[cnt]+Tf3[cnt-1])/2)/2
                # Tc = (Tf2[cnt]+Tf3[cnt])/2
            if g_msg_plt is True:
                print(fcnt, cnt, sign, mu[cnt]/nu.eV, Tf2[cnt], Tf3[cnt])
                if fcnt > 0:
                    print(fcnt, muc/nu.eV, Tc)
                plt.scatter(mu[:cnt+1]/nu.eV, Tf2[:cnt+1])
                plt.scatter(mu[:cnt+1]/nu.eV, Tf3[:cnt+1])
                plt.plot(mu[:cnt+1]/nu.eV, Tf2[:cnt+1])
                plt.plot(mu[:cnt+1]/nu.eV, Tf3[:cnt+1])
                plt.legend([r'f$_E$', r'f$_p$'])
                plt.xlabel(r'$\mu_c$(eV)')
                plt.ylabel(r'T$_c$(K)')
                plt.title(r'$\mu_e$(eV):{0:2.2f}'.format(mue/nu.eV))
                plt.show()
            if sign > 0:
                print('no solution', fcnt)
                return mu[cnt], -1, False, 0
            else:
                success = True
            return muc, Tc, success, sign

        def equations(p, mue):
            muc, Tc = p
            ret = (f1(muc/nu.eV, Tc, mue), f2(Tc, muc, mue))
            # print p[0]/nu.eV, p[1], ret
            return ret

        # let's solve for Tc first assuming muc=0
        if g_msg is True:
            print('solve muc Tc mue:{0:2.2f}'.format(mue/nu.eV))
        # flag for routine convergence
        success = True
        if self.EP is True:
            maxcnt = 5
            # set initial phonon temperature
            self.Tph = 300
        else:
            maxcnt = 1
        cnt = 0
        while (cnt < maxcnt):
            cnt += 1
            # estimate initial guess for muc and Tc
            # minmu = self.absb.Eg/2-140*sc.k*300
            minmu = self.absb.Eg/2-125*sc.k*300
            muc, Tc, success, sign = test_plot(minmu, 1*nu.eV)
            if Tc < 0:
                success = False
                return -1, -1, success
            if self.EP is True:
                ret, r = self.solveTph(muc, Tc)
                if r is True:
                    self.Tph = ret
                else:
                    print('no convergence Tph')
                if g_msg is True:
                    print("solveTph: muc:{0:f} Tc:{1:f}".format(muc/nu.eV, Tc))
                    print('Tph:{:f}'.format(self.Tph))
            if g_msg is True:
                print('guess muc:{0:g}(eV) Tc:{1:g}(K)'.format(muc/nu.eV, Tc))
        return muc, Tc, True

    def Jouthc(self, mue):
        """Extracted current at given mue
        return Jout(electron) at given mue (chemical potential at a contact)

        Returns
        -------
        float
            Current J
        bool
            If False, no solution was found
        """
        muc, Tc, success = self.solve_mu_Tc(mue)
        if success is True:
            self.absb.T = Tc
            J = self.Jout(muc)
        else:
            self.absb.T = self.rcnt.T
            # positive electron current -> negative current
            J = 1
        # electron current is negative
        # but to photogenerated current positive multiply by -1
        return -J, success

    def Vochc(self, fmaxP=False, dsp_msg=True):
        """solve open circuit voltage (Voc) of hot carrier solar cell
        mue should be in unit eV so that brentq gives accurate value

        if fmaxP is True skip scanning at lower bias like less than 30%
        for speed up process
        But cannot be used to display IV or PV curve

        Returns
        -------
        float
            Open circuit voltage Voc
        list
            chemical potential at a contact mue (J)
        list
            current flux J (A/m^2)
        list
            Power flux P (W/m^2)
        list
            Carrier temperature T (K)
        """
        def f(mue):
            ret, success = self.Jouthc(mue*nu.eV)
            # print ret
            return ret, success
        # if mue > resc.E then no current can flow
        maxmue = self.resc.E
        if fmaxP is True:
            if self.absb.Eg < .5*nu.eV:
                minmue = maxmue*.2
            elif self.absb.Eg < 1.0*nu.eV:
                minmue = maxmue*.3
            elif self.absb.Eg < 1.5*nu.eV:
                minmue = maxmue*.4
            elif self.absb.Eg < 1.7*nu.eV:
                minmue = maxmue*.5
            elif self.absb.Eg < 1.8*nu.eV:
                minmue = maxmue*.6
            elif self.absb.Eg < 2.1*nu.eV:
                minmue = maxmue*.7
        else:
            minmue = 0
        # use equal interval at every bandgap Eg
        # mue = np.linspace(minmue, maxmue, 7)
        mue = np.arange(minmue, maxmue, 2*nu.meV)
        if fmaxP is True:
            # test if there is a solution at short circuit condition
            mue = np.insert(mue, 0, 0)
        J = np.zeros(mue.size)
        P = np.zeros(mue.size)
        T = np.zeros(mue.size)
        if dsp_msg is True:
            print("Solving Voc")
        for cnt, _ in enumerate(mue):
            J[cnt], success = self.Jouthc(_)
            T[cnt] = self.absb.T
            P[cnt] = J[cnt]*2*_/nu.eV
            if dsp_msg is True:
                print('Count:{:d} Tc:{:.2f} mue(eV):{:.3f} Bias(V):{:.3f} J(A/m^2):{:.2f} P(W/m^2):{:.2f} '.format(cnt, self.absb.T, _/nu.eV, _/nu.eV*2, J[cnt], P[cnt]))
            if J[cnt] < 0:
                print('Stop iteration at mue(eV):{0:e}'.format(_/nu.eV))
                break
            if fmaxP is True:
                if cnt < 0:
                    pass
                elif P[cnt] < P[cnt-1]:
                    print(('fmaxP:Vochc Stop iteration at' +
                           'mue:{0:e} (eV)'.format(_/nu.eV)))
                    break
        # J is positive, then declines
        # print 'Vochc', J
        # optimization failure
        if success is False:
            print('Optimization failure')
            return 0, mue/nu.eV, np.zeros(mue.size), np.zeros(mue.size), np.zeros(mue.size)
        # if cnt reaches whole values, then no Voc
        if cnt == (mue.size-1):
            print('Voc is over the limit', cnt)
            return 0, mue/nu.eV, np.zeros(mue.size), np.zeros(mue.size), np.zeros(mue.size)
        elif cnt == 0:
            print('No Voc', cnt)
            return 0, np.zeros(mue.size), np.zeros(mue.size), np.zeros(mue.size), np.zeros(mue.size)
        ret = np.argmax(J < 0)
        ret1 = mue[ret-1]
        ret2 = mue[ret]
        ret = ret1
        return ret, mue[:cnt+1], J[:cnt+1], P[:cnt+1], T[:cnt+1]

    def Pouthc(self, mue):
        """return poewr output at given mue (a half is bias at contact in eV)
        Power = J x V (delta mu/eV)
        mue is measured from mid bandgap (symmetric two band model)
        delta mu = 2*mue

        Returns
        -------
        float
            power P (W/m^2)
        bool
            if False, solution was not found
        """
        J, success = self.Jouthc(mue)
        if success is False:
            print('Optimization failure')
            return -1
        P = J*2*mue/nu.eV
        return P, success

    def maxPouthc(self):
        """solve max power at given condition defined in hcsc

        Returns
        -------
        float
            maximum power
        float
            maximum carrier temperature
        """
        # open circuit voltage
        voc, mue, J, P, T = self.Vochc()
        maxP = max(P)
        # carrier temperature at maximum power point
        Tmpp = T[np.argmax(P)]
        # print "maxP: {0:g}".format(maxP)
        return maxP, Tmpp


def fmax(func_to_maximize, initial_guess=0):
    """return the x that maximizes func_to_maximize(x)
    a general function
    """
    func_to_minimize = lambda x: -func_to_maximize(x)
    return optimize.fmin(func_to_minimize, initial_guess, disp=False)[0]


# global variable
# if True, show message and plot during optimization process
g_msg = False
g_msg_plt = False


def init():
    """set up debus, fonts, error for numpy
    """
    # for debugs
    global g_msg
    global g_msg_plt
    global g_skip
    g_msg = False
    g_msg_plt = False
    # setup bigger fonts for plotting
    myutil.setup_fonts()
    # redirect_output_to_file()
    np.seterr(over='raise', divide='raise')


def main():
    init()
    scell = hcsc()
    scell.absb.Eg = 0.65*nu.eV
    ph_in = photon_in.Photon_in()
    ph_in.c = 10
    scell.shine(ph_in)
    scell.display_attributes()
    print("Jout:{:.3f}(A/m^2)".format(scell.Jouthc(0.2/2*nu.eV)[0]))
    print("Pout:{:.3f}(W/m^2)".format(scell.Pouthc(0.2/2*nu.eV)[0]))
    print("MaxPout:{:.3f}(W/m^2)".format(scell.maxPouthc()[0]))

    if scell.success is False:
        print('Something went wrong')
    return

if __name__ == '__main__':
    main()
