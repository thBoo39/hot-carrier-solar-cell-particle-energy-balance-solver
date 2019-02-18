"""test hcsc class function

test module to show IV curve, PV curve, Radiative recombination in terms of carrier temperature and a process to optimize energy extraction level at energy selective contact
"""
import numpy as np
from scipy import constants as sc
from matplotlib import pyplot as plt
import time
# custom imports
import hcsc
import photon_in
import nu


def test_Urec(scell):
    Tc = np.linspace(300, 1000)
    Jrec = np.zeros(Tc.size)
    for count, _ in enumerate(Tc):
        scell.absb.T = _
        Jrec[count] = scell.Urec(0)
    plt.plot(Tc, Jrec)
    plt.xlabel("Carrier temperature (K)")
    plt.ylabel("U by radiative rec (W/m^2)")
    plt.show()
    return


def IV(scell):
    """Plot IV characteristics

    parameters
    ----------
    scell : class hcsc

    Returns
    -------
    none
    """
    t = time.time()
    voc, mue, J, P, T = scell.Vochc()
    print("Time taken:", time.time()-t)
    print('Voc:{0:2.2f}'.format(voc/nu.eV))
    maxP = max(P)
    maxJ = max(J)
    print('maxP (W/m^2):{0:3.1f}'.format(maxP))
    print('maxJ (A/m^2):{0:3.1f}'.format(maxJ))
    plt.plot(2*mue/nu.eV, J)
    plt.ylim(0, max(J)*1.2)
    plt.xlabel("Bias (V)")
    plt.ylabel(r"Current density (A/m$^2$)")
    title_label = ('IV curve')
    plt.title(title_label)
    plt.show()
    return


def PV(scell):
    """Plot PV characteristics
    parameters
    ----------
    scell : class hcsc

    Returns
    -------
    none
    """
    voc, mue, J, P, T = scell.Vochc()
    print('Voc:{0:2.2f}'.format(voc))
    maxP = max(P)
    print("maxP: {0:g}".format(maxP))
    print('Tc:{0:4.1f}'.format(scell.absb.T))
    plt.plot(2*mue/nu.eV, P)
    plt.ylim(0, max(P)*1.2)
    plt.xlabel(r'Bias V (V)')
    plt.ylabel(r'Power output (W/m$^2$)')
    title_label = ('PV Curve')
    plt.title(title_label)
    plt.show()
    return


def maxPouthc_opt_ESC(scell, numEesc=50, g_skip=False):
    """solve maxpower at given Eg optimizing ESC E and Ew
    Energy level of ESC are between 50% of Eg to 300% of Eg

    Parameters
    ----------
    scell : class hcsc
    numEesc : int
        number of trials for ESC energy level
    g_skip : bool
        if g_skip is True, scan all possible values during optimization process possibly more accurate, but significantly increase process time

    Returns
    -------
    maxP : float
        W/m^2
    maxT : float
        K
    optEesc : float
        J
    optEwesc : float
        J
    """
    # this range is about right except for SB and low bandgap case
    minescE = scell.absb.Eg*.5
    maxescE = scell.absb.Eg*3.
    escE = np.linspace(minescE, maxescE, numEesc)
    # define incEesc Eesc increment
    # so that at each Eg, escE has same increment
    # numEesc is referenced to calculate one at Eg 2.0 eV
    # if 100 -> about 50 meV increment
    # 200 -> 25 meV
    # 400 -> 12.5 meV
    incEesc = (3-.5)*2.0*nu.eV/numEesc
    escE = np.arange(minescE, maxescE, incEesc)
    # escEw = np.linspace(.001, 1, 10)*nu.meV
    escEw = np.linspace(5, 30, 1)*nu.meV
    P = np.zeros((escE.size, escEw.size))
    T = np.zeros((escE.size, escEw.size))
    # cntPmax if encounter more than 2 times less than Pmax so far
    # skip iteration
    # there is tendency but don't trust this
    cntPmax = 0
    oldPmax = 0
    # if this flag is set, don't skip optimization
    if g_skip is True:
        skipcnt = 1000
    else:
        # if the calculated efficiency is lower than the max P so far
        # and reach this much number, skip optimization
        skipcnt = int(numEesc/20)
    for cnt1, _E in enumerate(escE):
        scell.resc.E = _E
        for cnt2, _Ew in enumerate(escEw):
            print('cnt1:{0:d}/{1:d} cnt2:{2:d}/{3:d}'.format(cnt1, len(escE), cnt2, len(escEw)))
            print(('deltaEesc:{0:1.2f}(eV) '.format(2*_E/nu.eV) +
                   'Ewesc:{0:3.3e}(meV)'.format(_Ew/nu.meV)))
            scell.resc.Ew = _Ew
            P[cnt1][cnt2], T[cnt1][cnt2] = scell.maxPouthc()
        # this is speed up routine
        # skip routine under the condition given
        # but not 100 % reliable
        if cnt1 > 0:
            # Power was calculated successfully then
            if (P[cnt1][0] < oldPmax) and (P[cnt1][0] > 0):
                # let's say if more than 5 times
                if cntPmax > skipcnt:
                    print('less than before. skip optimization')
                    break
                cntPmax += 1
            elif P[cnt1][0] > oldPmax:
                # new max and renew the counter
                oldPmax = P[cnt1][0]
                cntPmax = 0
    # maximum values for P and T
    # in the range of E and Ew
    maxP = max(P.reshape(P.size))
    maxPidx = np.unravel_index(np.argmax(P), P.shape)
    maxT = T[maxPidx]
    optEesc = escE[maxPidx[0]]
    optEwesc = escEw[maxPidx[1]]
    print("Max Power (W/m^2)", P)
    print("Carrier Temperature (K)", T)
    print("maxP:{0:g} (W/m^2)".format(maxP))
    print("T:{0:g}(K)".format(maxT))
    print(('dEesc:{0:1.2f}(eV) '.format(2*optEesc/nu.eV) +
           'Ewesc:{0:3.3e}(meV)'.format(optEwesc/nu.meV)))
    return maxP, maxT, optEesc, optEwesc


def main():
    hcsc.init()
    scell = hcsc.hcsc()
    ph_in = photon_in.Photon_in()  # AM15
    ph_in.c = 1  # 1 sun concentratoin
    scell.shine(ph_in)  # set light source

    # set parameters, examples
    # effective mass (defined in bulk.py)
    scell.absb.m_e = .1*sc.m_e
    # band gap (defined in bulk.py)
    scell.absb.Eg = 0.65*nu.eV
    # energy selective contact extraction energy
    scell.resc.E = scell.absb.Eg*1.1
    scell.resc.E = 1.35/2*nu.eV  # c 1000 QW
    scell.resc.Ew = 5*nu.meV

    scell.display_attributes()
    # IV(scell)  # show IV curve
    # PV(scell)  # show PV curve
    maxPouthc_opt_ESC(scell, 50)
    # Energy flux by radiative recombination in terms of carrier temperature
    # test_Urec(scell)
    return


if __name__ == '__main__':
    main()
