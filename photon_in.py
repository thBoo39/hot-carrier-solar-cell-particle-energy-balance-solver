import numpy as np
from scipy import constants as sc
from scipy import interpolate
from scipy import integrate
from scipy import optimize
import time
import os
import multiprocessing as mp
from matplotlib import pyplot as plt
# custom import
import nu


def NREL_AM15_spectra(n=2):
    """ NREL AM15 Spectra W/m^2/nm
    n=2 for AM15 spectra.
    """
    target_file = 'ASTMG173.csv'
    NREL_AM15_file = os.path.join(os.getcwd(), target_file)
    downloaded_array = np.genfromtxt(NREL_AM15_file,
                                     delimiter=',', skip_header=2)
    AM15 = downloaded_array[:, [0, n]]
    AM15[:, 0] *= 1e-9  # nm
    AM15[:, 1] *= nu.W * nu.m**-2 * nu.nm**-1
    return AM15


def radiative_recombination(V, Eg, T_cell=300):
    E_max = Eg+10*sc.k*T_cell
    if V*sc.e > E_max:
        return 0
    g = sc.e*2*np.pi/(sc.c**2*sc.h**3)
    fN = lambda E: E**2/(np.exp((E-sc.e*V)/(sc.k*T_cell))-1)
    # E = np.linspace(Eg, E_max)
    # plt.plot(E/sc.e, fN(E))
    # plt.show()
    return (g*integrate.quad(fN, Eg, E_max)[0])


def current_density(V, Eg):
    N_sun = photon_flux_from_NREL_data(Eg)
    N_env = photon_flux_from_env(Eg)
    N_rec = radiative_recombination(V, Eg)
    # print (1-F*C)*N_env, N_rec
    return sc.e*(C*N_sun+(1-F*C)*N_env-N_rec)


def V_mpp(Eg):
    f = lambda V: -V*current_density(V, Eg)
    return optimize.fmin(f, 0)[0]


def J_mpp(Eg):
    return current_density(V_mpp(Eg), Eg)


def P_max(Eg):
    V = V_mpp(Eg)
    return V*current_density(V, Eg)


class Photon_in(object):
    """ defining spectrum
    default spectrum is AM15 if n=2

    attributes
    ----------
    E_flux_raw : numpy array
    c : float
        concentratoin factor
    T_env : float
    f : float
    """
    def __init__(self,
                 n=2,
                 c_f=1,
                 E_flux_func=NREL_AM15_spectra,
                 T_env=300.0):
        self.E_flux_raw = E_flux_func(n)
        self.T_env = T_env
        self.c = c_f  # concentration factor
        self.f = 2.153e-5
        self.reset()

    def reset(self):
        self.num_data = np.size(self.E_flux_raw)
        self.wavelength_min = self.E_flux_raw[0, 0]
        self.wavelength_max = self.E_flux_raw[-1, 0]
        self.E_min = sc.h*sc.c/self.wavelength_max
        self.E_max = sc.h*sc.c/self.wavelength_min
        wav = self.E_flux_raw[:, 0]
        E_flux = self.E_flux_raw[:, 1]
        self.E_flux = interpolate.interp1d(wav, E_flux)
        # dlam = 1
        # dE = sc.h*sc.c/wav**2*dlam
        # E = sc.h*sc.c/wav
        # E_flux = self.E_flux_wv(sc.h*sc.c/E)*1
        # self.E_flux = interpolate.interp1d(E[::-1], E_flux[::-1])
        self.intensity = self.photon_intensity()  # intensity W/cm^-2/s
        self.flux = self.photon_flux()  # photons/cm^-2/s
        return

    def get_E_min(self):
        return self.E_min

    def photon_intensity(self, Eg=0):
        if Eg < self.E_min:
            Eg = self.E_min
        # return integrate.quad(self.E_flux, Eg, self.E_max)[0]*nu.cm**2
        fN = lambda E: self.E_flux(sc.h*sc.c/E)*sc.h*sc.c/E**2
        # return integrate.quad(fN, Eg, self.E_max)[0]*nu.cm**2
        return self.c*integrate.quad(fN, Eg, self.E_max)[0]

    def photon_fluxperE(self, E):
        # return self.E_flux(sc.h*sc.c/E)*sc.h*sc.c/E**3*nu.cm**2
        return self.c*self.E_flux(sc.h*sc.c/E)*sc.h*sc.c/E**3

    def photon_flux(self, Eg=0):
        if Eg < self.E_min:
            Eg = self.E_min
        fN = lambda E: self.E_flux(sc.h*sc.c/E)*sc.h*sc.c/E**3
        # fN = lambda E: self.E_flux(E)/E
#         return integrate.quad(fN, Eg, self.E_max)[0]*nu.cm**2
        ret = self.c*integrate.quad(fN, Eg, self.E_max)[0]
        if ret < 0:
            ret = float('nan')
        return ret

    def photon_flux_from_env(self, Eg=0):
        if Eg < self.E_min:
            Eg = self.E_min
        g = 2*np.pi/(sc.c**2*sc.h**3)
        fN = lambda E: E**2/(np.exp((E)/(sc.k*self.T_env)))
#         return (g*integrate.quad(fN, Eg, self.E_max)[0])*nu.cm**2
        return self.c*(g*integrate.quad(fN, Eg, self.E_max)[0])

    def plot_spectra(self):
        wv = np.linspace(self.wavelength_min, self.wavelength_max, 200)
        plt.plot(wv/nu.nm, self.E_flux(wv)/(nu.W/nu.m**2/nu.nm))
        plt.xlabel('wavelength (nm)')
        plt.ylabel(r'intensity (W/m$^2$/nm)')
        plt.show()
        return


class laser_diode(Photon_in):
    def __init__(self, P=1, wv_mu=781*nu.nm, wv_sig=10*nu.nm, T_env=300):
        self.E_flux_raw = self.spectra(P, wv_mu, wv_sig)
        self.T_env = T_env
        self.c = 1
        self.reset()

    def spectra(self, P=1, wv_mu=781*nu.nm, wv_sig=10*nu.nm):
        # produce Gaussian beam centered at wv_mu
        # W/m^-2/nm
        A = P/(np.sqrt(2*np.pi)*wv_sig)
        wv = np.arange(wv_mu-wv_sig*4, wv_mu+wv_sig*4, nu.nm)
        # dwv = wv[1]-wv[0]
        spectra = A*np.exp(-(wv-wv_mu)**2/(2*wv_sig**2))
        return np.column_stack((wv, spectra))

    def adjust_P(self, P, wv_mu=781*nu.nm, wv_sig=10*nu.nm):
        self.E_flux_raw = self.spectra(P, wv_mu, wv_sig)
        self.reset()
        return


# for multiprocessing
# find average energy flux at given energy gap
def worker(_):
    flux = g_ph_in.photon_flux(_)
    intn = g_ph_in.photon_intensity(_)
    E_ave = intn/flux
    return E_ave


# multiprocessing initializer
def initProcess(ph_in):
    global g_ph_in
    g_ph_in = ph_in


# global variable for multiprocessing
g_ph_in = Photon_in()


def main():
    ph_in = Photon_in()
    laser_in = laser_diode(1000*1000)
    E_min = ph_in.E_min
    ph_in.plot_spectra()
    print(laser_in.flux)
    print(radiative_recombination(0, 1.4*sc.e, 1000))
    laser_in.adjust_P(1000)
    laser_in.plot_spectra()
    print(ph_in.photon_fluxperE(2*nu.eV)*1e-3*sc.e)
    Eg = np.linspace(E_min, 3*sc.e)
    t = time.time()
    pool = mp.Pool(4, initializer=initProcess, initargs=(ph_in, ))
    E_ave = pool.imap(worker, Eg)
    print("Time taken: %.5g" % (time.time()-t))
    E_ave = np.fromiter(E_ave, np.float)
    print("Results"+'-'*40)
    print("E_ave", E_ave/sc.e)
    plt.plot(Eg/sc.e, E_ave/sc.e)
    plt.xlabel('Energy gap (eV)')
    plt.ylabel('Average Excess electron energy (eV)')
    plt.show()

if __name__ == "__main__":
    main()
