from __future__ import division
import numpy as np
from scipy import constants as sc
# custom module
import nu


class contact(object):
    def __init__(self, T=300, Ef=0*nu.eV):
        # carrier temperature (K)
        self.T = T
        # Fermi energy
        self.Ef = Ef


# modeled as resonant tunneling
class energy_selective_contact(contact):
    def __init__(self, T=300, Ef=0*nu.eV, E=1*nu.eV, Ew=1*nu.meV):
        contact.__init__(self, T, Ef)
        self.E = E  # energy level for tunneling (J)
        self.Ew = Ew  # tunneling energy width (J)
