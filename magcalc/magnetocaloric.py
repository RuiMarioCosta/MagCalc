import numpy as np
from scipy.optimize import fsolve
from magnetization import Magnetization
from lattice import Lattice
from electronic import Electronic


class Magnetocaloric(Magnetization, Lattice, Electronic):
    def __init__(self, phase, J, gJ, Tc, Nm, theta_D, N, F0):
        Magnetization.__init__(self, J, gJ, Tc, Nm)
        Lattice.__init__(self, theta_D, N)
        Electronic.__init__(self, F0)
        self.phase = phase  # Phase of material

    def get(self, arg):
        pass

    def entropy(self):
        pass

    def energy(self):
        pass

    def free_energy(self):
        pass


if __name__ == '__main__':
    mc = Magnetocaloric('phase1', 1, 1, 123, 1, 123, 1, 0)
