from constants import Constants


class Electronic(Constants):
    def __init__(self, F0):
        self.F0 = F0  # Internal Energies of Lattice, in eV

    def electronic_free_energy(self):
        return self.F0

    def electronic_energy(self):
        return self.F0