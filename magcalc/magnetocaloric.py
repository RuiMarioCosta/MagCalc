import numpy as np
import matplotlib.pyplot as plt
from magnetization import Magnetization
from lattice import Lattice
from electronic import Electronic


class Magnetocaloric(Magnetization, Lattice, Electronic):
    def __init__(self, J, gJ, Tc, Nm, theta_D, N, F0, phase='phase1'):
        Magnetization.__init__(self, J, gJ, Tc, Nm)
        Lattice.__init__(self, theta_D, N)
        Electronic.__init__(self, F0)
        self.phase = phase  # Phase of material

    def entropy(self, T, B):
        return self.magnetic_entropy(T, B) + self.lattice_entropy(T)

    def energy(self, T, B):
        return self.magnetic_energy(T, B) + self.lattice_energy(T) + self.electronic_energy()

    def free_energy(self, T, B):
        return self.magnetic_free_energy(T, B) + self.lattice_free_energy(T) + self.electronic_free_energy()

    def plot_free_energy(self, T, B):
        """Plots the free energy as a function of temperature or magnetic field.

        Parameters
        ----------
        T : scalar, array
            Array with the temperatures.
        B : scalar, array
            Magnetic fields.
        """
        plt.figure()

        T = np.asarray(T)
        B = np.asarray(B)
        free_energy = []
        if T.ndim == 1 and B.ndim == 0:
            for t in T:
                free_energy.append(self.magnetic_free_energy(t, B))
            plt.plot(T, free_energy, label="B=" + str(B) + "T")
            plt.xlabel("T (K)")
        elif B.ndim == 1 and T.ndim == 0:
            for b in B:
                free_energy.append(self.magnetic_free_energy(T, b))
            plt.plot(B, free_energy, label="T=" + str(T) + "K")
            plt.xlabel("B (T)")

        else:
            raise ValueError("Temperature and Magnetic Field need to be a combination of array and scalar")

        plt.title("Free Energy, $F$")
        plt.legend(loc=0, fontsize="small")
        plt.ylabel("$F(T, B)$")
        plt.show()


if __name__ == '__main__':
    mc = Magnetocaloric(1, 1, 123, 1, 123, 1, 1)
    print(mc.free_energy(1, 0))
    print(mc.free_energy(1, 1))
    print(mc.free_energy(1, -1))
    print(mc.free_energy(100000, 0))
    print(mc.free_energy(100000, 100000))
    print(mc.free_energy(300, 0))
    print(mc.free_energy(300, 100000))

    temperature = np.arange(1, 200, 1)
    mc.plot_magnetic_free_energy(T=temperature, B=0)
    magnetic_field = np.arange(-200, 200, 1)
    mc.plot_magnetic_free_energy(T=300, B=magnetic_field)
