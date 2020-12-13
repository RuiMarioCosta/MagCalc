import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from constants import Constants


class Magnetization(Constants):
    def __init__(self, J, gJ, Tc, Nm):
        self.J = J  # Total Angular Momentum, J
        self.gJ = gJ  # Landé g-factor
        self.Tc = Tc  # Curie Temperature, in Kelvin
        self.Nm = Nm  # Number of Magnetic Moments in Primitive Cell

        # Curie Constants divided by Vacuum Permeability (C/mu_0)
        self.const = (
            (self.mu_B ** 2.0)
            * self.Nm
            * (self.gJ ** 2.0)
            * self.J
            * (self.J + 1.0)
            / (3.0 * self.k_B)
        )
        self.lamb = (
            self.Tc / self.const
        )  # Value of the strength of the parameter of the Molecular Field

    def reduced_magnetization(self, T, B):
        """Calculates the reduced magnetization (Brillouin function).

        Parameters
        ----------
        T : scalar
            Temperature
        B : scalar
            Magnetic fields

        Returns
        --------
        y : scalar
            Reduced magnetization
        """

        # Function for the computation of sigma
        def brillouin(sigma, T, B):
            h = B / (self.lamb * self.Nm * self.gJ * self.mu_B * self.J)
            y = 3.0 * self.J / (self.J + 1.0) * (h + sigma) * self.Tc / T
            return (
                sigma
                - (2.0 * self.J + 1.0)
                / (2.0 * self.J * np.tanh((2.0 * self.J + 1.0) * y / (2.0 * self.J)))
                + 1.0 / (2.0 * self.J * np.tanh(y / (2 * self.J)))
            )

        starting_estimate = 0.5 if B >= 0 else -0.5
        sigma = optimize.fsolve(brillouin, starting_estimate, args=(T, B))

        return sigma[0]

    def magnetic_entropy(self, T, B):
        """Computes the magnetic entropy.

        Parameters
        ----------
        T : scalar
            Temperatures.
        B : scalar
            Magnetic fields.

        Returns
        -------
        y : scalar
            Magnetic entropy.
        """
        # relative field to the saturation magnetization
        h = B / (self.lamb * self.Nm * self.gJ * self.mu_B * self.J)

        sigma = self.reduced_magnetization(T, B)  # reduced magnetization

        y = (
            3.0 * self.J / (self.J + 1.0) * (h + sigma) * self.Tc / T
        )  # temporary variable

        a = np.sinh((2.0 * self.J + 1.0) * y / (2.0 * self.J))
        b = np.sinh(y / (2.0 * self.J))

        return self.k_B * self.Nm * (np.log(a / b) - sigma * y)

    def magnetic_free_energy(self, T, B):
        """Magnetic free energy of 1 spin.

        Parameters
        ----------
        T : scalar
            Temperatures.
        B : scalar
            Magnetic fields.

        Returns
        -------
        y : scalar
            Magnetic free energy.
        """
        h = B / (
            self.lamb * self.Nm * self.gJ * self.mu_B * self.J
        )  # relative field to the saturation magnetization
        sigma = self.reduced_magnetization(T, B)
        y = 3.0 * self.J / (self.J + 1.0) * (h + sigma) * self.Tc / T
        a = np.sinh((2.0 * self.J + 1.0) * y / (2.0 * self.J))
        b = np.sinh(y / (2.0 * self.J))

        return -T * self.Nm * self.k_B * np.log(a / b)

    def magnetic_energy(self, T, B):
        """Computes the magnetic energy.

        Parameters
        ---------
        T : scalar
            Temperatures.
        B : scalar
            Magnetic fields.

        Returns
        -------
        y : scalar
            Magnetic energy.
        """
        sigma = self.reduced_magnetization(T, B)

        return -self.gJ * self.mu_B * self.Nm * self.J * B * sigma - 3.0 * self.J / (
            self.J + 1.0
        ) * self.k_B * self.Nm * self.Tc * (sigma ** 2.0)

    def plot_reduced_magnetization(self, T, B):
        """Plots the reduced magnetization as a function of temperature or magnetic field.

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
        magnetization = []
        if T.ndim == 1 and B.ndim == 0:
            for t in T:
                magnetization.append(self.reduced_magnetization(t, B))
            plt.plot(T, magnetization, label="B=" + str(B) + "T")
            plt.xlabel("T (K)")
        elif B.ndim == 1 and T.ndim == 0:
            for b in B:
                magnetization.append(self.reduced_magnetization(T, b))
            plt.plot(B, magnetization, label="T=" + str(T) + "K")
            plt.xlabel("B (T)")
        else:
            raise ValueError("Temperature and Magnetic Field need to be a combination of array and scalar")

        plt.title("Reduced Magnetization, $\\sigma$")
        plt.legend(loc=0, fontsize="small")
        plt.ylim(-1.05, 1.05)
        plt.ylabel("$\\sigma(T, B)$")
        plt.show()

    def plot_magnetic_free_energy(self, T, B):
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

        plt.title("Magnetic Free Energy, $F$")
        plt.legend(loc=0, fontsize="small")
        plt.ylabel("$F(T, B)$")
        plt.show()


if __name__ == "__main__":
    mag = Magnetization(1, 1, 123, 1)
    # print(mag.reduced_magnetization(1, 0))
    # print(mag.reduced_magnetization(1, 1))
    # print(mag.reduced_magnetization(1, -1))
    # print(mag.reduced_magnetization(100000, 0))
    # print(mag.reduced_magnetization(100000, 100000))
    # print(mag.reduced_magnetization(300, 0))
    # print(mag.reduced_magnetization(300, 100000))
    # print()
    # print(mag.magnetic_entropy(1, 0))
    # print(mag.magnetic_entropy(1, 1))
    # print(mag.magnetic_entropy(1, -1))
    # print(mag.magnetic_entropy(100000, 0))
    # print(mag.magnetic_entropy(100000, 100000))
    # print(mag.magnetic_entropy(300, 0))
    # print(mag.magnetic_entropy(300, 100000))
    # print()
    # print(mag.magnetic_free_energy(1, 0))
    # print(mag.magnetic_free_energy(1, 1))
    # print(mag.magnetic_free_energy(1, -1))
    # print(mag.magnetic_free_energy(100000, 0))
    # print(mag.magnetic_free_energy(100000, 100000))
    # print(mag.magnetic_free_energy(300, 0))
    # print(mag.magnetic_free_energy(300, 100000))
    print()
    # temperature = np.arange(1, 200, 0.1)
    # mag.plot_reduced_magnetization(T=temperature, B=0)
    # magnetic_field = np.arange(0, 200, 0.1)
    # mag.plot_reduced_magnetization(T=300, B=magnetic_field)
    print()
    temperature = np.arange(1, 200, 1)
    mag.plot_magnetic_free_energy(T=temperature, B=0)
    magnetic_field = np.arange(-200, 200, 1)
    mag.plot_magnetic_free_energy(T=300, B=magnetic_field)
