import numpy as np
from scipy import integrate
from constants import Constants


class Lattice(Constants):
    def __init__(self, theta_D, N):
        self.theta_D = theta_D  # Debye temperatures, in Kelvin
        self.N = N  # Number of Atoms in Primitive Cell

    def lattice_entropy(self, T):
        """Computes the lattice entropy.

        Parameters
        ----------
        T : scalar
            Temperatures.

        Returns
        -------
        y : Lattice entropy.
        """
        integral, error = integrate.quad(
            lambda x: (x ** 3.0) / (np.exp(x) - 1.0), 0, self.theta_D / T
        )

        return (
            self.k_B
            * self.N
            * (
                12.0 * ((T / self.theta_D) ** 3.0) * integral
                - 3.0 * np.log(1.0 - np.exp(-self.theta_D / T))
            )
        )

    def lattice_free_energy(self, T):
        """Free Lattice Energy according to the Debye Model.

        Parameters
        ----------
        T : scalar
            Temperature.

        Returns
        --------
        y : scalar
            Free Lattice Energy
        """
        integral, error = integrate.quad(
            lambda x: (x ** 3.0) / (np.exp(x) - 1.0), 0, self.theta_D / T
        )

        return (
            self.k_B
            * self.N
            * (
                9.0 / 8.0 * self.theta_D
                - 3.0 * T * ((T / self.theta_D) ** 3.0) * integral
                + 3.0 * T * np.log(1.0 - np.exp(-self.theta_D / T))
            )
        )

    def lattice_energy(self, T):
        """Computes the lattice energy.

        Parameters
        ---------
        T : scalar
            Temperatures.

        Returns
        -------
        y : array
            Lattice energy.
        """
        integral, error = integrate.quad(
            lambda x: (x ** 3.0) / (np.exp(x) - 1.0), 0.0, self.theta_D / T
        )

        return self.k_B * self.N * (
            9.0 / 8.0 * self.theta_D
            + 9.0 * self.theta_D * ((T / self.theta_D) ** 4.0) * integral
        )


if __name__ == "__main__":
    lattice = Lattice(123, 1)
    print(lattice.lattice_entropy(100))
    print(lattice.lattice_free_energy(100))
    print(lattice.lattice_energy(300))
