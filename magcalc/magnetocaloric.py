import numpy as np
from scipy.optimize import fsolve


class Magnetocaloric:
    # Boltzmann Constant
    k_B = 8.6173324 * (10 ** (-5))  # eV K^-1
    # Bohr magneton
    mu_B = 5.7883818066 * (10 ** (-5))  # eV T^-1

    def __init__(self, phase, J, gJ, Tc, theta_D, F0, Nm, N):
        self.phase = phase  # Phase of material
        self.J = J  # Total Angular Momentum, J
        self.gJ = gJ  # Landé g-factor
        self.Tc = Tc  # Curie Temperature, in Kelvin
        self.theta_D = theta_D  # Debye temperatures, in Kelvin
        self.F0 = F0  # Internal Energies of Lattice, in eV
        self.Nm = Nm  # Number of Magnetic Moments in Primitive Cell
        self.N = N  # Number of Atoms in Primitive Cell

        # Curie Constants divided by Vacuum Permeability (C/mu_0)
        self.const = (self.mu_B ** 2.) * self.Nm * (self.gJ ** 2.) * self.J * (self.J + 1.) / (3. * self.k_B)
        self.lamb = self.Tc/self.const  # Value of the strength of the parameter of the Molecular Field

    def get(self, arg):
        pass

    # Function for the computation of sigma
    def B_J(self, sigma, T, B):
        h = B / (self.lamb * self.Nm * self.gJ * self.mu_B * self.J)
        y = 3. * self.J / (self.J + 1.) * (h + sigma) * self.Tc / T
        return sigma - (2. * self.J + 1.) / (2. * self.J * np.tanh((2. * self.J + 1.) * y / (2. * self.J))) + \
               1. / (2. * self.J * np.tanh(y / (2 * self.J)))

    def reduced_magnetization(self, T, B):
        """Brillouin function. Calculates the reduced magnetization.

        Parameters
        ----------
        T : scalar
            An array with the temperatures
        B : scalar
            An array with the magnetic fields

        Returns
        --------
        y : scalar
            An array with the values of the reduced magnetization
        """
        sigma = np.zeros_like(T)  # variable that stores the values

        is_negative = False  # For negative fields
        if B < 0:
            is_negative = True
            # Change sign, do calculations por positive magnetic field
            B = -1 * B
        sigma = fsolve(B_J, 0.5, args=(T, B, J, TC, lamb))
        if is_negative:
            sigma = -1. * sigma  # Change back the sign

        return sigma

    def magnetic_entropy(self, T, B):
        """Computes the magnetic entropy for 1 spin.

        Parameters
        ----------
        T : scalar, 2D array
            Temperatures
        B : scalar, 2D array
            Magnetic fields

        Returns
        -------
        y : scalar, 2D array
            Magnetic entropy
        """
        h = B / (self.lamb * self.Nm * self.gJ * self.mu_B * self.J)  # relative field to the saturation magnetization
        sigma = self.reduced_magnetization(T, B, self.J, self.gJ, self.TC, self.lamb, self.Nm)  # reduced magnetization

        y = 3. * J / (J + 1.) * (h + sigma) * TC / T  # temporary variable

        A = np.sinh((2. * J + 1.) * y / (2. * J))
        B = np.sinh(y / (2. * J))

        return k_B * (np.log(A / B) - sigma * y)

    def entropy(self):
        pass

    def energy(self):
        pass

    def free_energy(self):
        pass