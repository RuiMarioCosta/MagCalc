# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 17:03:14 2016

@author: Rui M. Costa
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import fmin

import free_energy as free

mu_B = 5.7883818066*(10**(-5))  # eV T^-1
k_B = 8.6173324*(10**(-5))  # eV K^-1


def Brillouin(T, B, J, gJ, TC, lamb, Nm):
    """Brillouin function. Calculates the reduced magnetization.

    Parameters
    ----------
    T : scalar, 2D array
        An array with the temperatures.
    B : scalar, 2D array
        An array with the magnetic fields.
    J : scalar
        Value of angular momentum.
    TC : scalar
        Value of Curie temperature.
    lamb : scalar
        Value of the strength of the parameter of the Molecular Field.
    gJ : scaclar
        Landé g-factor.
    Nm : scalar
        Number os spins.

    Returns
    --------
    y : scalar, array
        An array with the values of the reduced magnetization.
    """

    # Function for the computation of sigma
    def B_J(sigma, T, B, J, TC, lamb):
        h = B/(lamb*Nm*gJ*mu_B*J)
        y = 3.*J/(J+1.)*(h + sigma)*TC/T
        return sigma - (2.*J+1.)/(2.*J*np.tanh((2.*J+1.)*y/(2.*J))) + 1./(2.*J*np.tanh(y/(2*J)))

    sigma = np.zeros_like(T)  # variable that stores the values

    if np.isscalar(T) and np.isscalar(B):  # if T and B are scalars
        is_negative = False  # For negative fields
        if B < 0:
            is_negative = True
            # Change sign, do calculations por positive magnetic field
            B = -1*B

        sigma = fsolve(B_J, 0.5, args=(T, B, J, TC, lamb))
        if is_negative:
            sigma = -1.*sigma  # Change back the sign

    elif len(sigma.shape) == 2:  # if T and B are a 2D array
        B_range, T_range = sigma.shape
        # calculate reduced magnetization for each T and B
        for i in range(B_range):
            b = B[i, 0]

            is_negative = False  # For negative fields
            if b < 0:
                is_negative = True
                # Change sign, do calculations por positive magnetic field
                b = -1*b

            for j in range(T_range):
                # magnetization
                sig_sol = fsolve(B_J, 0.5, args=(T[0, j], b, J, TC, lamb))

                if is_negative:
                    # if B < 0, change signal back to negative
                    sigma[i, j] = -1.*sig_sol
                else:
                    sigma[i, j] = sig_sol  # if B >= 0
    else:
        raise Exception('The variables T and B should be both scalars or both 2D arrays.')
    return sigma



# Reduced Magnetization of the Stable Phase
def Brillouin_stable(T, B, J1, J2, TC1, TC2, lamb1, lamb2, theta_D1, theta_D2, F01, F02, gJ, Nm, N):
    """Returns the magnetization of the system with lowest free energy, i.e., the stable one.

    Parameters
    ----------
    T : scalar, 2D array
        Temperature.
    B : scalar, 2D array
        Magnetic fields.
    J1, J2 : scalar
        Value of angular momentum.
    TC1, TC2 : scalar
        Value of Curie temperatures.
    lamb1, lamb2 : scalar
        Value of the strength of the parameters of the Molecular Field.
    theta_D1, theta_D2 : scalar
        Debye temperatures of the material.
    F01, F02 : scalar
        Free eneries at 0 K.

    Returns
    ---------
    y : scalar, 2D array
        Magnetization corresponding to the stable phase.
    """
    sigma1 = Brillouin(T, B, J1, gJ, TC1, lamb1, Nm) # magnetization of phase 1
    sigma2 = Brillouin(T, B, J2, gJ, TC2, lamb2, Nm) # magnetization of phase 2

    F1 = free.F(T, B, J1, gJ, TC1, lamb1, Nm, theta_D1, N, F01) # free energy of phase 1
    F2 = free.F(T, B, J2, gJ, TC2, lamb2, Nm, theta_D2, N, F02) # free energy of phase 2

    F_cross_index2 = (F1 > F2).astype(int) # determine index where F1 > F2
    F_cross_index1 = (F1 < F2).astype(int) # determine index where F1 < F2

    sigma_stable = F_cross_index1*sigma1 + F_cross_index2*sigma2 # magnetization of stable phase

    return sigma_stable



# Magnetic Susceptibility
def Susceptibility(T, B, J, TC, lamb):
    sigma = Brillouin(T, B, J, TC, lamb)

    mu_0 = 4.*np.pi*(10**(-7))

    return mu_0*sigma/B



# Reduced Magnetization as a Function of Applied Field on Heating
def RedMag_heat(T, B, J1, TC1, lamb1, theta_D1, F01, J2, TC2, lamb2, theta_D2, F02, gJ, Nm, N):
    """Reduced magnetization as a function of magnetic field on heating.

    Parameters
    ---------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields

    Returns
    -------
    y : 2D array
        Reduced magnetization on heating.
    """
    sigmaheat = np.zeros_like(T)

    B_range, T_range = T.shape
    for i in range(B_range):
        sigma_guess = 1. #initial guess when heating, assuming it starts at low temperatures
        for j in range(T_range):
            sol = fmin(free.F_totstable_vs_M,
                       sigma_guess,
                       args=(T[0,j], B[i,0], J1, TC1, lamb1, theta_D1, F01, J2, TC2, lamb2, theta_D2, F02, gJ, Nm, N),
                       full_output=1, disp=0)
            sigma_guess = sol[0] # new guess = last magnetization
            sigmaheat[i,j] = sol[0] # value of free energy at local minimum

    return sigmaheat



# Reduced Magnetization as a Function of Applied Field on Cooling
def RedMag_cool(T, B, J1, TC1, lamb1, theta_D1, F01, J2, TC2, lamb2, theta_D2, F02, gJ, Nm, N):
    """Reduced magnetization as a function of magnetic field on cooling.

    Parameters
    ---------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields

    Returns
    -------
    y : 2D array
        Reduced magnetization on cooling.
    """
    sigmacool = np.zeros_like(T)

    B_range, T_range = T.shape
    for i in range(B_range):
        sigma_guess = 0. #initial guess when cooling, assuming it starts at high temperatures
        for j in range(T_range-1,-1,-1):
            sol = fmin(free.F_totstable_vs_M,
                       sigma_guess,
                       args=(T[0,j], B[i,0], J1, TC1, lamb1, theta_D1, F01, J2, TC2, lamb2, theta_D2, F02, gJ, Nm, N),
                       full_output=1, disp=0)
            sigma_guess = sol[0] + 0.001 # new guess = last magnetization
            sigmacool[i,j] = sol[0] # value of free energy at local minimum

    return sigmacool




if __name__ == "__main__":
    from variables import * # import variables to use in functions

    print Brillouin(5., 1., J1, gJ, TC1, lamb1, Nm)
    #print Brillouin(5., Delta_B, J1, gJ, TC1, lamb1, Nm)
    #print Brillouin(TT, BB, J1, gJ, TC1, lamb1, Nm)

    print Brillouin_stable(100., 1., J1, J2, TC1, TC2, lamb1, lamb2, theta_D1, theta_D2, F01, F02, gJ, Nm, N)
