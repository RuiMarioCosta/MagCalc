# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 21:08:56 2016

@author: Rui M. Costa
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import fmin

import magnetization as mag
import entropy as ent
import free_energy as free

mu_B = 5.7883818066*(10**(-5))  # eV T^-1
k_B = 8.6173324*(10**(-5))  # eV K^-1


def E_L(T, theta_D):
    """Computes the lattice energy.

    Parameters
    ---------
    T : array
        Temperatures.
    theta_D : scalar
        Debye temperature.

    Returns
    -------
    y : array
        Lattice energy.
    """
    # Function in the integral
    def f(x):
        return (x**3.)/(np.exp(x) - 1.)

    integral = np.zeros(T.shape)  # variable that stores the values

    for i, t in enumerate(T):  # calculate the integral for each temperature
        integral[i] = quad(f, 0., theta_D/t)[0]

    return k_B*(9./8.*theta_D + 9.*theta_D*((T/theta_D)**4.)*integral)


# Magnetic Energy of one Moment
def E_M(T, B, J, gJ, TC, lamb, Nm):
    """Computes the magnetic energy.

    Parameters
    ---------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.
    J : scalar
        Total angular momentum.
    TC : scalar
        Curie temperature.
    lamb : scalar
        Value of the strength of the parameter of the Molecular Field.

    Returns
    -------
    y : array
        Magnetic energy.
    """
    sigma = mag.Brillouin(T, B, J, TC, lamb)  # reduced magnetization

    return -gJ*mu_B*J*B*sigma - 3.*J/(J + 1.)*k_B*TC*(sigma**2.)


# Total Energy of Unit Cell
def E_tot(T, B, J, TC, lamb, theta_D, F0):

    el = np.ones(np.shape(T))*E_L(Delta_T, theta_D)

    return F0 + Nm*E_M(T, B, J, TC, lamb) + N*el


# Magnetic Energy as a function of Magnetization
def E_M_vs_M(sigma, T, B, J, TC, lamb):
    """Computes the magnetic energy as a function of the reduced magnetization.

    Parameters
    ---------
    sigma : array
        Reduced magnetization.
    T : scalar
        Temperatures.
    B : scalar
        Magnetic fields.
    J : scalar
        Total angular momentum.
    TC : scalar
        Curie temperature.
    lamb : scalar
        Value of the strength of the parameter of the Molecular Field.

    Returns
    -------
    y : array
        Magnetic energy.
    """
    def f(sigma, T, B, J, TC, lamb):
        A = np.sinh(3./(2.*(J+1.))*(B/(lamb*Nm*gJ*mu_B*J) + sigma)*TC/T)
        B = np.sinh(3.*(2.*J+1.)/(2.*(J+1.))*(B/(lamb*Nm*gJ*mu_B*J) + sigma)*TC/T)

        return (sigma**2.)/2. + (J+1.)/(3.*J)*T/TC*np.log(A/B)


    Tt, Bb = np.meshgrid(np.array([T]), np.array([B]))

    sigma0 = mag.Brillouin(Tt, Bb, J, TC, lamb) # reduced magnetization of average minimum

    h = Bb/(lamb*Nm*gJ*mu_B*J) # relative field to the saturation magnetization
    y = 3.*J/(J+1.)*(h + sigma0)*TC/Tt
    C = np.sinh((2.*J + 1.)*y/(2.*J))
    D = np.sinh(y/(2.*J))
    F0 = -k_B*T*np.log(C/D) # free energy of average magnetization

    F1 = 3./(2.*(J+1.))*(B/(lamb*Nm*gJ*mu_B*J) + sigma)*TC/T
    F2 = 3.*(2.*J+1.)/(2.*(J+1.))*(B/(lamb*Nm*gJ*mu_B*J) + sigma)*TC/T

    return (sigma**2.)/2. + 1./(2.*J)*(B/(lamb*Nm*gJ*mu_B*J) + sigma)*(1./np.tanh(F1) - (2.*J + 1.)/np.tanh(F2)) - f(sigma0, T, B, J, TC, lamb) + F0


if __name__ == "__main__":
    print 'energy.py'
