# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 17:42:01 2016

@author: Rui M. Costa
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelmax, argrelmin
from scipy.integrate import quad
from scipy.optimize import fmin

import magnetization as mag
import entropy as ent
import energy as ener

mu_B = 5.7883818066*(10**(-5)) # eV T^-1
k_B = 8.6173324*(10**(-5)) # eV K^-1


#==============================================================================
# Functions
#==============================================================================


def F_M(T, B, J, gJ, TC, lamb, Nm):
    """Magnetic free energy of 1 spin.

    Parameters
    ----------
    T : 2D array
        Temepratures.
    B : 2D array
        Magnetic fields.
    J : scalar
        Total angular momentum.
    gJ : scalar
        Landé g-factor.
    TC : scalar
        Curie temperature.
    lamb : scalar
        Value of the strength of the parameter of the Molecular Field.
    Nm : scalar
        Number of spins.

    Returns
    -------
    y : 2D array
        Magnetic free energy.
    """
    h = B/(lamb*Nm*gJ*mu_B*J) # relative field to the saturation magnetization
    sigma = mag.Brillouin(T, B, J, gJ, TC, lamb, Nm) # reduced magnetization
    y = 3.*J/(J+1.)*(h + sigma)*TC/T
    A = np.sinh((2.*J + 1.)*y/(2.*J))
    B = np.sinh(y/(2.*J))

    return -T*k_B*np.log(A/B) # magnetic free energy of 1 spin



# Lattice Free Energy
def F_L(T, theta_D):
    """Free Lattice Energy according to the Debye Model.

    Parameters
    ----------
    T : scalar, 1D array
        Temperature.
    theta_D : scalar
        Debye temperature of the material

    Returns
    --------
    y : scalar, array
        Free Lattice Energy
    """

    # Function in the integral
    def f(x):
        return (x**3.)/(np.exp(x) - 1.)


    integral = np.zeros_like(T) # variable that stores the values

    if np.isscalar(T): # if T is just a single Temperature
        integral = quad(f, 0., theta_D/T)[0]
    elif len(T.shape) == 1: #if T is an array of Temperatures
        for i,t in enumerate(T): # calculate the integral for each temperature
            integral[i] = quad(f, 0., theta_D/t)[0]
    else:
        raise Exception('The variables T should be a scalar or a 1D array.')


    return k_B*(9./8.*theta_D - 3.*T*((T/theta_D)**3.)*integral + 3.*T*np.log(1. - np.exp(-theta_D/T)))



# Total Free Energy
def F2(T, B, J, gJ, TC, lamb, Nm, theta_D, N, F0):
    """Total free energy as a functions of temperature and magnetic field in the
    unit cell.

    Parameters
    ---------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.
    J : scalar
        Total angular momentum.
    theta_D : scalar
        Debye temperature.
    F0 : scalar
        Electronic free enery at 0 K.
    lamb : scalar
        Value of the strength of the parameter of the Molecular Field.

    Returns
    -------
    y : 2D array
        Total free energy.
    """
    sigma = mag.Brillouin(T, B, J, TC, lamb) # reduced magnetization

    s_M = ent.S_M(T, B, J, TC, lamb) # magnetic entropy from MFT

    s_L = ent.S_L(T[0], theta_D) # lattice entropy from Debye model
    s_L = s_L*np.ones(np.shape(T)) # turn array into a matrix with equal rows for matrix multiplication

    f0 = F0*np.ones(T.shape)
    F_0 = Nm*(3.*J/(J + 1.)*k_B*TC) + f0 # offset of free energies to make F start at 0 with no magnetic field

    return F_0 + Nm*(-gJ*mu_B*J*B*sigma - 3.*J/(J + 1.)*k_B*TC*(sigma**2.)- T*s_M) - N*T*s_L



# Free Energy (faster calculation)
def F(T, B, J, gJ, TC, lamb, Nm, theta_D, N, F0):
    """Total free energy as a functions of temperature and magnetic field in the
    unit cell.

    Parameters
    ---------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.
    J : scalar
        Total angular momentum.
    theta_D : scalar
        Debye temperature.
    F0 : scalar
        Electronic free enery at 0 K.
    lamb : scalar
        Value of the strength of the parameter of the Molecular Field.

    Returns
    -------
    y : 2D array
        Total free energy.
    """
    f_M = F_M(T, B, J, gJ, TC, lamb, Nm) # magnetic free energy of 1 spin

    f_L = F_L(T[0], theta_D)*np.ones_like(T) # lattice free energy of 1 atom

    #h = B/(lamb*Nm*gJ*mu_B*J) # relative field to the saturation magnetization
    F_0 = Nm*(3.*J/(J + 1.)*(0. + 1.)*k_B*TC) - N*k_B*9./8.*theta_D  + F0 # offset of free energies to make F start at F0

    return  Nm*f_M + N*f_L + F_0



def transition_temp(T, B, J1, TC1, theta_D1, F01, lamb1, J2, TC2, theta_D2, F02, lamb2, gJ, Nm, N, *args):
    """Calculates the transition temperatures.

    Parameters
    ----------
    T : 2D array
            Temperature
    B : 2D array
            Applied Magnetic Field
    args : tuple
            Total free energies of both structures, (F1,F2)

    Returns
    ----------
    y : tuple
            Transition temperatures.
    """
    Ts = np.zeros_like(B[:,0])

    if not args: # if the free energies weren't calculated
        diff = F(T, B, J1, gJ, TC1, lamb1, Nm, theta_D1, N, F01) - F(T, B, J2, gJ, TC2, lamb2, Nm, theta_D2, N, F02)
    else: # if the free energies are passed as arguments
        diff = args[0] - args[1]

    for j in range(T.shape[0]):
        for i in range(T.shape[1] - 1):
            if diff[j,i] == 0. or diff[j,i] * diff[j,i + 1] < 0.: # crossover at i
                Ts[j] = T[0,i]

    return Ts


# Free energy of the Stable Phase
def F_stable(T, B, J1, J2, gJ, TC1, TC2, lamb1, lamb2, Nm, theta_D1, theta_D2, N, F01, F02):
    """Free energy following the stable phase.

    Parameters
    ----------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.

    Returns
    -------
    y : Stable free energy.
    """
    F1 = F(T, B, J1, gJ, TC1, lamb1, Nm, theta_D1, N, F01) # free energy of pahse 1
    F2 = F(T, B, J2, gJ, TC2, lamb2, Nm, theta_D2, N, F02) # free energy of phase 2

    return np.minimum(F1,F2)





# Magnetic Free Energy as a function of Magnetization
def F_M_vs_M(sigma, T, B, J, gJ, TC, lamb, Nm):
    """Magnetic Free Energy as a functio of Reduced Magnetization.

    Parameters
    ----------
    sigma : scalar, 2D array
        Reduced Magnetization.
    T : scalar, 2D array
        Temperature.
    B : scalar, 2D array
        Applied Magnetic Field.
    J : scalar
        Total Angular Momentum.
    TC : scalar
        Curie Temperature.
    lamb : scalar
        Value of the strength of the parameter of the Molecular Field.

    Returns
    --------
    y : scalar, array
        Magnetic Free Energy
    """

    def f(sigma, T, B, J, TC, lamb):
        A = np.sinh(3./(2.*(J+1.))*(B/(lamb*Nm*gJ*mu_B*J) + sigma)*TC/T)
        B = np.sinh(3.*(2.*J+1.)/(2.*(J+1.))*(B/(lamb*Nm*gJ*mu_B*J) + sigma)*TC/T)

        return (sigma**2.)/2. + (J+1.)/(3.*J)*T/TC*np.log(A/B)


    sigma0 = mag.Brillouin(T, B, J, gJ, TC, lamb, Nm) # reduced magnetization of absolute minimum

    h = B/(lamb*Nm*gJ*mu_B*J) # relative field to the saturation magnetization
    y = 3.*J/(J+1.)*(h + sigma0)*TC/T
    C = np.sinh((2.*J + 1.)*y/(2.*J))
    D = np.sinh(y/(2.*J))
    F0 = -k_B*T*np.log(C/D) # value of free energy at absolute minimum

    return  f(sigma, T, B, J, TC, lamb) - f(sigma0, T, B, J, TC, lamb) + F0


# Magnetic Free Energy of Stable Phase as a function of Magnetization
def F_Mstable_vs_M(sigma, T, B, J1, TC1, lamb1, J2, TC2, lamb2):
    """Stable magnetic free energy as a functio of reduced magnetization.

    Parameters
    ----------
    sigma : scalar, 1D array
        Reduced Magnetization.
    T : scalar, 1D array
        Temperature.
    B : scalar, 1D array
        Applied Magnetic Field.

    Returns
    --------
    y : scalar, array
        Stable magnetic free energy.
    """
    F1 = F_M_vs_M(sigma, T, B, J1, TC1, lamb1)
    F2 = F_M_vs_M(sigma, T, B, J2, TC2, lamb2)

    return np.minimum(F1,F2)



# Total Free Energy as a function of Magnetization
def F_tot_vs_M(sigma, T, B, J, gJ, TC, lamb, Nm, theta_D, N, F0):
    """Total free energy as a function of reduced magnetization.

    Parameters
    ----------
    sigma : 1D array
        Reduced magnetization, between -1 and 1.
    T : scalar, 1D array
        Temperature.
    B : scalar
        Applied magnetic field.
    J : scalar
        Total angular momentum.
    TC : scalar
        Curie temperature.
    lamb : scalar
        Value of the strength of the parameter of the Molecular Field.
    theta_D : scalar
        Debye temperature.
    F0 : scalar
        Electronic free energy.

    Returns
    --------
    y : scalar, array
        Total free energy.
    """

    fM = F_M_vs_M(sigma, T, B, J, gJ, TC, lamb, Nm) # magnetic free energy
    fL = F_L(T, theta_D)*np.ones(np.shape(fM)) # lattice free energy
    F_0 = Nm*(3.*J/(J + 1.)*1.*k_B*TC) - N*k_B*9./8.*theta_D  + F0 # offset of free energies to make F start at F0
    return fM*Nm + N*fL + F_0


# Total Free Energy of Stable Phase as a function of Magnetization
def F_totstable_vs_M(sigma, T, B, J1, TC1, lamb1, theta_D1, F01, J2, TC2, lamb2, theta_D2, F02, gJ, Nm, N, *args):
    """For every magnetization (sigma), computes the minimum between the 2 structures

    Parameters
    ----------
    sigma : array
            Reduced magnetization, from -1 to 1
    T : scalar
            Temperature
    B : scalar
            Applied magnetic field
    args : tuple
            Total free energies of both structures, (F1,F2)

    Returns
    --------
    y : array
            Array containing the values of the minimum free energy.
    """
    if not args:
        F1 = F_tot_vs_M(sigma, T, B, J1, gJ, TC1, lamb1, Nm, theta_D1, N, F01)
        F2 = F_tot_vs_M(sigma, T, B, J2, gJ, TC2, lamb2, Nm, theta_D2, N, F02)
    else:
        F1, F2 = args

    return np.minimum(F1,F2)


# Total Free Energy on Heating
def F_tot_stable_Heating(T, B, J1, TC1, lamb1, theta_D1, F01, J2, TC2, lamb2, theta_D2, F02, gJ, Nm, N):
    """Free energy of the stable phase following the metastable minimum on heating.

    Parameters
    ---------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.

    Returns
    -------
    y : 2D array
        Free energy on heating.
    """
    f_heat = np.zeros(np.shape(T))

    Delta_B = B[:,0]
    Delta_T = T[0]

    for i in range(len(Delta_B)):
        sigma_guess = 1. #initial guess when heating, assuming it starts at low temperatures
        for j in range(len(Delta_T)):
            sol = fmin(F_totstable_vs_M,
                       sigma_guess,
                       args=(Delta_T[j], Delta_B[i], J1, TC1, lamb1, theta_D1, F01, J2, TC2, lamb2, theta_D2, F02, gJ, Nm, N),
                       full_output=1, disp=0)
            sigma_guess = sol[0] # new guess = last magnetization
            f_heat[i,j] = sol[1] # value of free energy at local minimum

    return f_heat


# Total Free Energy on Cooling
def F_tot_stable_Cooling(T, B, J1, TC1, lamb1, theta_D1, F01, J2, TC2, lamb2, theta_D2, F02, gJ, Nm, N):
    """Free energy of the stable phase following the metastable minimum on cooling.

    Parameters
    ---------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.

    Returns
    -------
    y : 2D array
        Free energy on heating.
    """
    f_cool = np.zeros(np.shape(T))

    Delta_B = B[:,0]
    Delta_T = T[0]

    for i in range(len(Delta_B)):
        sigma_guess = 0. #initial guess when cooling, assuming it starts at high temperatures
        for j in range(len(Delta_T)-1,-1,-1):
            sol = fmin(F_totstable_vs_M,
                       sigma_guess,
                       args=(Delta_T[j], Delta_B[i], J1, TC1, lamb1, theta_D1, F01, J2, TC2, lamb2, theta_D2, F02, gJ, Nm, N),
                       full_output=1, disp=0)
            sigma_guess = sol[0] + 0.001 # new guess = last magnetization
            f_cool[i,j] = sol[1] # value of free energy at local minimum

    return f_cool


#==============================================================================


if __name__ == "__main__":
    from variables import *

    # print F_M(5., 1., J1, gJ, TC1, lamb1, Nm)
    # print F_M(5., Delta_B, J1, gJ, TC1, lamb1, Nm) # should raise and exception
    # print F_M(TT, BB, J1, gJ, TC1, lamb1, Nm)

    # print F_L(5., theta_D1)
    # print F_L(Delta_T, theta_D1)
    # print F_L(TT, theta_D1) # should raise and exception

    # print F(5., 1., J1, gJ, TC1, lamb1, Nm, theta_D1, N, F01)
    # print F(TT, BB, J1, gJ, TC1, lamb1, Nm, theta_D1, N, F01)

    # print transition_temp(TT, BB, J1, TC1, theta_D1, F01, lamb1, J2, TC2, theta_D2, F02, lamb2, gJ, Nm, N)

    # print F_stable(TT, BB, J1, J2, gJ, TC1, TC2, lamb1, lamb2, Nm, theta_D1, theta_D2, N, F01, F02)

    # print F_M_vs_M(sig, 5., 0., J1, gJ, TC1, lamb1, Nm)
    print F_M_vs_M(sig, TT, BB, J1, gJ, TC1, lamb1, Nm)


    # f_0T = F_totstable_vs_M(sig, T0, 0.)
    # f = F_totstable_vs_M(sig, T0, B0)

    # mins_0T = argrelmin(f_0T)[0]
    # mins = argrelmin(f)[0]
    # maxs_0T = argrelmax(f_0T)[0]
    # maxs = argrelmax(f)[0]

    # print mins_0T, f_0T[mins_0T], sig[mins_0T]
    # print mins , f[mins], sig[mins]
    # print maxs_0T
    # print maxs
