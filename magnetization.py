# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 17:03:14 2016

@author: Rui M. Costa
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import fmin

import free_energy as free

mu_B = 5.7883818066*(10**(-5)) # eV T^-1
k_B = 8.6173324*(10**(-5)) # eV K^-1


#==============================================================================
# Functions
#==============================================================================


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
        return sigma - (2.*J + 1.)/(2.*J*np.tanh((2.*J+1.)*y/(2.*J))) + 1./(2.*J*np.tanh(y/(2*J)))

    sigma = np.zeros_like(T) # variable that stores the values

    if sigma.shape == (): # if T and B are scalars
        is_negative = False # For negative fields
        if B < 0:
            is_negative = True
            B = -1*B # Change sign, do calculations por positive magnetic field

        sigma = fsolve(B_J, 0.5, args=(T, B, J, TC, lamb))
        if is_negative:
            sigma = -1.*sigma

    else: # if T and B are a 2D array
        B_range, T_range = sigma.shape
        for i in range(B_range): # calculate reduced magnetization for each temperature and magnetic field
            b = B[i,0]

            is_negative = False # For negative fields
            if b < 0:
                is_negative = True
                b = -1*b # Change sign, do calculations por positive magnetic field

            for j in range(T_range):
                sig_sol = fsolve(B_J, 0.5, args=(T[0,j], b, J, TC, lamb)) # magnetization

                if is_negative:
                    sigma[i,j] = -1.*sig_sol # if B < 0, change signal back to negative
                else:
                    sigma[i,j] = sig_sol # if B >= 0
    
    return sigma





# Reduced Magnetization of the Stable Phase
def Brillouin_stable(T, B, J1, J2, TC1, TC2, lamb1, lamb2, theta_D1, theta_D2,
                     F01, F02, gJ, Nm, N):
    """Computes the magnetization of the system considering the phase transition
    
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
    sigma1 = Brillouin(T, B, J1, gJ, TC1, lamb1, Nm) # magnetization of structure 1
    sigma2 = Brillouin(T, B, J2, gJ, TC2, lamb2, Nm) # magnetization of structure 2

    F1 = free.F(T, B, J1, gJ, TC1, lamb1, Nm, theta_D1, N, F01) # free energy of phase 1
    F2 = free.F(T, B, J2, gJ, TC2, lamb2, Nm, theta_D2, N, F02) # free energy of phase 2

    F_cross_index2 = (F1 > F2).astype(int) # determine index where F1 > F2
    F_cross_index1 = (F1 < F2).astype(int) # determine index where F1 < F2

    sigma_stable = F_cross_index1*sigma1 + F_cross_index2*sigma2 # magnetization of stable phase

    return sigma_stable



# Gaussian Distribution
def Gauss(T, mu, var): # T - Variable, mu - Average, var - Variance
    """Gaussian/normal distribution normalized to 1.
    
    Parameters
    ---------
    T : array
        Temperatures.
    mu : scalar
        Average temperature.
    var : scalar
        Variance of the temperature.    
        
    Returns
    -------
    y : array
        Gaussian/normal distribution.
    """
    return 1./np.sqrt(2*np.pi*var)*np.exp(-1.*((T - mu)**2.)/(2.*var))


# Reduced Magnetization with Gaussian Distribution of Tc's
def Brillouin_Gauss(T, B, J, TC, lamb, var):
    """Reduced magnetization for several Curie temperatures
    with a normal distribution.
    
    Parameters
    ---------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.
    J : scalar
        Total angular momentum.
    TC : scalar
        Value of the average Curie temperature.
    lamb : scalar
        Value of the strength of the parameter of the Molecular Field.
    
        
    Returns
    -------
    y : 2D array
        Reduced magnetization.
    """

    TCs = np.where(Gauss(T[0], TC, var) > 0.05*Gauss(TC, TC, var))[0] # index of temperatures

    TCs = T[0][TCs[0]:TCs[-1]+1] # temperatures

    sigma = np.zeros(np.shape(T))
    for tc in TCs:
        sigma += Gauss(tc, TC, var)*Brillouin(T, B, J, TC, lamb)

    return sigma


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
    
    
    if 1 == 0:
        print '\t Profiling...'
        import cProfile
        import Profiling as Prof
        import os
        file_name = "Profile_Output\ " + str(os.path.basename(__file__)) + '_Profile_Output'
        cProfile.runctx('Brillouin_stable(T, B)', {'Brillouin_stable':Brillouin_stable,'T':TT, 'B':BB}, {}, filename = file_name)
        Prof.save_profiling(file_name, sort='cumtime')
    

