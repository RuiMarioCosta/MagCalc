# -*- coding: utf-8 -*-
"""
Created on Tue Apr 05 00:02:29 2016

@author: Rui M. Costa
"""

import numpy as np


#==============================================================================
# Constants
#==============================================================================

# Boltzmann Constant
k_B = 8.6173324*(10**(-5))  # eV K^-1


# Bohr magneton
mu_B = 5.7883818066*(10**(-5)) # eV T^-1


#==============================================================================
# Parameters
#==============================================================================

# Gd5Si2Ge2 -> 1
# Tb5Si2Ge2 -> 2
# Er5Si4 -> 3

Material = 1

if Material == 1:    ## Gd5Si2Ge2

    # Total Angular Momentum, J
    J1 = 7/2.
    J2 = 7/2.

    # g-factor
    gJ = 2.


    # Curie Temperature, in Kelvin
    TC1 = 251.   # M
    TC2 = 308.   # O(I)


    # Debye temperatures, in Kelvin
    theta_D1 = 250. # M
    theta_D2 = 278. # O(I)


    # Free Energies at 0 K, in eV
    F01 = 0.36   # M
    F02 = 0.    # O(I)

    # Conversion from eV/K to J/Kg K to plot entropy changes.
    Conv = 24419.523 # Conversion from eV/K to J/Kg K.



if Material == 2:    ## Tb5Si2Ge2

    # Total Angular Momentum, J
    J1 = 6.
    J2 = 6.

    # g-factor
    gJ = 3/2.


    # Curie Temperature, in Kelvin
    TC1 = 112. # M
    TC2 = 200. # O(I)


    # Debye temperatures, in Kelvin
    theta_D1 = 153. # M
    theta_D2 = 170. # O(I)


    # Internal Energies of Lattice, in eV
    F01 = 0.43 # 0.11 # M         0.43
    F02 = 0. # O(I)


if Material == 3:    ## Er5Si4

    # Total Angular Momentum, J
    J1 = 15/2.
    J2 = 15/2.

    # g-factor
    gJ = 6/5.


    # Curie Temperature, in Kelvin
    TC1 =  30. # M
    TC2 =  38.5  # O(I)


    # Debye temperatures, in Kelvin
    theta_D1 = 200.   #405.91   # M
    theta_D2 = 175.   #385.     # O(I)


    # Internal Energies of Lattice, in eV
    F01 = 0. # M
    F02 = 0.143911006   # O(I)  0.143911006    0.114426506


if Material == 4:    ## Gd5Ge4

    # Total Angular Momentum, J
    J1 = 7/2.
    J2 = 7/2.

    # g-factor
    gJ = 2.


    # Curie Temperature, in Kelvin
    TC1 = 35.   # O(II)
    TC2 = 245.  # O(I)


    # Debye temperatures, in Kelvin
    theta_D1 = 250.  # O(II)
    theta_D2 = 278.  # O(I)


    # Free Energies at 0 K, in eV
    F01 = 0.   # O(II)
    F02 = 0.3  # O(I)


#------------------------------------------------------------------------------


# Number of Magnetic Moments in Primitive Cell
Nm = 20.

# Number of Atoms in Primitive Cell
N = 36.


# Initial, Final and Step Temperatures, in Kelvin
Ti = 250.
Tf = 350.
delta_T = 1.
Tf = Tf + delta_T

# Initial, Final and Step Applied Magnetic Fields, in Tesla
Bi = 0.
Bf = 5.
delta_B = 1.
Bf = Bf + delta_B


#==============================================================================
# Other Variables (used in the program)
#==============================================================================

# Temperature interval, in Kelvin
Delta_T = np.arange(Ti, Tf, delta_T)


# Temperature interval, in Kelvin
Delta_B = np.arange(Bi, Bf, delta_B)
#Delta_B = np.array([0., 2., 5.])

# Domain, Grid with Temperature and Magnetic Field Values
TT, BB = np.meshgrid(Delta_T, Delta_B)


# Curie Constants divided by Vacuum Permeability (C/mu_0)
C1 = (mu_B**2.)*Nm*(gJ**2.)*J1*(J1 + 1.)/(3.*k_B)
C2 = (mu_B**2.)*Nm*(gJ**2.)*J2*(J2 + 1.)/(3.*k_B)


# Parameter of the strength of the Molecular Field, lambda
lamb1 = TC1/C1
lamb2 = TC2/C2




