# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 14:37:46 2016

@author: Rui M. Costa
"""

import numpy as np
import matplotlib.pyplot as plt

import Magnetization as mag
import Entropy as ent
from Variables import *
import FreeEnergy as free

from scipy.signal import argrelmax, argrelmin, sawtooth
from scipy.integrate import odeint, quad, simps
from scipy.optimize import fmin, fsolve


#==============================================================================
# Functions
#==============================================================================


# Applied Magnetic Field
def B_applied(t, A, w, t0):
    """Applied magnetic field. Input signal.
    
    Parameters
    ----------
    t : array
        A sequence of time points for which to solve for y.  The initial
        value point should be the first element of this sequence.
    A : scalar
        Amplitude of the applied magnetic field.
    w : scalar
        Frequency of input signal.
    t0 : scalar
        Time delay until sinusoidal 


    Returns
    -------
    y : array
        Array containing the values of B_applied for each desired time in t.
    """
    #return A*np.sign(t+1.)
    return 0.5*(np.sign(t-t0) + 1.)*A*np.sin(w*(t - t0))
    #return sawtooth(w*t, width = 0.5)*0.5*(np.sign(t-0.5*np.pi/w) + 1.)*A
    #return A*np.sin(w*t - np.pi/2.) + A
    #return 0.5*(np.sign(1000. - t) + 1.)*A*np.sin(w*(t - t0))

# Applied Temperature
def T_applied(t, A, w, t0):
    """Applied temperature. Input signal.
    
    Parameters
    ----------
    t : array
        A sequence of time points for which to solve for y.  The initial
        value point should be the first element of this sequence.
    A : scalar
        Amplitude of the applied temperature.
    w : scalar
        Frequency of input signal.
    t0 : scalar
        Time delay until sinusoidal 


    Returns
    -------
    T_aplied : array
        Array containing the values of B_applied for each desired time in t.
    """
    return A*np.sign(t+1.)



def free_energy_intersections(sigma, T, B, *args):
    """Calculates the indices where the intersection of both free energies occurs.
    
    Parameters
    ----------
    sigma : array
            Reduced magnetization from -1 to 1
    T : scalar
            Temperature
    B : scalar
            Applied Magnetic Field
    args : tuple, optional
            (Total Free Energy of Structure 1, Total Free Energy of Structure 2)
    
    
    Returns
    ----------
    y : tuple
            Indices of intersection
    """
    roots = []
    
    if not args: # if the free energies weren't calculated
        diff = free.F_tot_vs_M(sigma, T, B, J1, TC1, lamb1, theta_D1, F01) - free.F_tot_vs_M(sigma, T, B, J2, TC2, lamb2, theta_D2, F02)
    else: # if the free energies are passed
        diff = args[0] - args[1]
        
        
    for i in range(len(diff) - 1):
        if diff[i] == 0. or diff[i] * diff[i + 1] < 0.:
            # crossover at i
            roots.append(i)
        
    return roots


    

def norm_factor(sigma, T, B, *args):
    """Calculates de normalizing factor for every minimum.
    $$\sum e^{-F(\sigma, T, B)}$$
    
    Parameters
    ----------
    sigma : array
            Reduced magnetization from -1 to 
    T : scalar
            Temperature
    B : scalar
            Applied magnetic Field
    args : tuple, optional
            (Total Free Energy of Structure 1, Total Free Energy of Structure 2)
    
    Returns
    -----------
    y : tuple
            Normalizing factors
    """
    
    if not args: # if the free energies weren't calculated
        fm = free.F_totstable_vs_M(sigma, T, B) # total free energy
    else: # if the free energies are passed
        fm = args[0] # total free energy
        
    
    dsigma = sigma[1] - sigma[0] # step
    
    min_indices = argrelmin(fm)[0] # indices of minimums
    
    if len(min_indices) == 1: # if there is only 1 minimum
        fm = np.exp(-fm/(k_B*T)) # reuse variable to build function to integrate
        return simps(fm, sigma, dsigma) # normalization factor is the sum of the whole domain
    
    elif len(min_indices) == 2: # if there are 2 minimums
        maximum = argrelmax(fm)[0] # maximum that separates the two states
        
        fm = np.exp(-fm/(k_B*T)) # reuse variable to build function to integrate
        int_left = simps(fm[:maximum], sigma[:maximum], dsigma) # left normalizing factor
        int_right = simps(fm[maximum:], sigma[maximum:], dsigma) # right normalizing factor
        return int_left, int_right
    
    elif len(min_indices) == 3: # if there are 3 minimums
        max_index_left, max_index_right = argrelmax(fm)[0] # find indices of maximums

        fm = np.exp(-fm/(k_B*T)) # reuse variable to build function to integrate
        int_left = simps(fm[:max_index_left], sigma[:max_index_left], dsigma) # left nomalizing factor
        int_right = simps(fm[max_index_right:], sigma[max_index_right:], dsigma) # right normalizing factor
        int_mid = simps(fm[max_index_left:max_index_right], sigma[max_index_left:max_index_right], dsigma) # middle normalizing factor
        
        return int_left, int_mid, int_right

    elif len(min_indices) == 4:
        max_index_left, max_index_mid, max_index_right = argrelmax(fm)[0] # find indices of maximums

        fm = np.exp(-fm/(k_B*T)) # reuse variable to build function to integrate
        int_left = simps(fm[:max_index_left], sigma[:max_index_left], dsigma) # left nomalizing factor
        int_right = simps(fm[max_index_right:], sigma[max_index_right:], dsigma) # right normalizing factor
        int_mid_left = simps(fm[max_index_left:max_index_mid], sigma[max_index_left:max_index_mid], dsigma) # middle left normalizing factor
        int_mid_right = simps(fm[max_index_mid:max_index_right], sigma[max_index_mid:max_index_right], dsigma) # middle right normalizing factor
        
        return int_left, int_mid_left, int_mid_right, int_right
        




def diff_eq(z, t):
    """System of differential equation that describe the evolution of the phase
    fractions.
    
    Parameters
    ----------
    z : array
        Fractions.
    t : array
        Time.
    """
    z2_left, z1_left, z1_mid, z1_right, z2_right = z # unpack z
    
    Bb = B_applied(t, B0, w, t0) # applied magnetic field
    Tt = T_applied(t, T0, w, t0) # applied temperature
    
    f1 = free.F_tot_vs_M(sig, Tt, Bb, J1, TC1, lamb1, theta_D1, F01)
    f2 = free.F_tot_vs_M(sig, Tt, Bb, J2, TC2, lamb2, theta_D2, F02)
    fm = free.F_totstable_vs_M(sig, Tt, Bb, f1, f2) # free energy function
    
    norm_factor_sums = norm_factor(sig, Tt, Bb, fm) # normalizing factors
    
    min_indices = argrelmin(fm)[0] # indices of minimums
    if len(min_indices) == 3: # if there are 3 minimums
        
        max_index_left, max_index_right = argrelmax(fm)[0] # find indices of maximums
     
        C2_minus = np.exp(-fm[max_index_left]/(k_B*Tt)) # left maximum
        C2_plus = np.exp(-fm[max_index_right]/(k_B*Tt)) # right maximum
        
        C3_minus, C3_mid, C3_plus = norm_factor_sums # normalizing factors
                
        p_minus_mid = C2_minus/C3_minus # probability of a "particle" in the left minimum to be on top of the left barrier
        p_mid_minus = C2_minus/C3_mid # probability of a "particle" in the middle minimum to be on top of the left barrier
        p_mid_plus = C2_plus/C3_mid # probability of a "particle" in the middle minimum to be on top of the right barrier
        p_plus_mid = C2_plus/C3_plus # probability of a "particle" in the right minimum to be on top of the right barrier
        
        
        if Tt >= TC1: # if it's above the Curie temperature then the middle minimum will be z1_mid
            z1_mid = 1. - z2_left - z2_right # constraint
            
            # system of differential equations (rate of change of phase fractions)
            dz2_left = -p_minus_mid*z2_left + z1_mid*p_mid_minus
            dz2_right = -p_plus_mid*z2_right + z1_mid*p_mid_plus
            dz1_left = 0.
            dz1_right = 0.
            dz1_mid = - dz2_left - dz2_right

            return [dz2_left, dz1_left, dz1_mid, dz1_right, dz2_right]
            
        elif Tt < TC1 and Bb > 0.: # if it's below the Curie temperature with B>0 then the middle minimum will be z1_right
            z1_right = 1. - z2_left - z2_right # constraint
            
            # system of differential equations (rate of change of phase fractions)
            dz2_left = -p_minus_mid*z2_left + z1_right*p_mid_minus
            dz2_right = -p_plus_mid*z2_right + z1_right*p_mid_plus
            dz1_left = 0.
            dz1_mid = 0.
            dz1_right = - dz2_left - dz2_right
                        
            return [dz2_left, dz1_left, dz1_mid, dz1_right, dz2_right]
            
        elif Tt < TC1 and Bb < 0.: #if it's below the Curie temperature with B<0 then the middle minimum will be z2_left
            z1_left = 1. - z2_left - z2_right # constraint
              
            # system of differential equations (rate of change of phase fractions)
            dz2_left = -p_minus_mid*z2_left + z1_left*p_mid_minus
            dz2_right = -p_plus_mid*z2_right + z1_left*p_mid_plus
            dz1_mid = 0.
            dz1_right = 0.
            dz1_left = - dz2_left - dz2_right
            
            return [dz2_left, dz1_left, dz1_mid, dz1_right, dz2_right]

    
    elif len(min_indices) == 1: # if there is only 1 minimum
        # if there is only one minimum there are no transitions between minimums
        return np.zeros_like(z)

    elif len(min_indices) == 2: # if there are 2 minimums     
        max_index = argrelmax(fm)[0] # index of maximum
        
        C2 = np.exp(-fm[max_index]/(k_B*Tt)) # maximum
        C3_left, C3_right = norm_factor_sums # normalizing factors
        
        p_left_right = C2/C3_left # probability of a "particle" in the left minimum to be on top of the barrier
        p_right_left = C2/C3_right # probability of a "particle" in the right minimum to be on top of the barrier
            
        intersections = free_energy_intersections(sig, Tt, Bb, f1 ,f2) # indices of intersections
        
        if intersections == []: # if the free energies don't intersect it means that the two minimums belong to the high magnetization states of the phase with highest Curie temperature
            z2_right = 1. - z2_left # constraint
            
            # system of differential equations (rate of change of phase fractions)
            dz2_left = -p_left_right*z2_left + p_right_left*z2_right
            dz2_right = - dz2_left
            dz1_left = 0.
            dz1_mid = 0.
            dz1_right = 0.
            
            return [dz2_left, dz1_left, dz1_mid, dz1_right, dz2_right]

        else:
            a, b = intersections
        
            if min_indices[0] < a and min_indices[1] > b: # if the minimums belong to the high magnetization states
                z2_right = 1. - z2_left # constraint
                
                # system of differential equations (rate of change of phase fractions)
                dz2_left = -p_left_right*z2_left + p_right_left*z2_right
                dz2_right = -dz2_left
                dz1_left = 0.
                dz1_mid = 0.
                dz1_right = 0.

                return [dz2_left, dz1_left, dz1_mid, dz1_right, dz2_right]
                
            elif (a < min_indices[0] < b) and min_indices[1] > b: # if the one minimum belongs to the low mag. state and the other to the positive high mag. state                
                z2_right = 1. - z1_mid # constraint
     
                # system of differential equations (rate of change of phase fractions)
                dz1_mid = -p_left_right*z1_mid + p_right_left*z2_right
                dz2_right = -dz1_mid
                dz2_left = 0.
                dz1_left = 0.
                dz1_right = 0.
                
                return [dz2_left, dz1_left, dz1_mid, dz1_right, dz2_right]
                
            elif min_indices[0] < a and (a < min_indices[1] < b): # if the one minimum belongs to the low mag. state and the other to the negative high mag. state
                z1_mid = 1. - z2_left # constraint
                
                # system of differential equations (rate of change of phase fractions)
                dz2_left = -p_left_right*z2_left + p_right_left*z1_mid
                dz1_mid = - dz2_left
                dz2_right = 0.
                dz1_left = 0.
                dz1_right = 0.
                
                return [dz2_left, dz1_left, dz1_mid, dz1_right, dz2_right]
    
    elif len(min_indices) == 4:
        z1_right = 1. - z2_left - z1_left - z2_right # constraint
        
        max_index_left, max_index_mid, max_index_right = argrelmax(fm)[0] # indices of maximums
        
        C2_left = np.exp(-fm[max_index_left]/(k_B*Tt)) # left maximum
        C2_mid = np.exp(-fm[max_index_mid]/(k_B*Tt)) # middle maximum
        C2_right = np.exp(-fm[max_index_right]/(k_B*Tt)) # right maximum
        
        C3_2minus, C3_1minus, C3_1plus, C3_2plus = norm_factor_sums # nomalizing factors
        
        p_2minus_1minus = C2_left/C3_2minus # probability of a "particle" in the left minimum of phase 2 to be on top of the left barrier
        p_1minus_2minus = C2_left/C3_1minus # probability of a "particle" in the left minimum of phase 1 to be on top of the left barrier
        p_1minus_1plus = C2_mid/C3_1minus # probability of a "particle" in the left minimum of phase 1 to be on top of the middle barrier
        p_1plus_1minus = C2_mid/C3_1plus # probability of a "particle" in the right minimum of phase 1 to be on top of the middle barrier
        p_1plus_2plus = C2_right/C3_1plus # probability of a "particle" in the right minimum of phase 1 to be on top of the right barrier
        p_2plus_1plus = C2_right/C3_2plus # probability of a "particle" in the right minimum of phase 2 to be on top of the right barrier

        # system of differential equations (rate of change of phase fractions)
        dz2_left = -p_2minus_1minus*z2_left + p_1minus_2minus*z1_left
        dz1_left = p_2minus_1minus*z2_left - p_1minus_2minus*z1_left - p_1minus_1plus*z1_left + p_1plus_1minus*z1_right
        dz2_right = p_1plus_2plus*z1_right - p_2plus_1plus*z2_right
        dz1_right =  -dz2_left - dz1_left - dz2_right
        dz1_mid = 0.
        
        return [dz2_left, dz1_left, dz1_mid, dz1_right, dz2_right]

    
    

def frac_and_M(t):
    """Computes the phase fractions and the magnetization.
    
    Parameters
    ----------
    t : array
        Time.
    
    Returns
    -------
    y : 2D array
        Phase fractions and magnetization.
    """
    # Set initial condition automatically, z0
    fm = free.F_totstable_vs_M(sig, T_applied(0, T0, w, t0), B_applied(0, B0, w, t0)) # free energy function
    min_indices = argrelmin(fm)[0] # indices of minimums
    if len(min_indices) == 2: # if there are 2 minimums put 50% in each minimum
        intersections = free_energy_intersections(sig, T0, B0) # check if there are intersections
        if intersections == []: # if the curves of free energy do not intersect then the 2 minimums belong to the phase with highest Curie temperature (phase 2)
            z0 = [0.5, 0., 0., 0., 0.5]
        else:
            a, b = intersections # if there are intersections
            if min_indices[0] < a and min_indices[1] > b: # if both minimus belong to the positive and negative ferromagnetic state of phase 2
                z0 = [0.5, 0., 0., 0., 0.5]
            elif (a < min_indices[0] < b) and min_indices[1] > b: # if one minimum is from the paramagnetic state of phase 1 and the other from the positive ferromagnetic state of phase 2
                z0 = [0., 0., 0.5, 0., 0.5]
            elif min_indices[0] < a and (a < min_indices[1] < b): # if one minimum is from the paramagnetic state of phase 1 and the other from the negative ferromagnetic state of phase 2
                z0 = [0.5, 0., 0.5, 0., 0.]
    elif len(min_indices) == 4: # if there are 4 minimums then they are all from positive and negative ferromagnetic states of both structures
        z0 = [0.25, 0.25, 0., 0.25, 0.25]
    elif len(min_indices) == 3: # if there are 3 minimums then, at B=0, they belong to the positive and negative ferromagnetic states of phase 2 and to the paramagnetic state of phase 1
        z0 = [0.25, 0., 0.5, 0., 0.25]
    else: # of there is only 1 minimum the the entire system is in the middle minimum assuming B = 0 T at the start
        z0 = [0., 0., 1., 0., 0.]
    
    z = odeint(diff_eq, z0, t, hmax=1., mxstep=5000000) # hmax in default is faster but wrong results may appear
    
    z2_left = z[:,0] # unpack values of the negative ferromagnetic state of phase 2
    z1_left = z[:,1] # unpack values of the negative ferromagnetic state of phase 1
    z1_mid = z[:,2] # unpack values of the paramagnetic state of phase 1
    z1_right = z[:,3] # unpack values of the positive ferromagnetic state of phase 1
    z2_right = z[:,4] # unpack values of the positive ferromagnetic state of phase 2

    B = B_applied(t, B0, w, t0) # applied magnetic field
    T = T_applied(t, T0, w, t0) # applied temperature

    M = np.zeros_like(t) # create variable that for the magnetization values

    for i in range(len(t)): # for every point in time
        fm = free.F_totstable_vs_M(sig, T[i], B[i]) # total stable free energy

        min_indices = argrelmin(fm)[0] # indices of minimums
        if len(min_indices) == 1: # if there is only 1 minimum
            intersections = free_energy_intersections(sig, T[i], B[i]) # indeces of intersections
            
            if T[i] >= TC2: # if the temperature if above the highest Curie temperature then the minimum is from the paramagnetic state of phase 1
                sig1_mid = sig[argrelmin(fm)[0]] # value of the magnetization at the given temperature and magnetic field
                M[i] = sig1_mid*z1_mid[i] # magnetization at index i
            elif T[i] < TC2 and B[i] > 0.: # if the temperature if below the highest Curie temperature and B>0 then the minimum is from the positive ferromagnetic state of phase 2
                sig2_right = sig[argrelmin(fm)[0]] # value of the magnetization at the given temperature and magnetic field
                M[i] = sig2_right*z2_right[i] # magnetization at index i
            elif T[i] < TC2 and B[i] < 0.: # if the temperature if below the highest Curie temperature and B<0 then the minimum is from the negative ferromagnetic state of phase 2
                sig2_left = sig[argrelmin(fm)[0]] # value of the magnetization at the given temperature and magnetic field
                M[i] = sig2_left*z2_left[i] # magnetization at index i
                
        elif len(min_indices) == 2: # if there are 2 minimums
            intersections = free_energy_intersections(sig, T[i], B[i]) # indices of intersections
            
            if intersections == []: # if there are nointersections
                sig2_left, sig2_right = sig[argrelmin(fm)[0]] # magnetizations of both minimums at the given T and B
                M[i] = sig2_right*z2_right[i] + sig2_left*z2_left[i] # magnetization at index i
                
            else:
                a, b = intersections
                if min_indices[0] < a and min_indices[1] > b: # if both minimus belong to the positive and negative ferromagnetic state of phase 2
                    sig2_left, sig2_right = sig[argrelmin(fm)[0]] # magnetizations of both minimums at the given T and B
                    M[i] = sig2_right*z2_right[i] + sig2_left*z2_left[i] # magnetization at index i
    
                elif (a < min_indices[0] < b) and min_indices[1] > b:
                    sig1_mid, sig2_right = sig[argrelmin(fm)[0]] # magnetizations of both minimums at the given T and B
                    M[i] = sig1_mid*z1_mid[i] + sig2_right*z2_right[i] # magnetization at index i
    
                elif min_indices[0] < a and (a < min_indices[1] < b):
                    sig2_left, sig1_mid = sig[argrelmin(fm)[0]] # magnetizations of both minimums at the given T and B
                    M[i] = sig1_mid*z1_mid[i] + sig2_left*z2_left[i] # magnetization at index i
        
        elif len(min_indices) == 3:
            if T[i] >= TC1: #if it's above the Curie temperature then it will only be z1_mid
                sig2_left, sig1_mid, sig2_right = sig[argrelmin(fm)[0]] # magnetizations of minimums at the given T and B
                M[i] = sig1_mid*z1_mid[i] + sig2_right*z2_right[i] + sig2_left*z2_left[i] # magnetization at index i
            
            elif T[i] < TC1 and B[i] > 0.: # if the temperature if below the highest Curie temperature and B>0 then the minimum is from the positive ferromagnetic state of phase 2
                sig2_left, sig1_right, sig2_right = sig[argrelmin(fm)[0]] # magnetizations of minimums at the given T and B
                M[i] = sig1_right*z1_right[i] + sig2_right*z2_right[i] + sig2_left*z2_left[i] # magnetization at index i
            
            elif T[i] < TC1 and B[i] < 0.: # if the temperature if below the highest Curie temperature and B<0 then the minimum is from the negative ferromagnetic state of phase 2
                sig2_left, sig1_left, sig2_right = sig[argrelmin(fm)[0]] # magnetizations of minimums at the given T and B
                M[i] = sig1_left*z1_left[i] + sig2_right*z2_right[i] + sig2_left*z2_left[i] # magnetization at index i
        
        elif len(min_indices) == 4: # if there are 4 minimums
            sig2_left, sig1_left, sig1_right, sig2_right = sig[argrelmin(fm)[0]] # magnetizations of minimums at the given T and B
            M[i] = sig2_left*z2_left[i] + sig1_left*z1_left[i] + sig1_right*z1_right[i] + sig2_right*z2_right[i] # magnetization at index i
        
        
    return z2_left, z1_left, z1_mid, z1_right, z2_right, M # return fractions and magnetization


#==============================================================================
# Plots
#==============================================================================


if __name__ == "__main__":
    do_prof = 0
    do_plots = 1
    do_hys = 1
    save = 1
    
    if do_prof == True:
        print '\t Profiling...'
        import cProfile
        import Profiling as Prof
        import os
        file_name = "Profile_Output\ " + str(os.path.basename(__file__)) + '_Profile_Output'
        cProfile.runctx('frac_and_M(Delta_t)', {'frac_and_M':frac_and_M,'Delta_t':Delta_t}, {}, filename = file_name)
        Prof.save_profiling(file_name, sort='cumtime')
    
    if do_plots == True:
        print 'Plotting...'
        
        fm = free.F_totstable_vs_M(sig, T0, B0)
    
        print '\nMinimums\n',argrelmin(fm)[0], sig[argrelmin(fm)],argrelmax(fm), sig[argrelmax(fm)]
        plt.figure()
        plt.plot(sig, fm, label='T='+str(T0)+'K, B='+str(B0)+'T')
        plt.plot(sig, free.F_totstable_vs_M(sig, T0, 0.), label='T='+str(T0)+'K, B= 0T')
        
        plt.legend(loc=0, fontsize='small')
        plt.xlabel('$\sigma$')
        plt.ylabel('$F^{Tot}$ (eV)')
        plt.xlim(-1.,1.)
        
        b_applied = B_applied(Delta_t, B0,  w, t0)
        t_applied = T_applied(Delta_t, T0,  w, t0)
        
        plt.figure()
        plt.plot(Delta_t, b_applied)
        plt.title('B_applied')
        plt.ylabel('B (T)')
        plt.xlabel('t')
        
        
        roots = free_energy_intersections(sig, T0, B0)
        print '\nFree Energy Intersections\n', roots, sig[roots]
        
        norm = norm_factor(sig, T0, B0)
        print '\nNormalization Factors\n', norm#, norm[0]+norm[1]+norm[2]#+norm[3]#-norm[4]
        
        
        
        if do_hys == True:
            print '\nComputing phase fractions (Diff. Eq.) ... '   
            y = frac_and_M(Delta_t)
            
            plt.figure()
            plt.plot(Delta_t, y[0], label='z2_left')
            plt.plot(Delta_t, y[1], label='z1_left')
            plt.plot(Delta_t, y[2], label='z1_mid')
            plt.plot(Delta_t, y[3], label='z1_right')
            plt.plot(Delta_t, y[4], label='z2_right')
            plt.title('Fractions')
            plt.ylim(-0.05,1.05)
            plt.legend(loc=0, fontsize='small')
            plt.ylabel('Fractions')
            plt.xlabel('t')
            
            
            plt.figure()
            plt.plot(Delta_t, y[5], label='M(t)')
            plt.title('Magnetization')
            plt.legend(loc=0, fontsize='small')
            plt.ylabel('$\sigma$')
            plt.xlabel('t')
            
            plt.figure()
            plt.plot(B_applied(Delta_t, B0,  w, t0), y[5], label='Hysteresis')
            plt.title('Hysteresis')
            #plt.ylim(-0.05,1.05)
            plt.legend(loc=0, fontsize='small')
            plt.ylabel('$\sigma$')
            plt.xlabel('B (T)')


            if save == True:
                import os
                if not os.path.exists("Hysteresis"):
                    os.makedirs("Hysteresis")
                zipped = zip(Delta_t, t_applied, b_applied, y[0], y[1], y[2], y[3], y[4], y[5])
                np.savetxt("Hysteresis\Hysteresis(T="+str(T0)+"K, B="+str(B0)+"T).txt", zipped, delimiter=',', header='Time, Temperature, Magnetic Field, z2_left, z1_left, z1_mid, z1_right, z2_right, Reduced Magnetization',)

