# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 17:04:06 2016

@author: Rui M. Costa
"""

from scipy.integrate import quad

from magcalc import free_energy as free, magnetization as mag

mu_B = 5.7883818066*(10**(-5))  # eV T^-1
k_B = 8.6173324*(10**(-5))  # eV K^-1


def S_M(T, B, J, gJ, TC, lamb, Nm):
    """Computes the magnetic entropy for 1 spin.

    Parameters
    ----------
    T : scalar, 2D array
        Temperatures.
    B : scalar, 2D array
        Magnetic fields.
    J : scalar
        Angular momentum.
    TC : scalar
        Curie temperature.
    lamb : scalar
        Value of the strength of the parameter of the Molecular Field.

    Returns
    -------
    y : scalar, 2D array
        Magnetic entropy.
    """
    h = B/(lamb*Nm*gJ*mu_B*J)  # relative field to the saturation magnetization

    sigma = mag.Brillouin(T, B, J, gJ, TC, lamb, Nm)  # reduced magnetization

    y = 3.*J/(J+1.)*(h + sigma)*TC/T  # temporary variable

    A = np.sinh((2.*J + 1.)*y/(2.*J))
    B = np.sinh(y/(2.*J))

    return k_B*(np.log(A/B) - sigma*y)


def S_L(T, theta_D):
    """Computes the lattice entropy for 1 atom.

    Parameters
    ----------
    T : array
        Temperatures.
    theta_D : scalar
        Debye temperature of the system.

    Returns
    -------
    y : Lattice entropy.
    """

    # Function in the integral
    def f(x):
        return (x**3.)/(np.exp(x) - 1.)

    integral = np.zeros_like(T)  # variable that stores the values

    if np.isscalar(T):  # if T is just a single Temperature
        integral = quad(f, 0., theta_D/T)[0]
    elif len(T.shape) == 1:  # if T is an array of Temperatures
        # calculate the integral for each temperature
        for i, t in enumerate(T):
            integral[i] = quad(f, 0., theta_D/t)[0]
    else:
        raise Exception('Variable T should be a scalar or a 1D array.')

    return k_B*(12.*((T/theta_D)**3.)*integral -
                3.*np.log(1. - np.exp(-theta_D/T)))


def S_tot(T, B, J1, J2, gJ, TC1, TC2, theta_D1, theta_D2, F01, F02,
          lamb1, lamb2, N, Nm):
    """Computes the entropy following the stable phase.

    Parameters
    ----------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.

    Returns
    -------
    y : 2D array
        Entropy of stable phase.
    """
    # total free energy of phase 1
    F1 = free.F(T, B, J1, gJ, TC1, lamb1, Nm, theta_D1, N, F01)
    # total free energy of phase 2
    F2 = free.F(T, B, J2, gJ, TC2, lamb2, Nm, theta_D2, N, F02)

    F_cross_index2 = (F1 > F2).astype(int)  # determines index where F1 > F2
    F_cross_index1 = (F1 < F2).astype(int)  # determines index where F1 < F2

    # total entropy of phase 1
    s_tot1 = Nm*S_M(T, B, J1, gJ, TC1, lamb1, Nm) + \
        N*S_L(T[0], theta_D1)*np.ones(np.shape(T))

    # total entropy of phase 2
    s_tot2 = Nm*S_M(T, B, J2, gJ, TC2, lamb2, Nm) + \
        N*S_L(T[0], theta_D2)*np.ones(np.shape(T))

    # entropy of stable phase
    s_tot = F_cross_index1*s_tot1 + F_cross_index2*s_tot2

    return s_tot


def S_M_tot(T, B, J1, J2, gJ, TC1, TC2, theta_D1, theta_D2, F01, F02, lamb1,
            lamb2, N, Nm):
    """Computes the magnetic entropy following the stable phase.

    Parameters
    ----------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.

    Returns
    -------
    y : 2D array
        Entropy of stable phase.
    """
    # total free energy of phase 1
    F1 = free.F(T, B, J1, gJ, TC1, lamb1, Nm, theta_D1, N, F01)
    # total free energy of phase 2
    F2 = free.F(T, B, J2, gJ, TC2, lamb2, Nm, theta_D2, N, F02)

    F_cross_index2 = (F1 > F2).astype(int)  # determines index where F1 > F2
    F_cross_index1 = (F1 < F2).astype(int)  # determines index where F1 < F2

    # magnetic entropy of phase 1
    s_tot1 = Nm*S_M(T, B, J1, gJ, TC1, lamb1, Nm)
    # magnetic entropy of phase 2
    s_tot2 = Nm*S_M(T, B, J2, gJ, TC2, lamb2, Nm)

    # magnetic entropy of stable phase
    s_tot = F_cross_index1*s_tot1 + F_cross_index2*s_tot2

    return s_tot


def S_L_tot(T, B, J1, J2, gJ, TC1, TC2, theta_D1, theta_D2, F01, F02, lamb1,
            lamb2, N, Nm):
    """Computes the lattice entropy following the stable phase.

    Parameters
    ----------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.

    Returns
    -------
    y : 2D array
        Entropy of stable phase.
    """
    # total free energy of phase 1
    F1 = free.F(T, B, J1, gJ, TC1, lamb1, Nm, theta_D1, N, F01)
    # total free energy of phase 2
    F2 = free.F(T, B, J2, gJ, TC2, lamb2, Nm, theta_D2, N, F02)

    F_cross_index2 = (F1 > F2).astype(int)  # determines index where F1 > F2
    F_cross_index1 = (F1 < F2).astype(int)  # determines index where F1 < F2

    s_tot1 = N*S_L(T[0], theta_D1)  # lattice entropy of phase 1
    s_tot2 = N*S_L(T[0], theta_D2)  # lattice entropy of phase 2

    s_tot = (F_cross_index1*s_tot1+F_cross_index2*s_tot2)*np.ones(np.shape(T))

    return s_tot


def S_M_vs_M(sigma, T, B, J, gJ, TC, lamb, Nm):
    """Computes the magnetic entropy as a function of magnetization for a given
    temperature and magnetic field.

    Parameters
    ----------
    sigma : scalar, 1D array
        Reduced magnetization.
    T : scalar
        Temperatures.
    B : scalar
        Magnetic fields.
    J : scalar
        Total angular momentum.
    TC : scalar
        Curie temperature
    lamb : scalar
        Value of the strength of the parameter of the Molecular Field.

    Returns
    -------
    y : array
        Magnetic entropy as a function of magnetization.
    """
    Ms = Nm*gJ*mu_B*J

    A = 3./2.*(2.*J + 1.)/(J + 1.)*TC/T*(B/(lamb*Ms) + sigma)
    B = 3./2./(J + 1.)*TC/T*(B/(lamb*Ms) + sigma)

    C1 = -1./(2.*J*T)*(B/(lamb*Ms) + sigma)*((2.*J + 1.)/np.tanh(A) -
                                             1./np.tanh(B))
    C2 = -(J + 1.)/(3.*J*TC)*np.log(np.sinh(B)/np.sinh(A))

    return C1 + C2


if __name__ == "__main__":
    from magcalc.variables import *

    # print S_M(5., 0., J1, gJ, TC1, lamb1, Nm)
    # print S_M(TT, BB, J1, gJ, TC1, lamb1, Nm)

    # print S_L(5., theta_D1)
    # print S_L(TT[0], theta_D1)

    # print S_tot(TT, BB, J1, J2, gJ, TC1, TC2, theta_D1, theta_D2, F01, F02,
    #             lamb1, lamb2, N, Nm)

    # print S_M_tot(TT, BB, J1, J2, gJ, TC1, TC2, theta_D1, theta_D2, F01, F02,
    #               lamb1, lamb2, N, Nm)

    # print S_L_tot(TT, BB, J1, J2, gJ, TC1, TC2, theta_D1, theta_D2, F01, F02,
    #               lamb1, lamb2, N, Nm)

    # print S_M_vs_M(1., 5., 0., J1, gJ, TC1, lamb1, Nm)
    # print S_M_vs_M(sig, 5., 0., J1, gJ, TC1, lamb1, Nm)
