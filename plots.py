# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 17:43:15 2016

@author: Rui M. Costa
"""

import numpy as np
import matplotlib.pyplot as plt
import os

import magnetization as mag
import entropy as ent
import free_energy as free
import energy as ener


# Boltzmann Constant
k_B = 8.6173324 * (10**(-5))  # eV K^-1
# Bohr magneton
mu_B = 5.7883818066 * (10**(-5))  # eV T^-1


def plt_M_vs_T(T, B, sigma_1, sigma_2, sigma_stable, Bstep=1, save=False):
    """Plots the magnetization of both structures as a function of temperature.

    Parameters
    ----------
    T : array
        Array with the temperatures.
    B : scalar, array
        Magnetic fields.
    sigma_1, sigma_2, sigma_stable : array, 2D array
        Arrays with the values of the respective magnetizations
    Bstep : int
        Index step when plotting for several magnetic fields.
    save : bool
        Save the data of this plot to a .txt file.
    """

    print '\t Brillouin as a function of Temperature'

    plt.figure()

    if B.shape == ():  # if B is a scalar
        plt.plot(T, sigma_1, label='1, B=' + str(B) + 'T')
        plt.plot(T, sigma_2, label='2, B=' + str(B) + 'T')
        plt.plot(T, sigma_stable, label='stable, B=' + str(B) + 'T')

        if save:
            if not os.path.exists("M_vs_T"):
                os.makedirs("M_vs_T")
            zipped = zip(T, sigma_1, sigma_2, sigma_stable)
            np.savetxt("M_vs_T\M_vs_T(" + str(B) + "T).txt",
                       zipped,
                       delimiter=',',
                       header='Temperature, Magnetization 1, Magnetization 2, \
                               Stable Magnetization',
                       )

    else:  # if B is 2D
        for i in range(0, len(B), Bstep):
            plt.plot(T, sigma_1[i], ':', label='1, B=' + str(B[i]) + 'T')
        plt.gca().set_color_cycle(None)
        for i in range(0, len(B), Bstep):
            plt.plot(T, sigma_2[i], '--', label='2, B=' + str(B[i]) + 'T')
        plt.gca().set_color_cycle(None)
        for i in range(0, len(B), Bstep):
            plt.plot(T, sigma_stable[i], label='stable, B=' + str(B[i]) + 'T')

            if save:
                if not os.path.exists("M_vs_T"):
                    os.makedirs("M_vs_T")
                zipped = zip(T, sigma_1[i], sigma_2[i], sigma_stable[i])
                np.savetxt("M_vs_T\M_vs_T(" + str(B[i]) + "T).txt",
                           zipped,
                           delimiter=',',
                           header='Temperature, Magnetization 1, \
                           Magnetization 2, Stable Magnetization',
                           )

    plt.title('Reduced Magnetizations, $\sigma$')
    plt.legend(loc=0, fontsize='small')
    plt.ylim(-1.05, 1.05)
    plt.ylabel('$\sigma(T, B)$')
    plt.xlabel('T (K)')


def plt_M_vs_B(T, B, sigma_1, sigma_2, sigma_stable, Tstep=1, save=False):
    """Plots the magnetization of both structures as a function of temperature.

    Parameters
    ----------
    T : scalar, array
        Array with the temperatures.
    B : array
        Magnetic fields.
    sigma_1, sigma_2, sigma_stable : array, 2D array
        Arrays with the values of the respective magnetizations
    Tstep : int
        Index step when plotting for several temperatures.
    save : bool
        Save the data of this plot to a .txt file.
    """

    print '\t Brillouin as a function of Magnetic Field'

    plt.figure()

    if T.shape == ():  # if T is a scalar
        plt.plot(B, sigma_1, label='1, T=' + str(T) + 'K')
        plt.plot(B, sigma_2, label='2, T=' + str(T) + 'K')
        plt.plot(B, sigma_stable, label='stable, T=' + str(T) + 'K')

        if save:
            if not os.path.exists("M_vs_B"):
                os.makedirs("M_vs_B")
            zipped = zip(B, sigma_1, sigma_2, sigma_stable)
            np.savetxt("M_vs_B\M_vs_B(" + str(T) + "K).txt",
                       zipped,
                       delimiter=',',
                       header='Magnetic Field, Magnetization 1, \
                       Magnetization 2, Stable Magnetization',
                       )

    else:  # if T is 2D
        for i in range(0, len(T), Tstep):
            plt.plot(B, sigma_1[:, i], ':', label='1, T=' + str(T[i]) + 'K')
        plt.gca().set_color_cycle(None)
        for i in range(0, len(T), Tstep):
            plt.plot(B, sigma_2[:, i], '--', label='2, T=' + str(T[i]) + 'K')
        plt.gca().set_color_cycle(None)
        for i in range(0, len(T), Tstep):
            plt.plot(B, sigma_stable[:, i],
                     label='stable, T=' + str(T[i]) + 'K')

            if save:
                if not os.path.exists("M_vs_B"):
                    os.makedirs("M_vs_B")
                zipped = zip(B, sigma_1[:, i], sigma_2[
                             :, i], sigma_stable[:, i])
                np.savetxt("M_vs_B\M_vs_B(" + str(T[i]) + "K).txt",
                           zipped,
                           delimiter=',',
                           header='Magnetic Field, Magnetization 1, \
                           Magnetization 2, Stable Magnetization',
                           )

    plt.title('Reduced Magnetizations, $\sigma$')
    plt.legend(loc=0, fontsize='small')
    plt.ylim(0, 1.05)
    plt.ylabel('$\sigma(T, B)$')
    plt.xlabel('B (T)')


def plt_M_vs_TB(T, B, sigma_stable,):
    """Plots the stable magnetization as a function of temperature and
    magnetic field.

    Parameters
    ----------
    T : array
        Array with the temperatures.
    B : array
        Magnetic fields.
    sigma_stable : 2D array
        Arrays with the values of the respective magnetizations
    """
    print '\t 2D plot of Brillouin_stable'

    plt.figure()
    plt.imshow(sigma_stable, aspect='auto', extent=(T[0], T[-1], B[-1], B[0]))
    plt.title('Stable Reduced Magnetizations, $\sigma$')
    plt.colorbar()
    plt.ylabel('B (T)')
    plt.xlabel('T (K)')


def plt_U_vs_T(T, B, sigma_1, sigma_2, Bstep=1, save=False,
               BB=None, J1=None, J2=None, gJ=None, TC1=None, TC2=None):
    """Plots the internal energy of both structures as a function of
    temperature.

    Parameters
    ----------
    T : array
        Array with the temperatures.
    B : scalar, array
        Magnetic fields.
    sigma_1, sigma_2 : 2D array
        Arrays with the values of the respective magnetizations
    """
    print '\t Internal Magnetic Energy'

    print type(J1), type(BB), type(TC1)

    u_M1 = -gJ * mu_B * J1 * BB * sigma_1 - 3. * \
        J1 / (J1 + 1.) * k_B * TC1 * (sigma_1**2.)
    u_M2 = -gJ * mu_B * J2 * BB * sigma_2 - 3. * \
        J2 / (J2 + 1.) * k_B * TC2 * (sigma_2**2.)

    plt.figure()
    for i in range(0, len(B), Bstep):
        plt.plot(T, u_M1[i], label='1, B=' + str(B[i]) + 'T')
        plt.plot(T, u_M2[i], label='2, B=' + str(B[i]) + 'T')

        if save:
            if not os.path.exists("U_vs_T"):
                os.makedirs("U_vs_T")
            zipped = zip(T, u_M1[i], u_M2[i])
            np.savetxt("U_vs_T\U_vs_T(" + str(B[i]) + "T).txt",
                       zipped, delimiter=',',
                       header='Temperature, Internal Energy 1, \
                       Internal Energy 2',
                       )

    plt.title('Internal Magnetic Energy')
    plt.legend(loc=0, fontsize='small')
    plt.ylabel('$U^M(T)$')
    plt.xlabel('T (K)')


def plt_M_hys_vs_T(T, B, sigma_heat, sigma_cool, Bstep=1, save=False):
    """Plots the hysteresis of magnetization as a function of temperature.

    Parameters
    ----------
    T : array
        Array with the temperatures.
    B : scalar, array
        Magnetic fields.
    sigma_heat, sigma_cool : 2D array
        Arrays with the values of the respective magnetizations
    """
    print '\t Magnetization Hysteresis as a function of Temperature'

    plt.figure()
    for i in range(0, len(B), Bstep):
        plt.plot(T, sigma_heat[i], label='Heating, T=' + str(B[i]) + 'T')
    plt.gca().set_color_cycle(None)
    for i in range(0, len(B), Bstep):
        plt.plot(T, sigma_cool[i], '--', label='Cooling, T=' + str(B[i]) + 'T')

        if save:
            if not os.path.exists("M_hys_vs_T"):
                os.makedirs("M_hys_vs_T")
            zipped = zip(T, sigma_heat[i], sigma_cool[i])
            np.savetxt("M_hys_vs_T\M_hys_vs_T(" + str(B[i]) + "T).txt",
                       zipped, delimiter=',',
                       header='Temperature, Magnetization Heating, \
                       Magnetization Cooling',
                       )

    plt.title('Reduced Magnetizations, $\sigma$')
    plt.legend(loc=0, fontsize='small')
    plt.ylim(0, 1.05)
    plt.ylabel('$\sigma(T, B)$')
    plt.xlabel('T (K)')


def plt_M_hys_vs_B(T, B, sigma_heat, sigma_cool, Tstep=1, save=False):
    """Plots the hysteresis of magnetization as a function of magnetic field.

    Parameters
    ----------
    T : scalar, array
        Array with the temperatures.
    B : array
        Magnetic fields.
    sigma_heat, sigma_cool : 2D array
        Arrays with the values of the respective magnetizations
    """
    print '\t Magnetization Hysteresis as a function of Magnetic Field'

    plt.figure()
    for j in range(0, len(T), Tstep):
        plt.plot(B, sigma_heat[:, j], label='Heating, T=' + str(T[j]) + 'K')
    plt.gca().set_color_cycle(None)
    for j in range(0, len(T), Tstep):
        plt.plot(B, sigma_cool[:, j], '--',
                 label='Cooling, T=' + str(T[j]) + 'K')

        if save:
            print save
            if not os.path.exists("M_hys_vs_B"):
                os.makedirs("M_hys_vs_B")
            zipped = zip(B, sigma_heat[j], sigma_cool[j])
            np.savetxt("M_hys_vs_B\M_hys_vs_B(" + str(T[j]) + "T).txt",
                       zipped, delimiter=',',
                       header='Magnetic Field, Magnetization Heating, \
                       Magnetization Cooling',
                       )

    plt.legend(loc=0, fontsize='small')
    plt.title('Magnetization on Heating and Cooling')
    plt.xlabel('$B (T)$')
    plt.ylabel('$\sigma(T, B)$')


def plt_M_hys_vs_TB(T, B, sigma_heat, sigma_cool):
    """Plots the magnetization on heating and cooling as a function of
    temperature
    and magnetic field.

    Parameters
    ----------
    T : array
        Array with the temperatures.
    B : array
        Magnetic fields.
    sigma_heat, sigma_cool, sigma_stable : 2D array
        Arrays with the values of the respective magnetizations
    """

    plt.figure()
    plt.imshow(sigma_heat, aspect='auto', extent=(T[0], T[-1], B[-1], B[0]))
    plt.title('Reduced Magnetization Heating, $\sigma$')
    plt.colorbar()
    plt.ylabel('B (T)')
    plt.xlabel('T (K)')

    plt.figure()
    plt.imshow(sigma_cool, aspect='auto', extent=(T[0], T[-1], B[-1], B[0]))
    plt.title('Reduced Magnetization Cooling, $\sigma$')
    plt.colorbar()
    plt.ylabel('B (T)')
    plt.xlabel('T (K)')


def plt_S_M_vs_T(T, B, ent_1, ent_2, Bstep=1, save=False):
    """Plots the magnetic entropy of both structures.

    Parameters
    ----------
    T : array
        Temperatures.
    B : scalar, array
        Magnetic fields.
    ent_1, ent_2 : array
        Entropies of structure 1 and 2 respectively.
    Bstep : int
        Index step when plotting for several magnetic fields.
    save : bool
        Save the data of this plot to a .txt file.
    """
    print '\t Magnetic Entropy'

    plt.figure()
    for i in range(0, len(B), Bstep):
        plt.plot(T, ent_1[i], label='1, B=' + str(B[i]) + 'T')
        plt.plot(T, ent_2[i], label='2, B=' + str(B[i]) + 'T')

    plt.title('Magnetic Entropies, $S^M$')
    plt.ylim(0, np.amax(ent_1) + 0.00005)
    plt.legend(loc=0)
    plt.ylabel('$S^M(T, B)$')
    plt.xlabel('T (K)')


def plt_S_L_vs_T(T, B, ent_1, ent_2, Bstep=1, save=False):
    """Plots the lattice entropy of both structures.

    Parameters
    ----------
    T : array
        Temperatures.
    B : scalar, array
        Magnetic fields.
    ent_1, ent_2 : array
        Entropies of structure 1 and 2 respectively.
    Bstep : int
        Index step when plotting for several magnetic fields.
    save : bool
        Save the data of this plot to a .txt file.
    """
    print '\t Lattice Entropy'

    plt.figure()
    plt.plot(T, ent_1, label='1')
    plt.plot(T, ent_2, label='2')
    plt.ylim(0, np.amax(ent_1) + 0.0001)
    plt.title('Lattice Entropies, $S^L$')
    plt.legend(loc=0)
    plt.ylabel('$S^L(T)$')
    plt.xlabel('T (K)')


def plt_S_tot_vs_T(T, B, s_tot, Bstep=1, save=False):
    """Plots the entropy of stable phase.

    Parameters
    ----------
    T : array
        Temperatures.
    B : scalar, array
        Magnetic fields.
    s_tot : array
        Entropy of stable phase.
    Bstep : int
        Index step when plotting for several magnetic fields.
    save : bool
        Save the data of this plot to a .txt file.
    """
    print '\t Entropy of stable phase'

    plt.figure()
    for i in range(0, len(B), Bstep):
        plt.plot(T, s_tot[i], label='$S^{Tot}$, B=' + str(B[i]) + 'T')

    plt.legend(loc=0)
    plt.title('Entropy, S')
    plt.xlabel('T (K)')
    plt.ylabel('S(B) (eV/K)')


def plt_DeltaS_vs_T(T, B, s_tot, s_M_tot, s_L_tot, Bstep=1, save=False,
                    conv=False):
    """Plots the entropy changes of stable phase.

    Parameters
    ----------
    T : array
        Temperatures.
    B : scalar, array
        Magnetic fields.
    s_tot, s_M_tot, s_L_tot : array
        Entropy of total, magnetic and lattice stable phase.
    Bstep : int
        Index step when plotting for several magnetic fields.
    save : bool
        Save the data of this plot to a .txt file.
    conv = bool
        Enables the conversion from eV/K to J/Kg K (or the conversion used in \
        Variables.py).
    """
    if conv is False:  # if the conversion is disabled
        Conv = 1.

    plt.figure()
    for i in range(1, len(B), Bstep):
        plt.plot(T,
                 Conv * (s_tot[i] - s_tot[0]),
                 label='$\Delta S^{Tot}($' + str(B[0]) + '$\longrightarrow$' +
                 str(B[i]) + ' T$)$')

    plt.gca().set_color_cycle(None)
    for i in range(1, len(B), Bstep):
        plt.plot(T,
                 Conv * (s_M_tot[i] - s_M_tot[0]),
                 '--',
                 label='$\Delta S^{M}($' + str(B[0]) + '$\longrightarrow$' +
                 str(B[i]) + ' T$)$')

    plt.gca().set_color_cycle(None)
    for i in range(1, len(B), Bstep):
        plt.plot(T,
                 Conv * (s_L_tot[i] - s_L_tot[0]),
                 ':',
                 label='$\Delta S^{L}($' + str(B[0]) + '$\longrightarrow$' +
                 str(B[i]) + ' T$)$')

        if save:
            if not os.path.exists("DeltaS_vs_T"):
                os.makedirs("DeltaS_vs_T")
            zipped = zip(T,
                         Conv * (s_tot[i] - s_tot[0]),
                         Conv * (s_M_tot[i] - s_M_tot[0]),
                         Conv * (s_L_tot[i] - s_L_tot[0]))
            np.savetxt("DeltaS_vs_T\DeltaS_vs_T(" + str(B[i]) + "T).txt",
                       zipped,
                       delimiter=',',
                       header='Temperature, Total Entropy Change, \
                       Magnetic Entropy Change, Lattice Entropy Change',
                       )

    plt.title('Entropy Change, $\Delta S$')
    plt.legend(loc=0, fontsize='small')
    plt.xlabel('T (K)')
    if conv is False:
        plt.ylabel('$\Delta S(B)$ (eV/K)')
    else:
        plt.ylabel('$\Delta S(B)$ (J/Kg K)')


def plt_max_DeltaS_vs_B(B, s_tot, save=False, conv=False):
    """Plots the maximum entropy changes as function of magnetic field.

    Parameters
    ----------
    B : scalar, array
        Magnetic fields.
    s_tot : array
        Entropy of total stable phase.
    save : bool
        Save the data of this plot to a .txt file.
    conv = bool
        Enables the conversion from eV/K to J/Kg K (or the conversion used in
        Variables.py).
    """
    print '\t Maximum Entropy Change as a Function of Applied Magnetic Field'

    s_B = np.zeros(np.shape(B))

    plt.figure()
    for i in range(len(B)):
        s_B[i] = np.amax(np.abs(s_tot[i] - s_tot[0]))

    if conv is False:
        Conv = 1.

    plt.plot(B, Conv * s_B)  # plot after 0 T

    plt.title('Maximum Entropy Change, $\Delta S$')
    plt.xlabel('B (T)')
    plt.ylabel('$\Delta S^{Max}(B)$ (eV/K)')  # J/Kg K

    if save:
        if not os.path.exists("max_DeltaS_vs_B"):
            os.makedirs("max_DeltaS_vs_B")
        zipped = zip(B, Conv * s_B)
        np.savetxt("max_DeltaS_vs_B\max_DeltaS_vs_B.txt",
                   zipped,
                   delimiter=',',
                   header='magnetic Field, Maximum Total Entropy Change',
                   )


def plt_S_M_vs_M(sigma, T, B, J1, J2, gJ, TC1, TC2, lamb1, lamb2, Nm,
                 Tstep=1, Bstep=1, save=False):
    print '\t Maximum Entropy Change as a Function of Applied Magnetic Field'

    for j in range(0, len(B), Bstep):
        plt.figure()
        for i in range(0, len(T), Tstep):
            f1 = ent.S_M_vs_M(sigma, T[i], B[j], J1, gJ, TC1, lamb1, Nm)
            plt.plot(sigma, f1, label='$S^M_1$, B=' + str(T[i]) + 'K')
        plt.gca().set_color_cycle(None)
        for i in range(0, len(T), Tstep):
            f2 = ent.S_M_vs_M(sigma, T[i], B[j], J2, gJ, TC2, lamb2, Nm)
            plt.plot(sigma, f2, '--', label='2, T=' + str(T[i]) + 'K')

            if save:
                if not os.path.exists("S_M_vs_M"):
                    os.makedirs("S_M_vs_M")
                zipped = zip(sigma, f1, f2)
                np.savetxt("S_M_vs_M\S_M_vs_M(" + str(T[i]) + "K," +
                           str(B[j]) + "T).txt",
                           zipped,
                           delimiter=',',
                           header='Reduced Magnetization, Entropy 1, \
                           Entropy 2',
                           )

        plt.title('Magnetic Entropy, $S^M$(B=' + str(B[j]) + ')')
        plt.xlabel('$\sigma$')
        plt.ylabel('$S^M$ (eV/K)')  # J/Kg K
        plt.legend(loc=0, fontsize='small')


def plt_F_vs_T(T, B, free_ener_1, free_ener_2, free_stable, Bstep=1,
               save=False):
    """Plots the free energies of both structures as functions of temperature.

    Parameters
    ----------
    T : array
        Temperatures.
    B : scalar, array
        Magnetic fields.
    free_ener_1, free_ener_2, free_stable : 2D array
        Free energies of structures 1 and 2 respectively.
    Bstep : int
        Index step when plotting for several magnetic fields.
    save : bool
        Save the data of this plot to a .txt file.
    """
    print '\t Free Energy as a function of Temperature'

    plt.figure()
    for i in range(0, len(B), Bstep):
        plt.plot(T, free_ener_1[i], '--', label='1, B=' + str(B[i]) + 'T')
    plt.gca().set_color_cycle(None)
    for i in range(0, len(B), Bstep):
        plt.plot(T, free_ener_2[i], ':', label='2, B=' + str(B[i]) + 'T')
        # plt.plot(T, free_stable[i], label='stable, B='+str(B[i])+'T')

        if save:
            if not os.path.exists("F_vs_T"):
                os.makedirs("F_vs_T")
            zipped = zip(T,  free_ener_1[i], free_ener_2[i], free_stable[i])
            np.savetxt("F_vs_T\F_vs_T(" + str(B[i]) + "T).txt",
                       zipped,
                       delimiter=',',
                       header='Temperature, Free Energy 1, Free Energy 2, \
                       Free Enrgy Stable',
                       )

    plt.legend(loc=0)
    plt.title('Free Energies, $\Delta F(T)$')
    plt.xlabel('T (K)')
    plt.ylabel('F(T) (eV)')


def plt_F_vs_B(T, B, free_ener_1, free_ener_2, free_stable, Tstep=1,
               save=False):
    """Plots the free energies of both structures as functions of magnetic
    field.

    Parameters
    ----------
    T : array
        Temperatures.
    B : array
        Magnetic fields.
    free_ener_1, free_ener_2, free_stable : 2D array
        Free energies of structures 1 and 2 respectively.
    Bstep : int
        Index step when plotting for several magnetic fields.
    save : bool
        Save the data of this plot to a .txt file.
    """
    print '\t Free Energy as a function of Magnetic Field'

    plt.figure()
    for i in range(0, len(T), Tstep):
        plt.plot(B, free_ener_1[:, i], label='1, T=' + str(T[i]) + 'K')
        plt.plot(B, free_ener_2[:, i], label='2, T=' + str(T[i]) + 'K')
        plt.plot(B, free_stable[:, i], label='stable, T=' + str(T[i]) + 'K')

    plt.legend(loc=0)
    plt.title('F vs B')
    plt.xlabel('B (T)')
    plt.ylabel('F (eV)')


def plt_transition_temp(T, B, save, J1, TC1, theta_D1, F01, lamb1,
                        J2, TC2, theta_D2, F02, lamb2, gJ, Nm, N, *args):
    """Plots the transition temperatures as a function of magnetic field.

    Parameters
    ----------
    T : 2D array
            Temperature
    B : 2D array
            Applied Magnetic Field
    """
    print '\t Transition Temperatures'

    Ts = free.transition_temp(T, B, J1, TC1, theta_D1, F01, lamb1,
                              J2, TC2, theta_D2, F02, lamb2, gJ, Nm, N, *args)

    if B.shape[0] == 1:  # if there is only value for the magnetic field
        print Ts

    else:  # if there is an array for the magnetic field
        plt.figure()
        plt.plot(B[:, 0], Ts)

        plt.legend(loc=0)
        plt.title('Transition temperatures')
        plt.xlabel('B (T)')
        plt.ylabel('$T_S$ (K)')

        if save:
            if not os.path.exists("Ts_vs_B"):
                os.makedirs("Ts_vs_B")
            zipped = zip(B[:, 0], Ts)
            np.savetxt("Ts_vs_B\Ts_vs_B.txt", zipped, delimiter=',',
                       header='Magnetic Field, Transition Temperatures',)


def plt_F_M_vs_M(sigma, T, B, J1, TC1, lamb1, J2, TC2, lamb2, gJ, Nm,
                 Tstep=1, Bstep=1, save=False):
    """Plots the magnetic free energy as a function of the reduced
    magnetization.

    Parameters
    ---------
    sigma : scalar, array
        Reduced magnetization.
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.
    Tstep : int
        Index step when plotting for several temperatures.
    Bstep : int
        Index step when plotting for several magnetic fields.
    save : bool
        Save the data of this plot to a .txt file.
    """
    print '\t Magnetic Free Energy as a function of Magnetization'

    Delta_T = T[0]
    Delta_B = B[:, 0]

    for j in range(0, len(Delta_B), Bstep):
        plt.figure()
        for i in range(0, len(Delta_T), Tstep):
            f1 = free.F_M_vs_M(sigma, Delta_T[i], Delta_B[
                               j], J1, gJ, TC1, lamb1, Nm)
            plt.plot(sigma, f1, ':', label='1, T=' + str(Delta_T[i]) + 'K')
        plt.gca().set_color_cycle(None)
        for i in range(0, len(Delta_T), Tstep):
            f2 = free.F_M_vs_M(sigma, Delta_T[i], Delta_B[
                               j], J2, gJ, TC2, lamb2, Nm)
            plt.plot(sigma, f2, '--', label='2, T=' + str(Delta_T[i]) + 'K')

        plt.legend(loc=0, fontsize='small')
        plt.title('Magnetic Free Energy as a function of Magnetization \
                  $F^M$, B=' + str(Delta_B[j]) + 'T')
        plt.xlabel('$\sigma$')
        plt.ylabel('F(T) (eV)')
        plt.xlim(-1, 1)


def plt_F_M_vs_T(T, B, J1, J2, TC1, TC2, lamb1, lamb2, gJ, Nm, Bstep=1,
                 save=False):
    """Plots the free energies of both structures as functions of temperature.

    Parameters
    ----------
    T : array
        Temperatures.
    B : scalar, array
        Magnetic fields.
    free_ener_1, free_ener_2, free_stable : 2D array
        Free energies of structures 1 and 2 respectively.
    Bstep : int
        Index step when plotting for several magnetic fields.
    save : bool
        Save the data of this plot to a .txt file.
    """
    print '\t Free Energy as a function of Temperature'

    free_ener_1 = free.F_M(T, B, J1, gJ, TC1, lamb1, Nm)
    free_ener_2 = free.F_M(T, B, J2, gJ, TC2, lamb2, Nm)

    Delta_T = T[0]
    Delta_B = B[:, 0]

    plt.figure()
    for i in range(0, len(Delta_B), Bstep):
        plt.plot(Delta_T, free_ener_1[i],
                 label='1, B=' + str(Delta_B[i]) + 'T')
        plt.plot(Delta_T, free_ener_2[i],
                 label='2, B=' + str(Delta_B[i]) + 'T')

        if save:
            if not os.path.exists("F_M_vs_T"):
                os.makedirs("F_M_vs_T")
            zipped = zip(Delta_T,  free_ener_1[i], free_ener_2[i])
            np.savetxt("F_M_vs_T\F_M_vs_T(" + str(Delta_B[i]) + "T).txt",
                       zipped,
                       delimiter=',',
                       header='Temperature, Free Energy 1, Free Energy 2',
                       )

    plt.legend(loc=0)
    plt.title('Magentic Free Energies, $\Delta F(T)$')
    plt.xlabel('T (K)')
    plt.ylabel('F(T) (eV)')


def plt_F_L_vs_T(T, theta_D1, theta_D2):
    """Plots the lattice free energy as a function of temperature.

    Parameters
    ----------
    T : array
        Temperatures.
    """
    print '\t Lattice Free Energy as a function of Temperature'

    f_L1 = free.F_L(T, theta_D1)
    f_L2 = free.F_L(T, theta_D2)

    plt.figure()
    plt.plot(T, f_L1, label='1')
    plt.plot(T, f_L2, label='2')
    # plt.ylim(0, np.amax(f_L1)+0.0001)
    plt.title('Lattice Free Energy, $F^L$')
    plt.legend(loc=0)
    plt.ylabel('$F^L(T)$')
    plt.xlabel('T (K)')


def plt_Ftot_vs_M(sigma, T, B, J1, J2, TC1, TC2, lamb1, lamb2, theta_D1,
                  theta_D2, F01, F02, gJ, Nm, N, Tstep=1, Bstep=1,
                  save=False):
    """Plots the total free energy of both structures as functions of
    magnetization.

    Parameters
    ----------
    sigma : array
        Reduced magnetizations.
    T : array
        Temperatures.
    B : array
        Magnetic fields.
    Tstep : int
        Index step when plotting for several temperatures.
    Bstep : int
        Index step when plotting for several magnetic fields.
    save : bool
        Save the data of this plot to a .txt file.
    """
    print '\t Total Free Energy as a function of Magnetization'

    for j in range(0, len(B), Bstep):
        plt.figure()
        for i in range(0, len(T), Tstep):
            f1 = free.F_tot_vs_M(sigma, T[i], B[j], J1, gJ, TC1, lamb1, Nm,
                                 theta_D1, N, F01)
            plt.plot(sigma, f1, ':', label='1, T=' + str(T[i]) + 'K')

        plt.gca().set_color_cycle(None)
        for i in range(0, len(T), Tstep):
            f2 = free.F_tot_vs_M(sigma, T[i], B[j], J2, gJ, TC2, lamb2, Nm,
                                 theta_D2, N, F02)
            plt.plot(sigma, f2, '--', label='2, T=' + str(T[i]) + 'K')

        plt.gca().set_color_cycle(None)
        for i in range(0, len(T), Tstep):
            f3 = free.F_totstable_vs_M(sigma, T[i], B[j], J1, TC1, lamb1,
                                       theta_D1, F01, J2, TC2, lamb2, theta_D2,
                                       F02, gJ, Nm, N)
            plt.plot(sigma, f3, label='3, T=' + str(T[i]) + 'K')

            if save:
                if not os.path.exists("Ftot_vs_M"):
                    os.makedirs("Ftot_vs_M")
                zipped = zip(sigma, f1, f2, f3)
                np.savetxt("Ftot_vs_M\Ftot_vs_M(" + str(T[i]) + "K," +
                           str(B[j]) + "T).txt",
                           zipped, delimiter=',',
                           header='Reduced Magnetization, Free Energy 1, \
                           Free Energy 2, Stable Free Energy',
                           )

        plt.legend(loc=0, fontsize='small')
        plt.title('Total Free Energy as a function of Magnetization \
                  $F^{Tot}=F^M+F^L+F^E$, B=' + str(B[j]) + 'T')
        plt.xlabel('$\sigma$')
        plt.ylabel('F(T) (eV)')
        plt.xlim(-1, 1)


def plt_F_heatcool_vs_T(T, B, f_heat, f_cool, Bstep=1, save=False):
    """Plots the free energy of the metastable minimum on heating and cooling
    as functions of temperature.

    Parameters
    --------
    T : array
        Temperatures.
    B : array
        Magnetic fields.
    f_heat, f_cool : array
        Free energies on heating and cooling.
    Bstep : int
        Index step when plotting for several magnetic fields.
    save : bool
        Save the data of this plot to a .txt file.
    """
    print '\t Total Free Energy on Heating and Cooling'

    plt.figure()
    for j in range(0, len(B), Bstep):
        plt.plot(T, f_heat[j], label='Heating, B=' + str(B[j]) + 'T')
    plt.gca().set_color_cycle(None)
    for j in range(0, len(B), Bstep):
        plt.plot(T, f_cool[j], '--', label='Cooling, B=' + str(B[j]) + 'T')

    plt.legend(loc=0, fontsize='small')
    plt.title('Total Free Energy on Heating and Cooling $F^{Tot}=F^M+F^L$')
    plt.xlabel('$T (K)$')
    plt.ylabel('F(T) (eV)')


def plt_F_heatcool_vs_B(T, B, f_heat, f_cool, Tstep=1, save=False):
    """Plots the free energy of the metastable minimum on heating and cooling
    as functions of temperature.

    Parameters
    --------
    T : array
        Temperatures.
    B : array
        Magnetic fields.
    f_heat, f_cool : array
        Free energies on heating and cooling.
    Tstep : int
        Index step when plotting for several temepratures.
    save : bool
        Save the data of this plot to a .txt file.
    """
    plt.figure()
    for j in range(0, len(T), Tstep):  # (25,len(Delta_T)-165,1)
        plt.plot(B, f_heat[:, j], label='Heating, T=' + str(T[j]) + 'K')
    plt.gca().set_color_cycle(None)
    for j in range(0, len(T), Tstep):  # (25,len(Delta_T)-165,1)
        plt.plot(B, f_cool[:, j], '--', label='Cooling, T=' + str(T[j]) + 'K')

    plt.legend(loc=0, fontsize='small')
    plt.title('Total Free Energy on Heating and Cooling $F^{Tot}=F^M+F^L$')
    plt.xlabel('$B (T)$')
    plt.ylabel('F(T) (eV)')


def plt_E_M_vs_T(T, B, em1, em2, Bstep=1):
    """Plots the magnetic energy as a function of temperature.

    Parameters
    ----------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.
    em1, em2 : 2D array
        Magnetic energies of phases 1 and 2.
    Bstep : int
        Index step when plotting for several magnetic fields.
    """
    print '\t Magnetic Energy as a function of Temperature'

    plt.figure()
    for i in range(0, len(B), Bstep):
        plt.plot(T, em1[i, :], label='$E^M_1$, B=' + str(B[i]) + 'T')
        plt.plot(T, em2[i, :], label='$E^M_2$, B=' + str(B[i]) + 'T')

    plt.legend(loc=0, fontsize='small')
    plt.title('Magnetic Energy')
    plt.xlabel('T (K)')
    plt.ylabel('$E^M$ (eV)')


def plt_E_L_vs_T(T, el1, el2):
    """Plots the lattice energy as a function of temperature.

    Parameters
    ----------
    T : 2D array
        Temperatures.
    el1, el2 : 2D array
        Lattice energies of phases 1 and 2.
    """
    print '\t Lattice Energy as a function of Temperature'

    plt.figure()
    plt.plot(T, el1, label='$E^L_1$')
    plt.plot(T, el2, label='$E^L_2$')

    plt.legend(loc=0, fontsize='small')
    plt.title('Lattice Energy')
    plt.xlabel('T (K)')
    plt.ylabel('$E^L$ (eV)')


def plt_Etot_vs_T(T, B, etot1, etot2, Bstep=1):
    """Plots the total energies as a functions of temperature.

    Parameters
    ----------
    T : 2D array
        Temperatures.
    B : 2D array
        Magnetic fields.
    em1, em2 : 2D array
        Magnetic energies of phases 1 and 2.
    Bstep : int
        Index step when plotting for several magnetic fields.
    """
    print '\t Total Energy as a function of Temperature'

    plt.figure()
    for i in range(0, len(B), Bstep):
        plt.plot(T, etot1[i, :], label='$E^T_1$, B=' + str(B[i]) + 'T')
        plt.plot(T, etot2[i, :], label='$E^T_2$, B=' + str(B[i]) + 'T')

    plt.legend(loc=0, fontsize='small')
    plt.title('Total Energy $E^{Tot}=E^M+E^L+E^E$')
    plt.xlabel('T (K)')
    plt.ylabel('E (eV)')


def plt_E_M_vs_M(sigma, T, B, Tstep=1, Bstep=1, save=False):
    """Plots the total free energy of both structures as functions of
    magnetization.

    Parameters
    ----------
    sigma : array
        Reduced magnetizations.
    T : array
        Temperatures.
    B : array
        Magnetic fields.
    Tstep : int
        Index step when plotting for several temperatures.
    Bstep : int
        Index step when plotting for several magnetic fields.
    save : bool
        Save the data of this plot to a .txt file.
    """
    print '\t Magnetic Energy as a function of Magnetization'

    for j in range(0, len(B), Bstep):
        plt.figure()
        for i in range(0, len(T), Tstep):
            plt.plot(sigma, ener.E_M_vs_M(sigma, T[i], B[j], J1, TC1, lamb1)[
                     0], ':', label='1, T=' + str(T[i]) + 'K')

        plt.gca().set_color_cycle(None)
        for i in range(0, len(T), Tstep):
            plt.plot(sigma, ener.E_M_vs_M(sigma, T[i], B[j], J2, TC2, lamb2)[
                     0], label='2, T=' + str(T[i]) + 'K')

        plt.legend(loc=0, fontsize='small')
        plt.xlim(-1, 1)
        plt.title('Magnetic Energy $E^M$, B=' + str(B[j]) + 'T')
        plt.xlabel('$\sigma$')
        plt.ylabel('E (eV)')


def plt_M(MvsT=0, MvsB=0, MvsTB=0, UvsT=0, M_hys_vs_T=0, save=0, TT=None,
          BB=None, J1=None, J2=None, TC1=None, TC2=None, lamb1=None,
          lamb2=None, Delta_T=None, Delta_B=None, theta_D1=None, theta_D2=None,
          gJ=None, F01=None, F02=None, Nm=None, N=None):
    """Menu for the magnetization plots.

    Parameters
    ----------
    MvsT : bool
        Plots the magnetization vs temperature.
    MvsB : bool
        Plots the magnetization vs magnetic field.
    MvsTB : bool
        Plots the magnetization vs temperature and magnetic field.
    UvsT : bool
        Plots the internal energy vs temperature.
    M_hys_vs_T : bool
        Plots the hysteresis of the magnetization vs tempeperature, magnetic
        field and both.
    """
    print 'Magnetization'
    if MvsT or MvsB or MvsTB or UvsT:
        # magnetization of phase 1
        sig_1 = mag.Brillouin(TT, BB, J1, gJ, TC1, lamb1, Nm)
        # magnetization of phase 2
        sig_2 = mag.Brillouin(TT, BB, J2, gJ, TC2, lamb2, Nm)
        # magnetization of stable phase
        sig_stable = mag.Brillouin_stable(TT, BB, J1, J2, TC1, TC2, lamb1,
                                          lamb2, theta_D1, theta_D2, F01, F02,
                                          gJ, Nm, N)
        if MvsT:  # plot of M vs T
            plt_M_vs_T(Delta_T, Delta_B, sig_1, sig_2,
                       sig_stable, Bstep=1, save=save)

        if MvsB:  # plot of M vs B
            plt_M_vs_B(Delta_T, Delta_B, sig_1, sig_2,
                       sig_stable, Tstep=1, save=save)

        if MvsTB:  # plot of M vs T and B
            plt_M_vs_TB(Delta_T, Delta_B, sig_stable)

        if UvsT:  # plot of U vs T
            plt_U_vs_T(Delta_T, Delta_B, sig_1, sig_2, Bstep=1, save=save,
                       BB=BB, J1=J1, J2=J2, gJ=gJ, TC1=TC1, TC2=TC2)

    if M_hys_vs_T:
        # magnetization on heating
        sig_heat = mag.RedMag_heat(TT, BB, J1, TC1, lamb1, theta_D1, F01, J2,
                                   TC2, lamb2, theta_D2, F02, gJ, Nm, N)
        # magnetization on cooling
        sig_cool = mag.RedMag_cool(TT, BB, J1, TC1, lamb1, theta_D1, F01, J2,
                                   TC2, lamb2, theta_D2, F02, gJ, Nm, N)

        # plot of magnetic hysteresis vs T
        plt_M_hys_vs_T(Delta_T, Delta_B, sig_heat,
                       sig_cool, Bstep=1, save=save)
        # plot of magnetic hysteresus vs B
        plt_M_hys_vs_B(Delta_T, Delta_B, sig_heat,
                       sig_cool, Tstep=1, save=save)
        # plot of magnetic hysteresis vs T and B
        plt_M_hys_vs_TB(Delta_T, Delta_B, sig_heat, sig_cool)

    return None


def plt_S(S_M_vs_T=0, S_L_VS_T=0, S_tot_vs_T=0, DeltaS_vs_T=0,
          max_DeltaS_vs_B=0, S_M_vs_M=0, save=0, TT=None, BB=None, J1=None,
          J2=None, TC1=None, TC2=None, lamb1=None, lamb2=None,
          Delta_T=None, Delta_B=None, theta_D1=None, theta_D2=None, gJ=None,
          F01=None, F02=None, Nm=None, N=None, sigma=None):
    """Menu for the entropy plots.

    Parameters
    ----------
    S_M_vs_T : bool
        Plots the magnetic entropy vs temperature.
    S_L_vs_T : bool
        Plots the lattice entropy vs temperature.
    S_tot_vs_T : bool
        Plots the total entropy vs temperature.
    DeltaS_vs_T : bool
        Plots the entropy variation vs temperature.
    max_DeltaS_vs_B : bool
        Plots the maximum entropy vatiation vs magnetic field.
    S_M_vs_M : bool
        Plots the magnetic entropy as a function of the reduced magnetization.
    """
    print 'Entropy'

    if S_M_vs_T:
        s_M1 = ent.S_M(TT, BB, J1, gJ, TC1, lamb1, Nm)
        s_M2 = ent.S_M(TT, BB, J2, gJ, TC2, lamb2, Nm)

        plt_S_M_vs_T(Delta_T, Delta_B, s_M1, s_M2, Bstep=1, save=save)

    if S_L_VS_T:
        s_L1 = ent.S_L(Delta_T, theta_D1)
        s_L2 = ent.S_L(Delta_T, theta_D2)

        plt_S_L_vs_T(Delta_T, Delta_B, s_L1, s_L2, Bstep=1, save=save)

    if S_tot_vs_T or DeltaS_vs_T or max_DeltaS_vs_B:
        s_tot = ent.S_tot(TT, BB, J1, J2, gJ, TC1, TC2, theta_D1, theta_D2,
                          F01, F02, lamb1, lamb2, N, Nm)

        if S_tot_vs_T:
            plt_S_tot_vs_T(Delta_T, Delta_B, s_tot, Bstep=1, save=save)

        if DeltaS_vs_T:
            print '1'
            s_M_tot = ent.S_M_tot(TT, BB, J1, J2, gJ, TC1, TC2, theta_D1,
                                  theta_D2, F01, F02, lamb1, lamb2, N, Nm)
            print '2'
            s_L_tot = ent.S_L_tot(TT, BB, J1, J2, gJ, TC1, TC2, theta_D1,
                                  theta_D2, F01, F02, lamb1, lamb2, N, Nm)
            print '3'
            plt_DeltaS_vs_T(Delta_T, Delta_B, s_tot, s_M_tot,
                            s_L_tot, Bstep=1, save=save, conv=0)
            print '4'
        if max_DeltaS_vs_B:
            print save
            plt_max_DeltaS_vs_B(Delta_B, s_tot, save=save, conv=0)

    if S_M_vs_M:
        plt_S_M_vs_M(sigma, Delta_T, Delta_B, J1, J2, gJ, TC1, TC2,
                     lamb1, lamb2, Nm, Tstep=10, Bstep=1, save=save)


def plt_F(FvsT=0, FvsB=0, trans_temp=0, F_M_vs_T=0, F_M_vs_M=0, F_L_vs_T=0,
          FtotvsM=0, Ftot_heatcool=0, save=0, TT=None, BB=None, J1=None,
          J2=None, TC1=None, TC2=None, lamb1=None, lamb2=None, Delta_T=None,
          Delta_B=None, theta_D1=None, theta_D2=None, gJ=None, F01=None,
          F02=None, Nm=None, N=None, sigma=None):
    """Menu for the free energy plots.

    Parameters
    ----------
    FvsT : bool
        Plots the free energy vs temperature.
    FvsB : bool
        Plots the free energy vs magnetic field.
    trans_temp : bool
        Plots the transtion temperatures vs magnetic field.
    F_M_vs_T : bool
        Plots the magnetic free energy vs temperature.
    F_M_vs_M : bool
        Plots the magnetic free energy vs reduced magnetization.
    F_L_vs_T : bool
        Plots the lattice free energy vs temperature.
    FtotvsM : bool
        Plots the total free energy as a function of the reduced magnetization.
    Ftot_heatcool : bool
        Plots the free energies following the metastable minimums on heating
        and cooling.
    """
    print 'Free Energies'

    if FvsT or FvsB:
        free_ener_1 = free.F(TT, BB, J1, gJ, TC1, lamb1, Nm, theta_D1, N, F01)
        free_ener_2 = free.F(TT, BB, J2, gJ, TC2, lamb2, Nm, theta_D2, N, F02)
        free_stable = free.F_stable(TT, BB, J1, J2, gJ, TC1, TC2, lamb1, lamb2,
                                    Nm, theta_D1, theta_D2, N, F01, F02)

        if FvsT:
            plt_F_vs_T(Delta_T, Delta_B, free_ener_1, free_ener_2,
                       free_stable, Bstep=1, save=save)

        if FvsB:
            plt_F_vs_B(Delta_T, Delta_B, free_ener_1, free_ener_2,
                       free_stable, Tstep=10, save=save)

        if trans_temp:
            plt_transition_temp(TT, BB, save, J1, TC1, theta_D1, F01, lamb1,
                                J2, TC2, theta_D2, F02, lamb2, gJ, Nm, N,
                                free_ener_1, free_ener_2)

    # if the free energies are not already calculated
    if trans_temp and (FvsT != True and FvsB != True):
        plt_transition_temp(TT, BB, save, J1, TC1, theta_D1, F01, lamb1, J2,
                            TC2, theta_D2, F02, lamb2, gJ, Nm, N)

    if F_M_vs_T:
        plt_F_M_vs_T(TT, BB, J1, J2, TC1, TC2, lamb1,
                     lamb2, gJ, Nm, Bstep=1, save=save)

    if F_M_vs_M:
        # plots the magnetic free energy vs magnetiation
        plt_F_M_vs_M(sigma, TT, BB, J1, TC1, lamb1, J2,
                     TC2, lamb2, gJ, Nm, Tstep=10, Bstep=1)

    if F_L_vs_T:
        plt_F_L_vs_T(Delta_T, theta_D1, theta_D2)

    if FtotvsM:
        plt_Ftot_vs_M(sigma, Delta_T, Delta_B, J1, J2, TC1, TC2, lamb1, lamb2,
                      theta_D1, theta_D2, F01, F02, gJ, Nm, N, Tstep=20,
                      Bstep=1, save=save)

    if Ftot_heatcool:
        f_heat = free.F_tot_stable_Heating(TT, BB, J1, TC1, lamb1, theta_D1,
                                           F01, J2, TC2, lamb2, theta_D2, F02,
                                           gJ, Nm, N)
        f_cool = free.F_tot_stable_Cooling(TT, BB, J1, TC1, lamb1, theta_D1,
                                           F01, J2, TC2, lamb2, theta_D2, F02,
                                           gJ, Nm, N)

        plt_F_heatcool_vs_T(Delta_T, Delta_B, f_heat,
                            f_cool, Bstep=1, save=save)

        plt_F_heatcool_vs_B(Delta_T, Delta_B, f_heat,
                            f_cool, Tstep=10, save=save)

    return None


def plt_E(EMvsT=0, ELvsT=0, EtotvsT=0, E_M_vs_M=0):
    """Menu for the energy plots.

    Parameters
    ----------
    EMvsT : bool
        Plots the magnetic energy vs temperature.
    ELvsT : bool
        Plots the lattice energy vs magnetic field.
    EtotvsT : bool
        Plots the total energy vs temperature.
    E_M_vs_M : bool
        Plots the magnetic energy vs reduced magnetization.
    """
    print 'Energy'

    if EMvsT:
        em1 = Nm * ener.E_M(TT, BB, J1, TC1, lamb1)
        em2 = Nm * ener.E_M(TT, BB, J2, TC2, lamb2)

        plt_E_M_vs_T(Delta_T, Delta_B, em1, em2, Bstep=1)

    if ELvsT:
        el1 = N * ener.E_L(Delta_T, theta_D1)
        el2 = N * ener.E_L(Delta_T, theta_D2)

        plt_E_L_vs_T(Delta_T, el1, el2)

    if EtotvsT:
        etot1 = ener.E_tot(TT, BB, J1, TC1, lamb1, theta_D1, F01)
        etot2 = ener.E_tot(TT, BB, J2, TC2, lamb2, theta_D2, F02)

        plt_Etot_vs_T(Delta_T, Delta_B, etot1, etot2, Bstep=1)

    if E_M_vs_M:
        plt_E_M_vs_M(sig, Delta_T, Delta_B, Tstep=30, Bstep=1, save=save)

    return None


if __name__ == "__main__":
    from variables import *

    print 'Plotting...\n'

    plt_M(MvsT=1, MvsB=0, MvsTB=0, UvsT=0, M_hys_vs_T=0, save=0, TT=TT, BB=BB,
          J1=J1, J2=J2, TC1=TC1, TC2=TC2, lamb1=lamb1, lamb2=lamb2,
          Delta_T=Delta_T, Delta_B=Delta_B, theta_D1=theta_D1,
          theta_D2=theta_D2, gJ=gJ, F01=F01, F02=F02, Nm=Nm, N=N)
    plt_S(S_M_vs_T=0, S_L_VS_T=0, S_tot_vs_T=0, DeltaS_vs_T=0,
          max_DeltaS_vs_B=0, S_M_vs_M=0, save=0)
    plt_F(FvsT=0, FvsB=0, trans_temp=0, F_M_vs_T=0, F_M_vs_M=0,
          F_L_vs_T=0, FtotvsM=0, Ftot_heatcool=0, save=0)
    plt_E(EMvsT=0, ELvsT=0, EtotvsT=0, E_M_vs_M=0)

    plt.show()
