# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 01:23:52 2016

@author: Utilizador
"""


def plt_susc():
    print 'Susceptibility  as a function of Temperature'
    
    susc_1 = mag.Susceptibility(TT, BB, J1, TC1, lamb1)
    susc_2 = mag.Susceptibility(TT, BB, J2, TC2, lamb2)
    
    plt.figure()
    for i in range(1, len(Delta_B), 10):
        plt.plot(Delta_T, 1./susc_1[i], label='1, B='+str(Delta_B[i])+'T')
        plt.plot(Delta_T, 1./susc_2[i], label='2, B='+str(Delta_B[i])+'T')
    
    plt.title('Inverse Susceptibility')
    plt.legend(loc=0)
    #plt.ylim(0,1.05)
    plt.ylabel('$1/ \chi(T, B)$')
    plt.xlabel('T (K)')

    return None
    
#==============================================================================
#     
#==============================================================================

# Total Free Energy Equilibrium
def F_eq(T, B):
    etot1 = ener.E_tot(T, B, J1, TC1, lamb1, theta_D1, F01)
    etot2 = ener.E_tot(T, B, J2, TC2, lamb2, theta_D2, F02)
    
    g = 1./(1. + np.exp(-1.*(etot2-etot1)/(k_B*T)))
        
    
#==============================================================================
#     f1 = F(T, B, J1, TC1, theta_D1, F01, lamb1)
#     f2 = F(T, B, J2, TC2, theta_D2, F02, lamb2)
#     
#     g = 1./(1. + np.exp(-1.*(f2-f1)/(k_B*T)))
#==============================================================================

    
    plt.figure()
    plt.plot(Delta_T, g[0])
    plt.title('$g \Rightarrow F=gF_1 + (1-g)F_2$, B='+str(Delta_B[0])+'T')
    plt.xlabel('T(K)')
    plt.ylabel('g')
    plt.show()
        
    F1 = F(T, B, J1, TC1, theta_D1, F01, lamb1)
    F2 = F(T, B, J2, TC2, theta_D2, F02, lamb2)
    
    return g*F1 + (1.-g)*F2
    

#==============================================================================
# 
#==============================================================================


    if FvsM == True:
        print '\t Magnetic Free Energy as a function of Magnetization'
        
        for j in range(len(Delta_B)):
            
            plt.figure()
            for i in range(0, len(Delta_T), 20):
                f1 = free.F_M_vs_M(sig, TT[0,i], BB[0,0], J1, TC1, lamb1)
                plt.plot(sig, f1,':', label='1, T='+str(Delta_T[i])+'K')
            plt.gca().set_color_cycle(None)            
            for i in range(0, len(Delta_T), 20):
                f2 = free.F_M_vs_M(sig, TT[0,i], BB[0,0], J2, TC2, lamb2)
                plt.plot(sig, f2,'--', label='2, T='+str(Delta_T[i])+'K')
 
            plt.legend(loc=0, fontsize='small')
            plt.title('Magnetic Free Energy as a function of Magnetization $F^M$, B='+str(Delta_B[j])+'T')
            plt.xlabel('$\sigma$')
            plt.ylabel('F(T) (eV)')
            plt.xlim(-1,1)
            
            
#==============================================================================
# 
#==============================================================================


















