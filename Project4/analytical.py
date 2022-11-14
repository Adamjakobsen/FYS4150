import numpy as np
import matplotlib.pyplot as plt
import random

#define constants 
T = 1 #J/k_B
J = 1
k_B = 1
N = 4
beta = 1/(k_B*T)
L = 2

#define equipartition function 
Z = 2*np.exp(8*beta*J) + 2*np.exp(-8*beta*J) + 12

#make functions for desired quantitites 

def av_energy():
    prefactor = -J/Z
    av_energy = prefactor * (16*np.exp(8*beta*J)-16*np.exp(-8*beta*J))
    return av_energy

def av_energy_per_spin():
    prefactor = -(4*J)/Z
    av_energy_per_spin = prefactor * (np.exp(8*J*beta) - np.exp(-8*beta*J))
    return av_energy_per_spin

def av_energy_sqrt():
    prefactor = (128*J**2)/Z
    av_energy_sqrt = prefactor*(np.exp(8*beta*J) + np.exp(-8*beta*J))
    return av_energy_sqrt

def av_energy_per_spin_sqrt():
    prefactor = (8*J**2)/Z
    av_energy_per_spin_sqrt = prefactor * (np.exp(8*J*beta) + np.exp(-8*beta*J))
    return av_energy_per_spin_sqrt

def av_magnet():
    prefactor = 1/Z
    av_magnet = prefactor*(8*np.exp(8*beta*J) + 16)
    return av_magnet 

def av_magnet_per_spin():
    prefactor = 1/Z
    av_magnet_per_spin = prefactor*(2*np.exp(8*beta*J) + 4)
    return av_magnet_per_spin
    
def av_magnet_sqrt():
    prefactor = 1/Z
    av_magnet_sqrt = prefactor*(32*np.exp(8*beta*J) + 32)
    return av_magnet_sqrt

def av_magnet_sqrt_per_spin():
    prefactor = 1/Z
    av_magnet_sqrt_per_spin = prefactor * (2*np.exp(8*beta*J) + 2)
    return av_magnet_sqrt_per_spin



exp_E = av_energy() 
exp_E_N = av_energy_per_spin()
exp_E2 = av_energy_sqrt()
exp_E2_N = av_energy_per_spin_sqrt()
exp_M = av_magnet() 
exp_M_N = av_magnet_per_spin()
exp_M2 = av_magnet_sqrt()
exp_M2_N = av_magnet_sqrt_per_spin()

#normalised specific heat capacity 
def specific_heat_cap():
    prefactor = 1/N * (1/(k_B*T**2))
    energy_diff = exp_E2 - (exp_E)**2
    C_v = prefactor*energy_diff
    return C_v

#normalised susceptibility 
def susceptibility():
    prefactor = 1/N * (1/(k_B*T))
    magnet_diff = exp_M2 - (abs(exp_M))**2
    chi = prefactor*magnet_diff
    return chi

C_v_analyt = specific_heat_cap()
chi_analyt = susceptibility()

data = np.array([exp_E_N, exp_M_N, C_v_analyt, chi_analyt])

#hacky way of saving 1D numpy array as a row and not a column
file = np.savetxt('analytical_values_L_2.txt', [data], delimiter = ' ')

"""Computing relative errors to compare analytical with numerical """

data_sim3 = np.loadtxt('quantities_L_2_mc_1000.txt') #output: avg_e << " " << avg_mabs << " " << Cv << " " << chi
# data_sim6 = np.loadtxt('quantities_L_2_mc_1000000.txt')
# data_sim8 = np.loadtxt('quantities_L_2_mc_100000000.txt')

#extract last values from the simulated data file to get the most accurate MC cycle output 
def relative_error(data):
    exp_E_N_mc = data[-1][0]
    exp_M_N_mc = data[-1][1]
    C_v_mc = data[-1][2]
    chi_mc = data[-1][3]
    
    delta_exp_E_N3 = abs((exp_E_N_mc - exp_E_N)/exp_E_N)
    delta_exp_M_N3 = abs((exp_M_N_mc - exp_M_N)/exp_M_N)
    delta_C_v3 = abs((C_v_mc - C_v_analyt)/C_v_analyt)
    delta_chi3 = abs((chi_mc - chi_analyt)/chi_analyt)
    
    return(delta_exp_E_N3, delta_exp_M_N3, delta_C_v3, delta_chi3)
    
mc_3 = relative_error(data_sim3)
# mc_6 = relative_error(data_sim6)
# mc_8 = relative_error(data_sim8)

print("Relative errors for MC = 10^3")
print("Relative error of exp. energy per spin = ", mc_3[0])
print("Relative error of exp. magnetisation per spin = " , mc_3[1])
print("Relative error of exp. C_v = ", mc_3[2])
print("Relative error of exp. chi = ", mc_3[3])


# print("Relative errors for MC = 10^6")
# print("Relative error of exp. energy per spin = ", mc_3[0])
# print("Relative error of exp. magnetisation per spin = " , mc_3[1])
# print("Relative error of exp. C_v = ", mc_3[2])
# print("Relative error of exp. chi = ", mc_3[3])

# print("Relative errors for MC = 10^8")
# print("Relative error of exp. energy per spin = ", mc_3[0])
# print("Relative error of exp. magnetisation per spin = " , mc_3[1])
# print("Relative error of exp. C_v = ", mc_3[2])
# print("Relative error of exp. chi = ", mc_3[3])

# data_sim5 = np.loadtxt('quantities_L_2_mc_1000.txt')
# mc_5 = relative_error(data_sim5)
# print("Relative errors for MC = 10^5")
# print("Relative error of exp. energy per spin = ", mc_3[0])
# print("Relative error of exp. magnetisation per spin = " , mc_3[1])
# print("Relative error of exp. C_v = ", mc_3[2])
# print("Relative error of exp. chi = ", mc_3[3])

