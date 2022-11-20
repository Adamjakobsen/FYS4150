import numpy as np
import matplotlib.pyplot as plt
import random
from pandas import read_csv

#define constants 
T = 1 #J/k_B
J = 1
k_B = 1
N = 4
beta = 1/(k_B*T)

#define equipartition function 
Z = 2*np.exp(8*beta*J) + 2*np.exp(-8*beta*J) + 12

#make functions for computing analytical thermodynamic quantities 

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
e = av_energy_per_spin()
exp_E2 = av_energy_sqrt()
e2 = av_energy_per_spin_sqrt()
exp_M = av_magnet() 
m = av_magnet_per_spin()
exp_M2 = av_magnet_sqrt()
m2 = av_magnet_sqrt_per_spin()

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

#save analytical values in a separate file 
data = np.array([e, m, C_v_analyt, chi_analyt])

#hacky way of saving 1D numpy array as a row and not a column
file = np.savetxt('analytical_values_L_2.txt', [data], delimiter = ' ')

"""Computing relative errors to compare analytical with numerical """
#no burn time implemented 
#import data using pandas since its .csv files 
data_qts_mc3 = read_csv('./2/qt_L2_A0_mc1000_burn0_t1.000.csv')
data_qts_mc6 = read_csv('./2/qt_L2_A0_mc1000000_burn0_t1.000.csv')
data_qts_mc8 = read_csv('./2/qt_L2_A0_mc100000000_burn0_t1.000.csv')

#turn csv into numpy array 
qts_mc3 = data_qts_mc3.to_numpy() #output: avg_e << " " << avg_mabs << " " << Cv << " " << chi
qts_mc6 = data_qts_mc6.to_numpy()
qts_mc8 = data_qts_mc8.to_numpy()

#extract last values from the simulated data file to get the most accurate MC cycle output 
def relative_error(data):
    exp_e_mc = data[-1][0]
    exp_m_mc = data[-1][1]
    C_v_mc = data[-1][2]
    chi_mc = data[-1][3]
    
    delta_exp_E_N3 = abs((exp_e_mc - e)/e)
    delta_exp_M_N3 = abs((exp_m_mc - m)/m)
    delta_C_v3 = abs((C_v_mc - C_v_analyt)/C_v_analyt)
    delta_chi3 = abs((chi_mc - chi_analyt)/chi_analyt)
    
    return(delta_exp_E_N3, delta_exp_M_N3, delta_C_v3, delta_chi3)

mc_3 = relative_error(qts_mc3)
mc_6 = relative_error(qts_mc6)
mc_8 = relative_error(qts_mc8)

print("Relative errors for MC = 10^3")
print("Relative error of exp. energy per spin = ", mc_3[0])
print("Relative error of exp. magnetisation per spin = " , mc_3[1])
print("Relative error of exp. C_v = ", mc_3[2])
print("Relative error of exp. chi = ", mc_3[3])

print("Relative errors for MC = 10^6")
print("Relative error of exp. energy per spin = ", mc_6[0])
print("Relative error of exp. magnetisation per spin = " , mc_6[1])
print("Relative error of exp. C_v = ", mc_6[2])
print("Relative error of exp. chi = ", mc_6[3])

print("Relative errors for MC = 10^8")
print("Relative error of exp. energy per spin = ", mc_8[0])
print("Relative error of exp. magnetisation per spin = " , mc_8[1])
print("Relative error of exp. C_v = ", mc_8[2])
print("Relative error of exp. chi = ", mc_8[3])


