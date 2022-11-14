import numpy as np
import matplotlib.pyplot as plt
import random

#define constants 
T = 2.4 #J/k_B
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
    magnet_diff = exp_M2 - (exp_M)**2
    chi = prefactor*magnet_diff
    return chi

C_v = specific_heat_cap()
chi = susceptibility()

data = np.array([exp_E_N, exp_E2_N, exp_M_N, exp_M2_N, C_v, chi])

#hacky way of saving 1D numpy array as a row and not a column
file = np.savetxt('analytical_values_L_2.txt', [data], delimiter = ' ')


data_read = np.loadtxt('analytical_values_L_2.txt')
#print out like this: 
# print(data_read[1])


