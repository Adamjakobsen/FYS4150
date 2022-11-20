# %%
import numpy as np 
import matplotlib.pyplot as plt
#import pandas
from pandas import read_csv
import statistics as st
import os 


#20x20 lattice 
#plotting average energy per spin as a function of MC cycles for two temperatures: T = 1 and T = 2.4

#import quantities for plots of quantities as a function of MC cycles 
data_t1_align = read_csv("./20/qt_L20_A1_mc1000000_burn0_t1.000.csv")
data_t24_align = read_csv("./20/qt_L20_A1_mc1000000_burn0_t2.400.csv")

data_t1_random = read_csv("./20/qt_L20_A0_mc1000000_burn0_t1.000.csv")
data_t24_random = read_csv("./20/qt_L20_A0_mc1000000_burn0_t2.400.csv")

#import energy per spin values with mc = 10^6, T = 1 and T = 2.4 and burn in  = 1% of mc cycles with random configurations
#for histograms
epsilon_t1_burn1 = read_csv("./20/epsilons_L20_A0_mc1000000_burn1_t1.000.csv")
epsilon_t24_burn1 = read_csv("./20/epsilons_L20_A0_mc1000000_burn1_t2.400.csv")


# %%
#turn pandas dataframe into numpy array
data_t1_align = data_t1_align.to_numpy()
data_t24_align = data_t24_align.to_numpy()

data_t1_random = data_t1_random.to_numpy()
data_t24_random = data_t24_random.to_numpy()

#epsilon is for histogramas
epsilon_t1_burn1 = epsilon_t1_burn1.to_numpy()
epsilon_t24_burn1 = epsilon_t24_burn1.to_numpy()

#Plot energies per spin as a function of MC cycles 
#import quantities for T = 1, aligned
e_1_align = data_t1_align[:,0] #average energy per spin 
m_1_align = data_t1_align[:,1] #average magnetisation per spin 

#import quantities for T=2.4, aligned
e_24_align = data_t24_align[:,0] #average energy per spin 
m_24_align = data_t24_align[:,1] #average magnetisation per spin 

#import quantities for T = 1, random
e_1_random = data_t1_random[:,0] #average energy per spin
m_1_random = data_t1_random[:,1] #average magnetisation per spin

#import quantities for T = 2.4, random
e_24_random = data_t24_random[:,0] #average energy per spin
m_24_random = data_t24_random[:,1] #average magnetisation per spin


#plot average energies per spin as a function of MC cycles for T = 1 and T = 2.4, aligned and random
fig, ax = plt.subplots(1, 2, figsize=(15,5))
#T = 1, aligned and random
ax[0].plot(e_1_align, label = "$L = 20, T = 1 \; [J/k_{B}]$ aligned config", color = 'purple')
ax[0].plot(e_1_random, label = "$L = 20, T = 1 \; [J/k_{B}]$ random config", color = 'red')
ax[0].set_xscale("log")
ax[0].grid()
ax[0].legend() 
ax[0].set_xlabel("MC cycles")
ax[0].set_ylabel("$\langle \epsilon \\rangle \;[J]$")
#T = 2.4, aligned and random 
ax[1].plot(e_24_align, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ aligned config",  color = 'purple')
ax[1].plot(e_24_random, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ random config",  color = 'red')
ax[1].set_xscale("log")
ax[1].legend()
ax[1].grid()
ax[1].set_xlabel("MC cycles")
ax[1].set_ylabel("$\langle \epsilon \\rangle \;[J]$")
fig.savefig("./figs/energy_per_spin_L20.pdf") 

# %%
#plot average magnetisations per spin as a function of MC cycles for T = 1 and T = 2.4, aligned and random
fig, ax = plt.subplots(1, 2, figsize=(15,5))
#T = 1, aligned and random
ax[0].plot(m_1_align, label = "$L = 20, T = 1 \; [J/k_{B}]$ aligned config", color = 'purple')
ax[0].plot(m_1_random, label = "$L = 20, T = 1 \; [J/k_{B}]$ random config", color = 'red')
ax[0].set_xscale("log")
ax[0].grid()
ax[0].legend() 
ax[0].set_xlabel("MC cycles")
ax[0].set_ylabel("$\langle |m| \\rangle \;[1]$")

#T = 2.4, aligned and random
ax[1].plot(m_24_align, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ aligned config",  color = 'purple')
ax[1].plot(m_24_random, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ random config",  color = 'red')
ax[1].set_xscale("log")
ax[1].legend()
ax[1].grid()
ax[1].set_xlabel("MC cycles")
ax[1].set_ylabel("$\langle |m| \\rangle \;[1]$")
fig.savefig("./figs/magnetisation_per_spin_L20.pdf")

#plot histograms of energy per spin for T = 1 and T = 2.4 to visualise probability distributions
#define fig and ax for suplots 
fig, ax = plt.subplots(1, 2, figsize=(15,5))
a = ax[0].hist(epsilon_t1_burn1, bins = 300, density = True, color = 'red', label = "$L = 20, T = 1 \; [J/k_{B}]$, burn-in = 1%")
# ax1.hist(epsilon_t24_burn0, bins = 200, density = True, color = 'red', label = "$L = 20, T = 1 \; [J/k_{B}]$, burn-in = 0")
b = ax[1].hist(epsilon_t24_burn1, bins = 300, density = True, color = 'red', label = "$L = 20, T = 2.4 \; [J/k_{B}]$, burn-in = 1%")
ax[0].legend()
ax[1].legend()
fig.savefig("./figs/histograms_L20.pdf")







