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
data_t1_align = read_csv("../data/20/all_qt_L20_A1_mc1000000_burn0_lt1.000_ut1.000.csv")
data_t24_align = read_csv("../data/20/all_qt_L20_A1_mc1000000_burn0_lt2.400_ut2.400.csv")

data_t1_random = read_csv("../data/20/all_qt_L20_A0_mc1000000_burn0_lt1.000_ut1.000.csv")
data_t24_random = read_csv("../data/20/all_qt_L20_A0_mc1000000_burn0_lt2.400_ut2.400.csv")

#import energy per spin values with mc = 10^6, T = 1 and T = 2.4 and burn in  = 1% of mc cycles with random configurations
#for histograms
epsilon_t1_burn1 = read_csv("../data/20/epsilons_L20_A0_mc1000000_burn1_lt1.000_ut1.000.csv")
epsilon_t24_burn1 = read_csv("../data/20/epsilons_L20_A0_mc1000000_burn1_lt2.400_ut2.400.csv")


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
#plt.rcParams['lines.markersize'] = 50
plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'serif'
plt.tight_layout()
plt.subplots_adjust(left=0.1, bottom=0.17, right=0.97, top=0.97, wspace=0.2, hspace=0.2)

#T = 1, aligned and random
ax[0].plot(e_1_align, label = "$L = 20, T = 1 \; [J/k_{B}]$ aligned config", color = 'purple', lw = 2.2)
ax[0].plot(e_1_random, label = "$L = 20, T = 1 \; [J/k_{B}]$ random config", color = 'red', lw = 2.2)
ax[0].set_xscale("log")
ax[0].legend() 
ax[0].set_xlabel("MC cycles")
ax[0].set_ylabel("$\langle \epsilon \\rangle \;[J]$")
ax[0].tick_params(axis='both', which='major', labelsize=14)

#T = 2.4, aligned and random 
ax[1].plot(e_24_align, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ aligned config",  color = 'purple', lw = 2.2)
ax[1].plot(e_24_random, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ random config",  color = 'red', lw = 2.2)
ax[1].set_xscale("log")
ax[1].legend()
ax[1].set_xlabel("MC cycles")
ax[1].set_ylabel("$\langle \epsilon \\rangle \;[J]$")
ax[1].tick_params(axis='both', which='major', labelsize=14)

ax[1].yaxis.label.set_size(15)
ax[0].yaxis.label.set_size(15)
ax[1].xaxis.label.set_size(15)
ax[0].xaxis.label.set_size(15)



fig.savefig("../figs/energy_per_spin_L20.pdf") 

# %%
#plot average magnetisations per spin as a function of MC cycles for T = 1 and T = 2.4, aligned and random
fig, ax = plt.subplots(1, 2, figsize=(15,5))
#plt.rcParams['lines.markersize'] = 50
plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'serif'
plt.tight_layout()
plt.subplots_adjust(left=0.1, bottom=0.17, right=0.97, top=0.97, wspace=0.2, hspace=0.2)

#T = 1, aligned and random
ax[0].plot(m_1_align, label = "$L = 20, T = 1 \; [J/k_{B}]$ aligned config", color = 'purple', lw = 2.2)
ax[0].plot(m_1_random, label = "$L = 20, T = 1 \; [J/k_{B}]$ random config", color = 'red', lw = 2.2)

ax[0].set_xscale("log")
ax[0].legend() 
ax[0].set_xlabel("MC cycles")
ax[0].set_ylabel("$\langle |m| \\rangle \;$")
ax[0].tick_params(axis='both', which='major', labelsize=14)


#T = 2.4, aligned and random
ax[1].plot(m_24_align, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ aligned config",  color = 'purple', lw = 2.2)
ax[1].plot(m_24_random, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ random config",  color = 'red', lw = 2.2)
ax[1].set_xscale("log")
ax[1].legend()
ax[1].set_xlabel("MC cycles")
ax[1].set_ylabel("$\langle |m| \\rangle \;$")
ax[1].tick_params(axis='both', which='major', labelsize=14)

ax[1].yaxis.label.set_size(15)
ax[0].yaxis.label.set_size(15)
ax[1].xaxis.label.set_size(15)
ax[0].xaxis.label.set_size(15)


fig.savefig("../figs/magnetisation_per_spin_L20.pdf")

#plot histograms of energy per spin for T = 1 and T = 2.4 to visualise probability distributions
#define fig and ax for suplots 
fig, ax = plt.subplots(1, 2, figsize=(15,5))
plt.rcParams['lines.markersize'] = 10
plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'serif'

plt.tight_layout()
plt.subplots_adjust(left=0.1, bottom=0.17, right=0.97, top=0.97, wspace=0.2, hspace=0.2)

a = ax[0].hist(epsilon_t1_burn1, bins = 5, density = True, color = 'purple', label = "$L = 20, T = 1 \; [J/k_{B}]$, burn-in = 1%")
# ax1.hist(epsilon_t24_burn0, bins = 200, density = True, color = 'red', label = "$L = 20, T = 1 \; [J/k_{B}]$, burn-in = 0")
b = ax[1].hist(epsilon_t24_burn1, bins = 300, density = True, color = 'purple', label = "$L = 20, T = 2.4 \; [J/k_{B}]$, burn-in = 1%")
ax[0].set_xlabel("$\epsilon \; [J]$")
ax[0].set_ylabel("$p(\epsilon)$")
ax[1].set_xlabel("$\epsilon \; [J]$")
ax[1].set_ylabel("$p(\epsilon)$")
ax[0].tick_params(axis='both', which='major', labelsize=14)
ax[1].tick_params(axis='both', which='major', labelsize=14)
ax[1].yaxis.label.set_size(15)
ax[0].yaxis.label.set_size(15)
ax[1].xaxis.label.set_size(15)
ax[0].xaxis.label.set_size(15)


ax[0].legend()
ax[1].legend()
fig.savefig("../figs/histograms_L20.pdf")