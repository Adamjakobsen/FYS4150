# %%
import numpy as np 
import matplotlib.pyplot as plt
#import pandas
from pandas import read_csv
import statistics as st

# %% [markdown]
# # Import data 
# 

# %%
#20x20 lattice 
#plotting average energy per spin as a function of MC cycles for two temperatures: T = 1 and T = 2.4

#import quantities 

data_t1_align = read_csv("./20/qt_L20_A1_mc1000000_burn0_t1.000.csv")
data_t24_align = read_csv("./20/qt_L20_A1_mc1000000_burn0_t2.400.csv")

data_t1_random = read_csv("./20/qt_L20_A0_mc1000000_burn0_t1.000.csv")
data_t24_random = read_csv("./20/qt_L20_A0_mc1000000_burn0_t2.400.csv")

#import energy per spin values with mc = 10^6, T = 1 and T = 2.4 and burn in  = 1% of mc cycles with random configurations
epsilon_t1_burn1 = read_csv("./20/epsilons_L20_A0_mc1000000_burn1_t1.000.csv")
epsilon_t24_burn1 = read_csv("./20/epsilons_L20_A0_mc1000000_burn1_t2.400.csv")






# %%
#turn pandas dataframe into numpy array
data_t1_align = data_t1_align.to_numpy()
data_t24_align = data_t24_align.to_numpy()

data_t1_random = data_t1_random.to_numpy()
data_t24_random = data_t24_random.to_numpy()

epsilon_t1_burn1 = epsilon_t1_burn1.to_numpy()
epsilon_t24_burn1 = epsilon_t24_burn1.to_numpy()



# %% [markdown]
# # Plot energies per spin as a function of MC cycles 

# %%
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
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,5))

N = 100
# markersize = 5
# e_1_align = e_1_align[::N]
# e_24_align = e_24_align[::N]

# e_1_random= e_1_random[::N]
# e_24_random = e_24_random[::N]

ax1.plot(e_1_align, label = "$L = 20, T = 1 \; [J/k_{B}]$ aligned config", color = 'purple')
ax1.plot(e_1_random, label = "$L = 20, T = 1 \; [J/k_{B}]$ random config", color = 'red')
ax1.set_xscale("log")
ax1.grid()
ax1.legend() 
#set limit for x axis
# ax1.set_ylim(-2, -1.8)
ax1.set_xlabel("MC cycles")
ax1.set_ylabel("$\langle \epsilon \\rangle \;[J]$")

#plotting energies per spin for T = 2.4, aligned and random 
# extract every 50th value for plotting
ax2.plot(e_24_align, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ aligned config",  color = 'purple')
ax2.plot(e_24_random, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ random config",  color = 'red')
ax2.set_xscale("log")
ax2.legend()
ax2.grid()
# ax2.set_ylim(-1.31, -1.05)
ax2.set_xlabel("MC cycles")
ax2.set_ylabel("$\langle \epsilon \\rangle \;[J]$")
fig.savefig("energy_per_spin_L20.pdf")

# %%
#plot average magnetisations per spin as a function of MC cycles for T = 1 and T = 2.4, aligned and random
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,5))


ax1.plot(m_1_align, label = "$L = 20, T = 1 \; [J/k_{B}]$ aligned config", color = 'purple')
ax1.plot(m_1_random, label = "$L = 20, T = 1 \; [J/k_{B}]$ random config", color = 'red')
ax1.set_xscale("log")
ax1.grid()
ax1.legend() 
ax1.set_xlabel("MC cycles")
ax1.set_ylabel("$\langle |m| \\rangle \;[1]$")

#plotting magnetisations per spin for T = 2.4, aligned and random 
ax2.plot(m_24_align, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ aligned config",  color = 'purple')
ax2.plot(m_24_random, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ random config",  color = 'red')
ax2.set_xscale("log")
ax2.legend()
ax2.grid()
ax2.set_xlabel("MC cycles")
ax2.set_ylabel("$\langle |m| \\rangle \;[1]$")
fig.savefig("magnetisation_per_spin_L20.pdf")


# %% [markdown]
# ## Plot energies per spin for a zoomed in area (not necessary)

# %%
# #plot average energies per spin as a function of MC cycles for T = 1 and T = 2.4, aligned and random
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,5))

# N = 100
# markersize = 5
# e_1_align = e_1_align[::N]
# e_24_align = e_24_align[::N]

# e_1_random= e_1_random[::N]
# e_24_random = e_24_random[::N]

# ax1.plot(e_1_align, label = "$L = 20, T = 1 \; [J/k_{B}]$ aligned config", color = 'purple')
# ax1.plot(e_1_random, label = "$L = 20, T = 1 \; [J/k_{B}]$ random config", color = 'red')
# ax1.set_xscale("log")
# ax1.grid()
# ax1.legend() 
# #set limit for x axis
# ax1.set_ylim(-2, -1.8)
# ax1.set_xlabel("MC cycles")
# ax1.set_ylabel("$\langle \epsilon \\rangle \;[J]$")

# #plotting energies per spin for T = 2.4, aligned and random 
# # extract every 50th value for plotting
# ax2.plot(e_24_align, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ aligned config",  color = 'purple')
# ax2.plot(e_24_random, label = "$L = 20, T = 2.4 \; [J/k_{B}]$ random config",  color = 'red')
# ax2.set_xscale("log")
# ax2.legend()
# ax2.grid()
# ax2.set_ylim(-1.31, -1.05)
# ax2.set_xlabel("MC cycles")
# ax2.set_ylabel("$\langle \epsilon \\rangle \;[J]$")

# %% [markdown]
# # Plot histograms to approximate the probability function for energy per spin for T = 1 and T = 2.4

# %%

#define fig and ax for suplots 
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,5))
a = ax1.hist(epsilon_t1_burn1, bins = 300, density = True, color = 'red', label = "$L = 20, T = 1 \; [J/k_{B}]$, burn-in = 1%")
# ax1.hist(epsilon_t24_burn0, bins = 200, density = True, color = 'red', label = "$L = 20, T = 1 \; [J/k_{B}]$, burn-in = 0")
b = ax2.hist(epsilon_t24_burn1, bins = 400, density = True, color = 'red', label = "$L = 20, T = 2.4 \; [J/k_{B}]$, burn-in = 1%")
ax1.legend()
ax2.legend()





