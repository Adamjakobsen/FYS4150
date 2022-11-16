import numpy as np 
import matplotlib.pyplot as plt

#20x20 lattice 
#plotting average energy per spin as a function of MC cycles for two temperatures: T = 1 and T = 2.4

#import quantities 
#need to import data for L=20, mc1000000, burn0 and t1
#need to import data for L=20, mc1000000, burn0 and t2.4
data_t1 = np.loadtxt("./20/qt_L20_A0_mc1000_burn0_t1.000.txt")
data_t24 = np.loadtxt("./20/qt_L20_A0_mc1000_burn0_t1.000.txt")

#import quantities for T = 1
E_N1 = data_t1[:,0] #average energy per spin 
M_N1 = data_t1[:,1] #average magnetisation per spin 
C_v1 = data_t1[:,2] #specific heat capacity
chi1 = data_t1[:,3] #susceptibility

#import quantities for T=2.4
E_N24 = data_t24[:,0] #average energy per spin 
M_N24 = data_t24[:,1] #average magnetisation per spin 
C_v24 = data_t24[:,2] #specific heat capacity
chi24 = data_t24[:,3] #susceptibility

#creating normalised histograms of generated epsilon samples for T = 1 and T = 2.4
#have burn in of 10% of total MC cycles
epsilon = np.loadtxt("./problem6/20/epsilon_L20_A0_mc1000_burn0_t1.000.txt")
epsilon = epsilon #normalising by number of spins
plt.hist(epsilon, bins = 100, density = True, grid = True)
plt.show()

