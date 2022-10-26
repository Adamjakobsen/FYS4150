import numpy as np
import matplotlib.pyplot as plt


 

n_particles=100

#Plotting the number of particles inside the trap as a function of frequencies

#Interactions of
data1=np.loadtxt("../Data/RK4_perturbation_100_50000_off0.100000.txt")
data2=np.loadtxt("../Data/RK4_perturbation_100_50000_off0.400000.txt")
data3=np.loadtxt("../Data/RK4_perturbation_100_50000_off0.700000.txt")
n1=data1[:,1]
n2=data2[:,1]
n3=data3[:,1]

omega=data1[:,0]

plt.plot(omega,(n_particles-n1)/n_particles,linewidth=0.6, label="f=0.1")
plt.plot(omega,(n_particles-n2)/n_particles,linewidth=0.6, label="f=0.4")
plt.plot(omega,(n_particles-n3)/n_particles,linewidth=0.6, linestyle='dashed', label="f=0.7")
plt.xlabel(r'$\omega_{V}$ [Mhz]')
plt.ylabel('Particle fraction inside trap')
plt.legend()
plt.savefig("../Fig/timedep_100p_interaction_off.pdf")
plt.show()

#Interactions of
data1=np.loadtxt("../Data/RK4_perturbation_100_50000_on0.100000.txt")
data2=np.loadtxt("../Data/RK4_perturbation_100_50000_on0.400000.txt")
data3=np.loadtxt("../Data/RK4_perturbation_100_50000_on0.700000.txt")
n1=data1[:,1]
n2=data2[:,1]
n3=data3[:,1]

omega=data1[:,0]

plt.plot(omega,(n_particles-n1)/n_particles,linewidth=0.6, label="f=0.1")
plt.plot(omega,(n_particles-n2)/n_particles,linewidth=0.6, label="f=0.4")
plt.plot(omega,(n_particles-n3)/n_particles,linewidth=0.6, linestyle='dashed', label="f=0.7")
plt.xlabel(r'$\omega_{V}$ [Mhz]')
plt.ylabel('Particle fraction inside trap')
plt.legend()
plt.savefig("../Fig/timedep_100p_interaction_on.pdf")
plt.show()

# Fine grained plot
data_fine=np.loadtxt("../Data/Fine_perturbation_100_50000_on0.400000.txt")
data_fine_off=np.loadtxt("../Data/Fine_perturbation_100_50000_off0.400000.txt")
n_fine=data_fine[:,1]
n_fine_off=data_fine_off[:,1]
omega_fine=data_fine[:,0]
plt.plot(omega_fine,(n_particles-n_fine)/n_particles,linewidth=0.6, label="With interaction ")
plt.plot(omega_fine,(n_particles-n_fine_off)/n_particles,linewidth=0.6, label="Without interaction")
plt.xlabel(r'$\omega_{V} [Mhz]$')
plt.ylabel('Particle fraction inside trap')
plt.legend()
plt.savefig("../Fig/Fine_timedep_100p_interaction_0.4.pdf")
plt.show()