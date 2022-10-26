import numpy as np
import matplotlib.pyplot as plt


 

n_particles=100

data1=np.loadtxt("../Data/RK4_perturbation_100_50000_off0.100000.txt")
data2=np.loadtxt("../Data/RK4_perturbation_100_50000_off0.400000.txt")
data3=np.loadtxt("../Data/RK4_perturbation_100_50000_off0.700000.txt")
n1=data1[:,1]
n2=data2[:,1]
n3=data3[:,1]

omega=data1[:,0]

plt.plot(omega,(n_particles-n1)/n_particles,linewidth=0.5, label="f=0.1")
plt.plot(omega,(n_particles-n2)/n_particles,linewidth=0.5, label="f=0.4")
plt.plot(omega,(n_particles-n3)/n_particles,linewidth=0.5, linestyle='dashed', label="f=0.7")
plt.xlabel(r'$\omega$')
plt.ylabel('Particle fraction inside trap')
plt.legend()
plt.show()