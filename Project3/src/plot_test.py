import numpy as np
import matplotlib.pyplot as plt


#Euler's method
data = np.loadtxt('../positions_Euler.txt')

n_particles =2
n_timesteps = int(data[:,0].size/n_particles)

x_E = np.zeros((n_particles, n_timesteps))
y_E = np.zeros((n_particles, n_timesteps))
z_E = np.zeros((n_particles, n_timesteps))

for i in range(n_particles):
    x_E[i,:] = data[i*n_timesteps:(i+1)*n_timesteps,0]
    y_E[i,:] = data[i*n_timesteps:(i+1)*n_timesteps,1]
    z_E[i,:] = data[i*n_timesteps:(i+1)*n_timesteps,2]

    plt.plot(x_E[i,:], y_E[i,:], label=f'Particle {i+1}')

plt.legend()
plt.show()

#RK4
data = np.loadtxt('../positions_rk4.txt')
print(np.shape(data))
n_particles =2
n_timesteps = int(data[:,0].size/n_particles)

x_R = np.zeros((n_particles, n_timesteps))
y_R = np.zeros((n_particles, n_timesteps))
z_R = np.zeros((n_particles, n_timesteps))

for i in range(n_particles):
    x_R[i,:] = data[i*n_timesteps:(i+1)*n_timesteps,0]
    y_R[i,:] = data[i*n_timesteps:(i+1)*n_timesteps,1]
    z_R[i,:] = data[i*n_timesteps:(i+1)*n_timesteps,2]

    plt.plot(x_R[i,:], y_R[i,:], label=f'Particle {i+1}')

plt.legend()
plt.show()


plt.plot(x_R[1,:], y_R[1,:], label=f'Particle {2}')
plt.legend()
plt.show()