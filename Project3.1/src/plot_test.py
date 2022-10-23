import numpy as np
import matplotlib.pyplot as plt



#Taking command line input argument: method
import sys
method = sys.argv[1]
n_particles = 2

if method=='Euler':
    #Euler's method
    data = np.loadtxt('../positions_Euler.txt')

    n_timesteps = int(data[:,0].size/n_particles)


    x_E = np.zeros((n_particles, n_timesteps))
    y_E = np.zeros((n_particles, n_timesteps))
    z_E = np.zeros((n_particles, n_timesteps))

    for i in range(n_particles):
        x_E[i,:] = data[i*n_timesteps:(i+1)*n_timesteps,0]
        y_E[i,:] = data[i*n_timesteps:(i+1)*n_timesteps,1]
        z_E[i,:] = data[i*n_timesteps:(i+1)*n_timesteps,2]

        plt.plot(x_E[i,:], y_E[i,:], label=f'Particle {i+1}')

    plt.axis('equal')
    plt.legend()
    plt.show()

if method=='RK4':
    #RK4 method
    data = np.loadtxt('../positions_rk4.txt')

    n_timesteps = int(data[:,0].size/n_particles)


    x_RK4 = np.zeros((n_particles, n_timesteps))
    y_RK4 = np.zeros((n_particles, n_timesteps))
    z_RK4 = np.zeros((n_particles, n_timesteps))

    for i in range(n_particles):
        x_RK4[i,:] = data[i*n_timesteps:(i+1)*n_timesteps,0]
        y_RK4[i,:] = data[i*n_timesteps:(i+1)*n_timesteps,1]
        z_RK4[i,:] = data[i*n_timesteps:(i+1)*n_timesteps,2]

        plt.plot(x_RK4[i,:], y_RK4[i,:], label=f'Particle {i+1}')

    plt.axis('equal')
    plt.legend()
    plt.show()    




