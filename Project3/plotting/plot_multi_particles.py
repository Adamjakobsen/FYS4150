import numpy as np
import matplotlib.pyplot as plt
from Utils import *


#RK4
n_steps=50000
n_particles=10
filename=f"RK4_{n_particles}_{n_steps}off.txt"

utils=Utils(n_steps,50,10)
t,x,y,z,vx,vy,vz=utils.get_data(filename)

for i in range(n_particles):
    plt.plot(x[i,:],y[i,:],label=f"Particle {i+1}")

plt.show()



