import numpy as np
import matplotlib.pyplot as plt
from Utils import *


#RK4
n_steps=32000
n_particles=2
filename=f"RK4_{n_particles}_{n_steps}on.txt"
#filename="RK4_2_50000.txt"
utils=Utils(n_steps,50,2)
t,x,y,z,vx,vy,vz=utils.get_data(filename)

#Particle one positions
x_particle1=x[0,:]
y_particle1=y[0,:]
z_particle1=z[0,:]
#Particle one velocities
vx_particle1=vx[0,:]
vy_particle1=vy[0,:]
vz_particle1=vz[0,:]

#Particle two positions
x_particle2=x[1,:]
y_particle2=y[1,:]
z_particle2=z[1,:]

#Particle two velocities
vx_particle2=vx[1,:]
vy_particle2=vy[1,:]
vz_particle2=vz[1,:]

#Plot both particles in xy plane
plt.plot(x_particle1,y_particle1,label="Particle 1")

plt.plot(x_particle2,y_particle2,label="Particle 2")

plt.scatter(x_particle2[0],y_particle2[0],c="g")
plt.scatter(x_particle1[0],y_particle1[0],c="g",label="Initial position")
plt.scatter(x_particle2[-1],y_particle2[-1],c="r")
plt.scatter(x_particle1[-1],y_particle1[-1],c="r",label="End position")
plt.xlabel(r"x [$\mu$m]")
plt.ylabel(r"y [$\mu$m]")
plt.legend()
plt.savefig("../Fig/RK4_2_particles_xy_plane_on.pdf")
plt.show()

#Plot both particles phase space

plt.plot(x_particle1,vx_particle1,label="Particle 1")
plt.plot(x_particle2,vx_particle2,label="Particle 2")
plt.axis('equal')
plt.xlabel(r"x [$\mu$m]")
plt.ylabel(r"v_x [$\mu$m/$\mu$s]")
plt.legend()
plt.savefig("../Fig/RK4_2_particles_phase_x_on.pdf")

plt.close()

plt.plot(z_particle1,vz_particle1,label="Particle 1")
plt.plot(z_particle2,vz_particle2,label="Particle 2")
#plt.axis('equal')
plt.xlabel(r"z [$\mu$m]")
plt.ylabel(r"v_z [$\mu$m/$\mu$s]")
plt.legend()
plt.savefig("../Fig/RK4_2_particles_phase_z_on.pdf")
plt.close()

#make 3D plot of particles in (x,y,z) space
ax = plt.figure().add_subplot(projection='3d')

ax.plot(x_particle1, y_particle1, z_particle1, label='Particle 1')
ax.plot(x_particle2, y_particle2, z_particle2, label='Particle 2')

ax.scatter(x_particle2[0],y_particle2[0],z_particle2[0],c="g")
ax.scatter(x_particle1[0],y_particle1[0],z_particle1[0],c="g",label="Initial position")
ax.scatter(x_particle2[-1],y_particle2[-1],z_particle2[-1],c="r")
ax.scatter(x_particle1[-1],y_particle1[-1],z_particle1[-1],c="r",label="End position")

ax.set_xlabel(r"x [$\mu$m]")
ax.set_ylabel(r"y [$\mu$m]")
ax.set_zlabel(r"z [$\mu$m]")

ax.legend()
#plt.savefig("../Fig/RK4_2_particles_3D_off.pdf")
# plt.close()
plt.show()

