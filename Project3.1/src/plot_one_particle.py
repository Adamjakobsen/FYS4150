import numpy as np
import matplotlib.pyplot as plt
from Utils import *

N_list=[4000,8000,16000,32000]
n_particles=1
filenames=[f"RK4_{n_particles}_{n_steps}.txt" for n_steps in N_list ]

############### X(t),Y(t) and Z(t) with TK4 ################

for i in range(len(N_list)):
    n_steps=N_list[i]
    filename=filenames[i]
    utils=Utils(n_steps,50)
    t,x,y,z,vx,vy,vz=utils.get_data(filename)
    x_analytical,y_analytical,z_analytical=utils.analytical()
    plt.plot(t,x_analytical,label="Analytical")
    plt.plot(t,x,linestyle='dashed',label=f"Numerical n = {N_list[i]}")
    plt.legend()
    plt.xlabel(r"Time [$\mu$s]")
    plt.ylabel(r"x(t) [$\mu$m]")
    plt.savefig(f"../Fig/RK4_1_particle_x(t)_{n_steps}_steps.pdf")
    plt.close()


####Y(t)########

for i in range(len(N_list)):
    n_steps=N_list[i]
    filename=filenames[i]
    utils=Utils(n_steps,50)
    t,x,y,z,vx,vy,vz=utils.get_data(filename)
    x_analytical,y_analytical,z_analytical=utils.analytical()
    plt.plot(t,y_analytical,label="Analytical")
    plt.plot(t,y,linestyle='dashed',label=f"Numerical n = {N_list[i]}")
    plt.savefig(f"../Fig/RK4_1_particle_y(t)_{n_steps}_steps.pdf")
    plt.legend()
    plt.xlabel(r"Time [$\mu$s]")
    plt.ylabel(r"y(t) [$\mu$m]")
    plt.close()


####Z(t)########

for i in range(len(N_list)):
    n_steps=N_list[i]
    filename=filenames[i]
    utils=Utils(n_steps,50)
    t,x,y,z,vx,vy,vz=utils.get_data(filename)
    x_analytical,y_analytical,z_analytical=utils.analytical()
    plt.plot(t,z_analytical,label="Analytical")
    plt.plot(t,z,linestyle='dashed',label=f"Numerical n = {N_list[i]}")
    plt.xlabel(r"Time [$\mu$s]")
    plt.ylabel(r"z(t) [$\mu$m]")
    plt.legend()
    plt.savefig(f"../Fig/RK4_1_particle_z(t)_{n_steps}_steps.pdf")
    plt.close()





    
########Relative error########

#RK4
filenames=[f"RK4_{n_particles}_{n_steps}.txt" for n_steps in N_list ]

for i in range(len(N_list)):
    n_steps=N_list[i]
    filename=filenames[i]
    utils=Utils(n_steps,50)
    t,x,y,z,vx,vy,vz=utils.get_data(filename)
    x_analytical,y_analytical,z_analytical=utils.analytical()
    r_analytical=np.sqrt(x_analytical**2 + y_analytical**2 + z_analytical**2)
    r=np.sqrt(x**2 + y**2 + z**2)
    plt.plot(t,utils.relative_error(r_analytical,r),label=f"{n_steps} steps")
    
plt.xlabel(r"Time [$\mu$s]")
plt.ylabel("Relative error")
plt.yscale("log")
plt.legend()
plt.savefig("../Fig/RK4_relative_error.pdf")
plt.close()


#Euler
filenames=[f"Euler_{n_particles}_{n_steps}.txt" for n_steps in N_list ]
for i in range(len(N_list)):
    n_steps=N_list[i]
    filename=filenames[i]
    utils=Utils(n_steps,50)
    t,x,y,z,vx,vy,vz=utils.get_data(filename)
    x_analytical,y_analytical,z_analytical=utils.analytical()
    r_analytical=np.sqrt(x_analytical**2 + y_analytical**2 + z_analytical**2)
    r=np.sqrt(x**2 + y**2 + z**2)
    plt.plot(t,utils.relative_error(r_analytical,r),label=f"{n_steps} steps")
plt.xlabel(r"Time [$\mu$s]")
plt.ylabel("Relative error")
plt.yscale("log")
plt.legend()
plt.savefig("../Fig/Euler_relative_error.pdf")
plt.close()