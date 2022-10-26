import numpy as np
import matplotlib.pyplot as plt
from Utils import *

N_list=[4000, 8000, 16000, 32000]
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
    plt.plot(t,x, linestyle= "dashed", label=f"Numerical n = {N_list[i]}")
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
    plt.plot(t,y, linestyle='dashed',label=f"Numerical n = {N_list[i]}")
    plt.legend()
    plt.xlabel(r"Time [$\mu$s]")
    plt.ylabel(r"y(t) [$\mu$m]")
    plt.savefig(f"../Fig/RK4_1_particle_y(t)_{n_steps}_steps.pdf")
    plt.close()


####Z(t)########

for i in range(len(N_list)):
    n_steps=N_list[i]
    filename=filenames[i]
    utils=Utils(n_steps,50)
    t,x,y,z,vx,vy,vz=utils.get_data(filename)
    x_analytical,y_analytical,z_analytical=utils.analytical()
    plt.plot(t,z_analytical,label="Analytical")
    plt.plot(t,z, linestyle='dashed', label=f"Numerical n = {N_list[i]}")
    #plt.plot(t,z, label=f"Numerical n = {N_list[i]}")
    plt.xlabel(r"Time [$\mu$s]")
    plt.xlabel(r"Time [$\mu$s]")
    plt.ylabel(r"z(t) [$\mu$m]")
    plt.legend()
    plt.savefig(f"../Fig/RK4_1_particle_z(t)_{n_steps}_steps.pdf")
    plt.close()


# *** added pahse space in z plot
#########phase space in z ###############  
for i in range(len(N_list)):
    n_steps=N_list[i]
    filename=filenames[i]
    utils=Utils(n_steps,50)
    t,x,y,z,vx,vy,vz=utils.get_data(filename)
    x_analytical,y_analytical,z_analytical =utils.analytical()
    #plt.plot(vz,z_analytical,label="Analytical")
    plt.plot(vz,z,label=f"Numerical n = {N_list[i]}")
    plt.xlabel(r"Time [$\mu$s]")
    plt.ylabel(r"z(t) [$\mu$m]")
    plt.legend()
    plt.savefig(f"../Fig/RK4_1_particle_z(t)vz_{n_steps}_steps.pdf")
    plt.close()


"""   
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
"""

#########RELATIVE ERROR#############
    #RK4
filenames=[f"RK4_{n_particles}_{n_steps}.txt" for n_steps in N_list ]
for i in range(len(N_list)):
    n_steps=N_list[i]
    filename=filenames[i]
    print(filename)
    utils=Utils(n_steps,50)
    t,x,y,z,vx,vy,vz=utils.get_data(filename)
    x_analytical,y_analytical,z_analytical=utils.analytical()
    r_analytical=np.array([x_analytical, y_analytical,  z_analytical])
    r_num =np.array([x, y, z])
    # r=np.sqrt(x**2 + y**2 + z**2)
    plt.plot(t,utils.relative_error(r_analytical, r_num),label=f"{n_steps} steps")
plt.xlabel(r"Time [$\mu$s]")
plt.ylabel("Relative error")
plt.yscale("log")
plt.legend()
plt.savefig("../Fig/RK4_relative_error_.pdf")
plt.show()

    #Euler
filenames=[f"Euler_{n_particles}_{n_steps}.txt" for n_steps in N_list ]
for i in range(len(N_list)):
    n_steps=N_list[i]
    filename=filenames[i]
    print(filename)
    utils=Utils(n_steps,50)
    t,x,y,z,vx,vy,vz=utils.get_data(filename)
    x_analytical,y_analytical,z_analytical=utils.analytical()
    r_analytical=np.array([x_analytical, y_analytical,  z_analytical])
    r_num =np.array([x, y, z])
    # r=np.sqrt(x**2 + y**2 + z**2)
    plt.plot(t,utils.relative_error(r_analytical, r_num),label=f"{n_steps} steps")
plt.xlabel(r"Time [$\mu$s]")
plt.ylabel("Relative error")
plt.yscale("log")
plt.legend()
plt.savefig("../Fig/Euler_relative_error_.pdf")
plt.show()

#######Error convergence rate RK4########
filenames=[f"RK4_{n_particles}_{n_steps}.txt" for n_steps in N_list ]
h_list=[1/n_steps for n_steps in N_list]
max_error_list=[]
for i in range(len(N_list)):
    n_steps=N_list[i]
    filename=filenames[i]
    utils=Utils(n_steps,50)
    t,x,y,z,vx,vy,vz=utils.get_data(filename)
    x_analytical,y_analytical,z_analytical=utils.analytical()
    r_analytical=np.sqrt(x_analytical**2 + y_analytical**2 + z_analytical**2)
    r=np.sqrt(x**2 + y**2 + z**2)
    #Estimating the error convergence rate
    max_error_list.append(np.max(np.abs(r_analytical-r)))
r_err = 0
print(f'max errors: {max_error_list}')
print(f'h values: {h_list}')

for i in range(1,4):
    r_err += (1/3)*(np.log10(max_error_list[i]/max_error_list[i-1])/np.log10(h_list[i]/h_list[i-1]))

print(f"Error convergence rate RK4: {r_err}")


#######Error convergence rate Forward Euler########
filenames=[f"Euler_{n_particles}_{n_steps}.txt" for n_steps in N_list ]
h_list=[1/n_steps for n_steps in N_list]
max_error_list=[]
for i in range(len(N_list)):
    n_steps=N_list[i]
    filename=filenames[i]
    utils=Utils(n_steps,50)
    x_analytical,y_analytical,z_analytical=utils.analytical()
    t,x,y,z,vx,vy,vz=utils.get_data(filename)
    r_analytical=np.sqrt(x_analytical**2 + y_analytical**2 + z_analytical**2)
    r=np.sqrt(x**2 + y**2 + z**2)
    #Estimating the error convergence rate
    max_error_list.append(np.max(np.abs(r_analytical-r)))
r_err = 0
for i in range(1,4):
    r_err += (1/3)*(np.log10(max_error_list[i]/max_error_list[i-1])/np.log10(h_list[i]/h_list[i-1]) )

print(f"Error convergence rate Forward Euler: {r_err}")
