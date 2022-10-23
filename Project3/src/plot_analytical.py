from re import M
import numpy as np
import matplotlib.pyplot as plt



class Utils():

    def __init__(self, n_timesteps, end_time,n_paraticles=1,method="RK4"):
        """
        Class with utility functions for plotting and calculating errors

        Args:
            method (str, optional): Numerical method to use. Defaults to "RK4".
            n_timesteps (int): Number of timesteps to use
            end_time (float): End time of the simulation
            n_paraticles (int, optional): Number of particles simulated. Defaults to 1.       
       
        """
        #constants:
        self.v0 = 25
        self.B0 = 9.65 * 10
        self.V0 = 2.41*10**6 
        self.m = 40.078 
        self.d = 500 
        self.q = 1 
        self.w_z = np.sqrt(2*self.q*self.V0 / (self.m * self.d**2))
        self.w0 = self.q*self.B0 / self.m

        #Args:
        self.n_timesteps=n_timesteps
        self.end_time=end_time
        self.n_particles=n_paraticles
        self.method=method

        
    
    def analytical(self):


        # initial conditions
        z0 = 20
        y0 = 0
        x0 = 20

        t = np.linspace(0, self.end_time, self.n_timesteps+1)

        w_p = (self.w0 + np.sqrt(self.w0**2 - 2*self.w_z**2)) / 2
        w_m = (self.w0 - np.sqrt(self.w0**2 - 2*self.w_z**2)) / 2

        A_p = (self.v0 + w_m*x0)/(w_m - w_p)
        A_m = -(self.v0 + w_p*x0)/(w_m - w_p)

        self.x_analytical = A_p*np.cos(-w_p*t) + A_m*np.cos(-w_m*t)
        self.y_analytical = A_p*np.sin(-w_p*t) + A_m*np.sin(-w_m*t)
        self.z_analytical = z0 * np.cos(self.w_z * t)

        return self.x_analytical, self.y_analytical, self.z_analytical
    
    def relative_error(self,analytical,numerical):
        #Analytical solution passes through zero several times, so we shift both the numerical and analytical solution in order to avoid dividing by zero

        max_diff=np.max(analytical) - np.min(analytical)
        analytical += max_diff + 1
        numerical += max_diff + 1

        return np.abs((analytical-numerical)/analytical)

    def get_data(self,filename):
        data = np.loadtxt(f'../Data/{filename}')
        n=self.n_timesteps + 1
        t=np.linspace(0,self.end_time,n)

        if self.n_particles==1:
            x = data[:,0]
            y = data[:,1]
            z = data[:,2]
            
            vx = data[:,3]
            vy = data[:,4]
            vz = data[:,5]
        
        else:
            x=np.zeros((self.n_particles,n))
            y=np.zeros((self.n_particles,n))
            z=np.zeros((self.n_particles,n))
            vx=np.zeros((self.n_particles,n))
            vy=np.zeros((self.n_particles,n))
            vz=np.zeros((self.n_particles,n))

            for i in range(self.n_particles):
                first = i*n
                last=i*n + self.n_timesteps

                x[i,:] = data[first:last,0]
                y[i,:] = data[first:last,1]
                z[i,:] = data[first:last,2]
                
                vx[i,:] = data[first:last,3]
                vy[i,:] = data[first:last,4]
                vz[i,:] = data[first:last,5]
        return t,x,y,z,vx,vy,vz
            



if __name__ == '__main__':

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
    plt.savefig("../Fig/Euler_relative_error_.pdf")
    plt.close()
    


    
    




