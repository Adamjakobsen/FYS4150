#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "./include/Particle.hpp"
#include "./include/PenningTrap.hpp"



///////// MAKING A PENNING TRAP AND MAKING IT GLOBALLY ACCESSIBLE /////////

//define your physical system's variables here:
//use atomic mass units, micrometers, microseconds and element.charge e 
double B0 = 9.65 * pow(10,1);
double V0 = 2.41 * pow(10, 6);  
double d = 500; //micrometers
double q = 1;//.602 * pow(10,-19); //double check with the A-team
double m = 40.078; //calcium ion's atomic mass in atomic mass units







// main function takes input arguments from the command line: Method, n_timesteps, total_time
int main(int argc, char*argv[] )
{   
    //Random seed
    arma::arma_rng::set_seed(0);

    //If the user doesn't input the correct number of arguments, print an error message and exit
    if (argc != 6)
    {
        std::cout << "Error: Incorrect number of arguments. Please input the method, number of timesteps, and total time, perturbation on/off and number of particles." << std::endl;
        return 1;
    }
    
    //Defining number of particles from input
    int n_particles = atoi(argv[5]);

    
    
    // Define time step and number of time steps
    double total_time= atoi(argv[3]);
    int N = atoi(argv[2]);
    double dt = total_time/N;
    //Defining frequency linspace 
    arma::vec frequency_list = arma::linspace(1.0,1.8,80);
    //Definin amplitude values 
    arma::vec amplitude_list = {0.4};

    //Option to turn interactions on and off
    std::string inter=argv[4];

    for (int f_i=0;f_i< 3;f_i++)

    {   std::cout << "Amplitude: " << amplitude_list(f_i) << std::endl;
        // Format parameters
        std::ofstream outfile;
        int width = 25;
        int prec = 15;
        // make path to output file with the method, number of timesteps, and total time
        std::string path_str = "./Data/Fine_perturbation_"+std::to_string(n_particles)+"_" + std::to_string(N) +"_" + inter+std::to_string(amplitude_list.at(f_i)) +".txt";

        //const char *path=path_str;
        outfile.open(path_str);
        for (int omega_i=0; omega_i<100; omega_i++)
        {   std::cout << "Freqquency:"<< frequency_list(omega_i)  << std::endl;
            
            PenningTrap PT = PenningTrap(B0, V0, d);
            //Turning on perturbation
            PT.perturbation=true;
            for (int i=0; i<n_particles; i++)
        {
            //Defining positions and velocities from input
            arma::vec r = arma::vec(3).randn()*0.1*d; //micrometers
            //std::cout << "r = " << r << std::endl;
            arma::vec v = arma::vec(3).randn()*0.1*d; //micrometers per microsec

            //Defining particles and adding them to the Penning Trap
            Particle p = Particle(q, m, r, v);

            PT.add_particle(p);
            PT.particles.at(i).time=0;
        }
    if (inter == "on")
    {
        PT.toggle_interaction(true);
    }
    else if (inter == "off")
    {
        PT.toggle_interaction(false);
    }
    else
    {
        std::cout << "Error: Please input 'on' or 'off' for particle-particle interaction." << std::endl;
        return 1;
    }
            //initialize time
            PT.time = 0.0;
            //std::cout << "frequency = " << frequency_list.at(omega_i) << std::endl;
            
            
            PT.omega=frequency_list.at(omega_i);
            //std::cout << "IN loopsies omega = " << PT.omega << std::endl;
            //Setting amplitude

            PT.f=amplitude_list.at(f_i);                     
                       
            // Use rk4 to evolve and write to file


            double time =0.0;
            int particles_out=0;

            for (int i = 0; i < N +1; i++) 
            {       

                    PT.evolve_RK4(dt);
                    PT.time=time;

                    time=dt*i;
    
            }
            for (int j=0;j<n_particles;j++)
            {   
                if (PT.particles.at(j).q==0.0)
                {
                    particles_out+=1;

                }
            }

            outfile << frequency_list.at(omega_i)<<"  " << particles_out << std::endl;
                
        
            
        

         }
         outfile.close();
    }



    return 0;
}
