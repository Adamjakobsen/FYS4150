#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "./include/Particle.hpp"
#include "./include/PenningTrap.hpp"

void prettyprint(arma::vec armavec, std::string title, std::vector<std::string> labels);
void RK4_data(int n_particles,int n_timesteps,double total_time, std::string inter, double frequency, double amplitude);
void Euler_data(int n_particles,int n_timesteps,double total_time);

///////// MAKING A PENNING TRAP AND MAKING IT GLOBALLY ACCESSIBLE /////////

//define your physical system's variables here:
//use atomic mass units, micrometers, microseconds and element.charge e 
double B0 = 9.65 * pow(10,1);
double V0 = 2.41 * pow(10, 6);  
double d = 500; //micrometers
double q = 1;//.602 * pow(10,-19); //double check with the A-team
double m = 40.078; //calcium ion's atomic mass in atomic mass units



PenningTrap PT = PenningTrap(B0, V0, d);



// main function takes input arguments from the command line: Method, n_timesteps, total_time
int main(int argc, char*argv[] )
{   

    //If the user doesn't input the correct number of arguments, print an error message and exit
    if (argc != 6)
    {
        std::cout << "Error: Incorrect number of arguments. Please input the method, number of timesteps, and total time, perturbation on/off and number of particles." << std::endl;
        return 1;
    }
    
    //Defining number of particles from input
    int n_particles = atoi(argv[5]);

    for (int i=0; i<n_particles; i++)
    {
        //Defining positions and velocities from input
        arma::vec r = arma::vec(3).randn()*0.1*d; //micrometers
        std::cout << "r = " << r << std::endl;
        arma::vec v = arma::vec(3).randn()*0.1*d; //micrometers per microsec

        //Defining particles and adding them to the Penning Trap
        Particle p = Particle(q, m, r, v);
        
        PT.add_particle(p);
        PT.particles.at(i).time=0;
    }
    //Option to turn interactions on and off
    std::string inter=argv[4];
    std::cout << inter <<"|" << "Type: " << typeid(inter).name() << std::endl;
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
    //Turning on perturbation
    PT.perturbation=true;
    //Setting frequency
    double frequency = 0.6;
    PT.omega=frequency;
    //Setting amplitude
    double amplitude = 0.1;
    PT.f=amplitude;


    // Define time step and number of time steps
    double total_time= atoi(argv[3]);
    int N = atoi(argv[2]);
    double dt = total_time/N;
    // If methis is Rk4, call RK4_data
    if (std::string(argv[1]) == "RK4")
    {
        RK4_data(n_particles,N,total_time,inter,frequency,amplitude);
    }
    // If method is Euler, call Euler_data
    else if (std::string(argv[1]) == "Euler")
    {
        Euler_data(n_particles,N,total_time);
    }
    // If method is not Euler or RK4, print error message
    else
    {
        std::cout << "Error: Method must be either Euler or RK4" << std::endl;
    }






    return 0;
}

//Function for producing data with RK4

void RK4_data(int n_particles,int n_timesteps,double total_time,std::string inter, double frequency, double amplitude)
{
    // Define time step and number of time steps
    double dt = total_time/n_timesteps;
    std::cout << "dt = " << dt << "total time:" << total_time << "timesteps: "<< n_timesteps<<std::endl;
    
    // Use rk4 to evolve and write to file

    // Format parameters
    std::ofstream outfile;
	int width = 25;
	int prec = 15;
    // make path to output file with the method, number of timesteps, and total time
    std::string path_str = "./Data/RK4_perturbation_"+std::to_string(n_particles)+"_" + std::to_string(n_timesteps) +"_" + inter +std::to_string(frequency)+"_" +std::to_string(amplitude)+"_" + ".txt";

    //const char *path=path_str;
    outfile.open(path_str);

    double time =0.0;

    for (int i = 0; i < n_timesteps +1; i++) 
    {       

            PT.evolve_RK4(dt);
            PT.time=time;
            
            time=dt*i;
            std::cout << "======= IN RK4 =======" << std::endl;
            std::cout << "time PT: " << PT.time << std::endl;
            std::cout << "time RK4: " << time << std::endl;
        for (int j = 0; j < n_particles; j++)
                    {
                
                    outfile << 
                    std::setw(width) << std::setprecision(prec) <<PT.particles.at(j).r.at(0) << 
                    std::setw(width) << std::setprecision(prec) << PT.particles.at(j).r.at(1) << 
                    std::setw(width) << std::setprecision(prec)<< PT.particles.at(j).r.at(2) <<  
                    
                    std::setw(width) << std::setprecision(prec) << PT.particles.at(j).v.at(0) << 
                    std::setw(width) << std::setprecision(prec) << PT.particles.at(j).v.at(1) << 
                    std::setw(width) << std::setprecision(prec) << PT.particles.at(j).v.at(2) <<
                    std::setw(width) << std::setprecision(prec) << j << std::endl;
                    }
            
        
            }
        
        
    outfile.close();



        

}

//Function for producing data with Euler

void Euler_data(int n_particles,int n_timesteps,double total_time)
{
    // Define time step and number of time steps
    double dt = total_time/n_timesteps;
    
    // Use rk4 to evolve and write to file

    // Format parameters
    std::ofstream outfile;
    int width = 25;
    int prec = 15;

        // make path to output file with the method, number of timesteps, and total time
    std::string path_str = "./Data/Euler_" +std::to_string(n_particles)+"_"+ std::to_string(n_timesteps) + ".txt";
    outfile.open(path_str);

        for (int i = 0; i < n_timesteps +1; i++) 
    {   
            PT.evolve_RK4(dt);
            
        for (int j = 0; j < n_particles; j++)
        {
            
                outfile << 
                std::setw(width) << std::setprecision(prec) <<PT.particles.at(j).r.at(0) << 
                std::setw(width) << std::setprecision(prec) << PT.particles.at(j).r.at(1) << 
                std::setw(width) << std::setprecision(prec)<< PT.particles.at(j).r.at(2) <<  
                
                std::setw(width) << std::setprecision(prec) << PT.particles.at(j).v.at(0) << 
                std::setw(width) << std::setprecision(prec) << PT.particles.at(j).v.at(1) << 
                std::setw(width) << std::setprecision(prec) << PT.particles.at(j).v.at(2) <<  
                std::endl;
                
            
            }
    }
    outfile.close();


}

//nice function to print stuff out : enables print of xyz position and velocity
void prettyprint(arma::vec armavec, std::string title, std::vector<std::string> labels)
{
    if (labels.size() != armavec.n_elem)
    {
        std::cout << "Size mismatch." << std::endl;
        return;
    }

    std::cout << title << std::endl;

    for (int i = 0; i < labels.size(); ++i)
        std::cout << "  " << labels.at(i) << ": " << armavec.at(i) << std::endl;
}
