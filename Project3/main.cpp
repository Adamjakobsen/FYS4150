#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "./include/Particle.hpp"
#include "./include/PenningTrap.hpp"

void prettyprint(arma::vec armavec, std::string title, std::vector<std::string> labels);

int main()
{
    //define your positions and velocities here 
    //Anders defines his vectors using strings but I don't like it so I changed it to this 
    //It basically does the same thing so it should not matter hehe 
    arma::vec r1 = arma::vec(std::vector<double> { 20, 0, 20 }); //micrometers 
    arma::vec v1 = arma::vec(std::vector<double> { 0, 25, 0 }); //micrometers per microsec
    

    arma::vec r2 = arma::vec(std::vector<double> { 25, 25, 0 }); //micrometers 
    arma::vec v2 = arma::vec(std::vector<double> { 0, 40, 5 }); //micrometers per microsec
    

    //define your physical system's variables here:
    //use atomic mass units, micrometers, microseconds and element.charge e 
    double B0 = 9.65 * pow(10,1);
    double V0 = 2.41 * pow(10, 6);  
    double d = 500; //micrometers
    double q = 1;//.602 * pow(10,-19); //double check with the A-team
    double m = 40.078; //calcium ion's atomic mass in atomic mass units


    //remember to define correct elementary charge value here!!
    //here, I create object of type "Particle" using our class, which I call "particle1", and assign attributes to it such as 
    //charge, mass, position and velocity!
    Particle particle1 = Particle(q, m, r1, v1);
    //I do the same thing here, just for the second particle 
    Particle particle2 = Particle(q, m, r2, v2);
    PenningTrap PT = PenningTrap(B0, V0, d);


    // particle1.print_attributes();
    // particle2.print_attributes();
    PT.add_particle(particle1);
    PT.add_particle(particle2);


    // arma::vec ext_electric_field = PT.external_E_field(particle.r);
    // arma::vec ext_magnetic_field = PT.external_B_field(particle.r);

    // Define time step and number of time steps
    double dt = 0.01; // microseconds
    int N = 100000; // number of time steps
    int n_particles = PT.particles.size();
    // Use rk4 to evolve and write to file
    // Format parameters

    std::ofstream outfile;
	int width = 18;
	int prec = 10;

    
    outfile.open("positions_rk4.txt");

    for (int j = 0; j < n_particles; j++)
    {
    for (int i = 0; i < N; i++) 
        {
        
        PT.evolve_RK4(dt);
        
            outfile << 
            std::setw(width) << std::setprecision(prec) <<PT.particles.at(j).r.at(0) << 
            std::setw(width) << std::setprecision(prec) << PT.particles.at(j).r.at(1) << 
            std::setw(width) << std::setprecision(prec)<< PT.particles.at(j).r.at(2) <<  std::endl;
        
        }
    }
    
    outfile.close();
    
    outfile.open("positions_Euler.txt");
    dt=0.001;
    N=10000;
    

    for (int j = 0; j < n_particles; j++)
    {
    //outfile << "particle " << j << std::endl;
    for (int i = 0; i < N; i++) 
    {
        
        PT.evolve_forward_Euler(dt);
        
            outfile << 
            std::setw(width) << std::setprecision(prec) <<PT.particles.at(j).r.at(0) << 
            std::setw(width) << std::setprecision(prec) << PT.particles.at(j).r.at(1) << 
            std::setw(width) << std::setprecision(prec)<< PT.particles.at(j).r.at(2) <<  std::endl;
        
    }
    }
    outfile.close();
    
    
    for (int i = 0; i < PT.particles.size(); ++i)
    {
        PT.particles.at(i).print_attributes();
        std::cout << "Particle cpp has run" << std::endl;
        prettyprint(PT.particles.at(i).r, "Particle " + std::to_string(i + 1) + "'s positions", std::vector<std::string>{"x", "y", "z"});
        prettyprint(PT.particles.at(i).v, "Particle " + std::to_string(i + 1) + "'s velocities", std::vector<std::string>{"x", "y", "z"});
    }

    // std::cout << "External E-field" << "\n" << electric_field;
    // prettyprint(ext_electric_field, "External E-field", std::vector<std::string> { "x", "y", "z" });
    // prettyprint(ext_magnetic_field, "External B-field", std::vector<std::string> { "x", "y", "z" });
    
    arma::vec force_ij = PT.force_particle(0,1);
    prettyprint(force_ij, "Force_ij", std::vector<std::string> { "x", "y", "z" });

    // evolve with rk4 from PenningTrap class with dt = 0.01 microseconds
    
    
    

    

    return 0;
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
