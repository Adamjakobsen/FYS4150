#include <iostream>
#include <cmath>

#include <armadillo>
#include "Particle.hpp"
#include "PenningTrap.hpp"

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
    double q = 1.602 * pow(10,-19); //double check with the A-team


    //remember to define correct elementary charge value here!!
    //here, I create object of type "Particle" using our class, which I call "particle1", and assign attributes to it such as 
    //charge, mass, position and velocity!
    Particle particle1 = Particle(1.6, 2.5, r1, v1);
    //I do the same thing here, just for the second particle 
    Particle particle2 = Particle(1.6, 2.5, r2, v2);
    PenningTrap PT = PenningTrap(B0, V0, d);


    // particle1.print_attributes();
    // particle2.print_attributes();
    PT.add_particle(particle1);
    PT.add_particle(particle2);

    // arma::vec ext_electric_field = PT.external_E_field(particle.r);
    // arma::vec ext_magnetic_field = PT.external_B_field(particle.r);

    for (int i = 0; i < PT.particles.size(); ++i)
    {
        PT.particles.at(i).print_attributes();
        prettyprint(PT.particles.at(i).r, "Particle " + std::to_string(i + 1) + "'s positions", std::vector<std::string>{"x", "y", "z"});
        prettyprint(PT.particles.at(i).v, "Particle " + std::to_string(i + 1) + "'s velocities", std::vector<std::string>{"x", "y", "z"});
    }

    // std::cout << "External E-field" << "\n" << electric_field;
    // prettyprint(ext_electric_field, "External E-field", std::vector<std::string> { "x", "y", "z" });
    // prettyprint(ext_magnetic_field, "External B-field", std::vector<std::string> { "x", "y", "z" });
 

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
