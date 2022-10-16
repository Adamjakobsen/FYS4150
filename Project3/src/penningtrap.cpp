#include "PenningTrap.hpp"

// Constructor 
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;

    std::vector<Particle> particles_init; 
    particles = particles_init;

}

//Add particle to the trap 
void PenningTrap::add_particle(Particle p_in)
{
    particles.push_back(p_in);
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r)
{
    double x = r(0);
    double y = r(1);
    double z = r(2);

    //define the potential given in the exercise
    double prefactor_value = V0/(d*d); //not sure if we even need this definition if we already know what the prefactor value is
    
    //derive with respect to all three components to find e-field 
    return arma::vec(std::vector<double> { prefactor_value * x, prefactor_value * y, -prefactor_value * 2. * z });
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r)
{
//include B-field here with z-direction, given by B-field strength 
    return arma::vec(std::vector<double> { 0, 0, B0});
}

//Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
{

}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
    
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j)
{

    //define empty 3 dimentional vector for our force
    arma::vec force_ij = arma::vec(3, arma::fill::zeros);

    //define some constants and other variables that will be used 
    double k_e = 1.38935333;
    double charges = particles.at(i).q * particles.at(j).q;

    //calculate the distance from particle_i to particle_j and find norm of the force
    arma::vec distance_ij = particles.at(i).r - particles.at(j).r;
    double distance_ij_norm = arma::norm(distance_ij, 1);

    //calculate the force on particle_i from particle_j
    force_ij = (k_e*charges * distance_ij ) / (pow(distance_ij_norm,3));
    return force_ij;
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i)
{

}

 // The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i)
{

}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i)
{

}

