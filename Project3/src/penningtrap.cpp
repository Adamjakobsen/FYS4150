#include "PenningTrap.hpp"

// Constructor 
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;

    std::vector<Particle> particles_init; 
    particles = particles_init;
    // particles.at(0)
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r)
{
    double x = r(0);
    double y = r(1);
    double z = r(2);

    //define the potential given in the exercise
    double prefactor = V0/(d*d);
    
    //derive with respect to all three components to find e-field
    arma::vec E_field = { prefactor*x, prefactor*y, -prefactor*2.*z };
    return E_field;
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r)
{
//include B-field here with z-direction, given by B-field strength 
//basically define an arma::vector here for it 
    // return arma::vec:: external_B_field = {}
}

//Add particle to the trap 
void PenningTrap::add_particle(Particle p_in)
{
    
}

//Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
{

}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
    
}

