#include "../include/PenningTrap.hpp"

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
// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
{
    //PenningTrap PT = PenningTrap(B0, V0, d);

    

    
    double q=particles.at(0).q;
    double m=particles.at(0).m;
    //Lets first test with one particle in the xternal field
    
    arma::vec r = particles.at(0).r;
    arma::vec v = particles.at(0).v;
    arma::vec ext_force = q*external_E_field(r) + q*arma::cross(v,external_B_field(r) );
    arma::vec k1= (ext_force/m) * dt;

    particles.at(0).v = v + k1*1/2.*dt;
    particles.at(0).r = r + k1*1/2.*dt;
    r = particles.at(0).r;
    v = particles.at(0).v;

    arma::vec k2 = 1/m*( q*external_E_field(r) + q*arma::cross(v,external_B_field(r) ) )*dt;

    particles.at(0).v = v + k2*1/2.*dt;
    particles.at(0).r = r + k2*1/2.*dt;
    r = particles.at(0).r;
    v = particles.at(0).v;

    arma::vec k3 = 1/m*(q*external_E_field(r) + q*arma::cross(v,external_B_field(r) ))*dt;

    particles.at(0).v = v + k3*dt;
    particles.at(0).r = r + k3*dt;
    r = particles.at(0).r;
    v = particles.at(0).v;

    arma::vec k4 = 1/m*(q*external_E_field(r) + q*arma::cross(v,external_B_field(r) ))*dt;

    particles.at(0).v = v + 1/6.*dt*(k1 + 2*k2 + 2*k3 + k4);
    particles.at(0).r = r + 1/6.*dt*(k1 + 2*k2 + 2*k3 + k4);
    
    /* 
    arma::vec r = particles.at(0).r;
    arma::vec v = particles.at(0).v;

    arma::vec ext_force = q*external_E_field(r) + q*arma::cross(v,external_B_field(r) );
    arma::vec k_v1 = ext_force/m * dt;
    arma::vec k_r1 = v * dt;

    arma::vec k_v2 = (q*external_E_field(r + k_r1*1/2..) + q*arma::cross(v + k_v1*1/2..,external_B_field(r + k_r1*1/2..) ))*dt/m;
    arma::vec k_r2 = (v + k_v1*1/2..) * dt;

    arma::vec k_v3 = (q*external_E_field(r + k_r2*1/2..) + q*arma::cross(v + k_v2*1/2..,external_B_field(r + k_r2*1/2..) ))*dt/m;
    arma::vec k_r3 = (v + k_v2*1/2..) * dt;

    arma::vec k_v4 = (q*external_E_field(r + k_r3) + q*arma::cross(v + k_v3,external_B_field(r + k_r3) ))*dt/m;
    arma::vec k_r4 = (v + k_v3) * dt;

    particles.at(0).v = v + 1/6.*dt*(k_v1 + 2.*k_v2 + 2.*k_v3 + k_v4);
    particles.at(0).r = r + 1/6.*dt*(k_r1 + 2.*k_r2 + 2.*k_r3 + k_r4);
    */

    /* 
    //define the k's for the RK4 method
    double kr_1 = dt * r;
    double kv_1 = dt * v;

    double mid_r = r + 0.5 * kr_1;
    double mid_v = v + 0.5 * kv_1;

    double kr_2 = dt * mid_r;
    double kv_2 = dt * mid_v;

    double mid_r2 = r + 0.5 * kr_2;
    double mid_v2 = v + 0.5 * kv_2;

    double kr_3 = dt * mid_r2;
    double kv_3 = dt * mid_v2;

    double r3 = r + kr_3;
    double v3 = v + kv_3;

    double kr_4 = dt * r3;
    double kv_4 = dt * v3;

    double position_updated = r + (1/6)*(kr_1 + 2*kr_2 + 2*kr_3 + kr_4);
    double velocity_updated = v + (1/6)*(kv_1 + 2*kv_2 + 2*kv_3 + kv_4);
     */ 
    //return position_updated, velocity_updated;

}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
    double q=particles.at(0).q;
    double m=particles.at(0).m;
    //Lets first test with one particle in the xternal field

    arma::vec r = particles.at(0).r;
    arma::vec v = particles.at(0).v;
    arma::vec ext_force = q*external_E_field(r) + q*arma::cross(v,external_B_field(r) );
    arma::vec a = ext_force/m;


    //Evolve using foward euler
    particles.at(0).r=r+v*dt;
    particles.at(0).v= v + a*dt;

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

