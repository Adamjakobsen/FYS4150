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

    
    for (int i = 0; i < particles.size(); i++){
        
        double q=particles.at(i).q;
        double m=particles.at(i).m;
        //Lets first test with one particle in the xternal field
        
        arma::vec r = particles.at(i).r;
        arma::vec v = particles.at(i).v;
        arma::vec F = total_force(i);
        arma::vec k1= (F/m) * dt;

        particles.at(i).v = v + k1*1/2.*dt;
        particles.at(i).r = r + k1*1/2.*dt;
        r = particles.at(i).r;
        v = particles.at(i).v;

        F = total_force(i);
        arma::vec k2 = (F/m)*dt;

        particles.at(i).v = v + k2*1/2.*dt;
        particles.at(i).r = r + k2*1/2.*dt;
        r = particles.at(i).r;
        v = particles.at(i).v;

        F = total_force(i);
        arma::vec k3 = (F/m)*dt;

        particles.at(i).v = v + k3*dt;
        particles.at(i).r = r + k3*dt;
        r = particles.at(i).r;
        v = particles.at(i).v;

        F = total_force(i);
        arma::vec k4 = (F/m)*dt;

        particles.at(i).v = v + 1/6.*dt*(k1 + 2*k2 + 2*k3 + k4);
        particles.at(i).r = r + 1/6.*dt*(k1 + 2*k2 + 2*k3 + k4);
    }

}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{   
    for (int i = 0; i < particles.size(); i++)
    {
        double q=particles.at(i).q;
        double m=particles.at(i).m;
        

        arma::vec r = particles.at(i).r;
        arma::vec v = particles.at(i).v;
        arma::vec F = total_force(i);
        arma::vec a = F/m;


        //Evolve using foward euler
        particles.at(i).r=r+v*dt;
        particles.at(i).v= v + a*dt;
    }
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

    //define empty 3 dimentional vector for our force
    arma::vec force_external = arma::vec(3, arma::fill::zeros);

    //define some constants and other variables that will be used 
    double q = particles.at(i).q;
    double m = particles.at(i).m;
    arma::vec r = particles.at(i).r;
    arma::vec v = particles.at(i).v;

    //calculate the force on particle_i from the external fields
    force_external = q*external_E_field(r) + q*arma::cross(v,external_B_field(r) );
    return force_external;

}

 // The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i)
{
    
        //define empty 3 dimentional vector for our force
        arma::vec force_particles = arma::vec(3, arma::fill::zeros);
    
        //define some constants and other variables that will be used 
        double q = particles.at(i).q;
        double m = particles.at(i).m;
        arma::vec r = particles.at(i).r;
        arma::vec v = particles.at(i).v;
    
        int n_particles = particles.size();
        //calculate the force on particle_i from all the other particles
        for (int j = 0; j < n_particles; j++){
            if (j != i){
                force_particles += force_particle(i,j);
            }
        }
        return force_particles;

}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i)
{
    arma::vec force_total = arma::vec(3, arma::fill::zeros);
    force_total = total_force_external(i) + total_force_particles(i);
    return force_total;

}


