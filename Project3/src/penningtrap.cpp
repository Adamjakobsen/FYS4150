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
     std::cout << "External E-field is:" << arma::vec(std::vector<double> { prefactor_value * x, prefactor_value * y, -prefactor_value * 2. * z })<< std::endl;
    return arma::vec(std::vector<double> { prefactor_value * x, prefactor_value * y, -prefactor_value * 2. * z });
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r)
{
//include B-field here with z-direction, given by B-field strength
    std::cout << "External B-field is:" << arma::vec(std::vector<double> { 0, 0, B0})<< std::endl;
    return arma::vec(std::vector<double> { 0, 0, B0});
}

//Evolve the system one time step (dt) using Runge-Kutta 4th order

void PenningTrap::evolve_RK4(double dt)
{
     
// created a copy of particles to avoid changing the original vector
std::vector<Particle> particles_copy = particles;

    double q=particles.at(0).q;
    double m=particles.at(0).m;
 
    arma::vec k1_r = arma::vec(3);
    arma::vec k1_v = arma::vec(3);
    arma::vec k2_r = arma::vec(3);
    arma::vec k2_v = arma::vec(3);
    arma::vec k3_r = arma::vec(3);
    arma::vec k3_v = arma::vec(3);
    arma::vec k4_r = arma::vec(3);
    arma::vec k4_v = arma::vec(3);
    //Total force on particle i is total_force
    for (int i=0; i<particles.size(); i++)
    {
        k1_r = particles.at(i).v*dt;
        k1_v = 1/m * total_force(i)*dt;

        particles.at(i).r = particles.at(i).r + 0.5*k1_r(i)*dt;
        particles.at(i).v = particles.at(i).v + 0.5*k1_v(i)*dt;

        k2_r = dt * particles.at(i).r + 0.5*k1_r(i);
        k2_v = dt * total_force(i)/m;

        particles.at(i).r = particles.at(i).r + 0.5*k2_r(i)*dt;
        particles.at(i).v = particles.at(i).v + 0.5*k2_v(i)*dt;

        k3_r = dt * particles.at(i).r + 0.5*k2_r(i);
        k3_v = dt * total_force(i)/m;

        particles.at(i).r = particles.at(i).r + k3_r(i)*dt;
        particles.at(i).v = particles.at(i).v + k3_v(i)*dt;

        k4_r = dt * particles.at(i).r + k3_r(i);
        k4_v = dt * total_force(i)/m;

        particles.at(i).r = particles_copy.at(i).r + 1/6.*(k1_r + 2*k2_r + 2*k3_r + k4_r );
        particles.at(i).v = particles_copy.at(i).v + 1/6.*(k1_v + 2*k2_v + 2*k3_v + k4_v );
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

    std::cout << "Total external force is:" << force_external << std::endl;
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
         std::cout << "Total force from particles is:" << force_particles << std::endl;
        return force_particles;

}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i)
{
    arma::vec force_total = arma::vec(3, arma::fill::zeros);
    force_total = total_force_external(i) + total_force_particles(i);
    std::cout << "The total force is" << force_total << std::endl; 
    return force_total;

}


