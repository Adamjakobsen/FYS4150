#include "../include/PenningTrap.hpp"
#include <cmath>

// Constructor 
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;

    std::vector<Particle> particles_init; 
    particles = particles_init;
    double omega;
    double f;
    double time;

}
//update time
void PenningTrap::update_time(double time_in)
{
    PenningTrap::time = time_in;
}
// Set frequency
void PenningTrap::set_frequency(double frequency_in)
{
    PenningTrap::omega = frequency_in;
}

void PenningTrap::set_amplitude(double f_in)
{
    PenningTrap::f = f_in;
}

//Add particle to the trap 
void PenningTrap::add_particle(Particle p_in)
{
    particles.push_back(p_in);
}
// Method for turning interaction on or off
void PenningTrap::toggle_interaction(bool interaction_in)
{
    PenningTrap::interaction = interaction_in;

}
//Method for turning on perturbation
void PenningTrap::toggle_perturbation(bool perturbation_in)
{
    PenningTrap::perturbation = perturbation_in;
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

// External E-field with perturbation
arma::vec PenningTrap::external_E_field_perturbed(arma::vec r, double t_in, double omega_in,double f_in)
{
    double x = r(0);
    double y = r(1);
    double z = r(2);

    //define the potential given in the exercise
    double prefactor_value = V0*(1. + f_in*cos(omega*t_in))/(d*d); //not sure if we even need this definition if we already know what the prefactor value is
    
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
     
// created a copy of particles to avoid changing the original vector
    std::vector<Particle> particles_copy = particles;

    
    double q=particles.at(0).q;
    double m=particles.at(0).m;
    
    double n_particles = particles.size();
    std::vector<arma::vec> k1_r(n_particles);
    std::vector<arma::vec> k1_v(n_particles);
    std::vector<arma::vec> k2_r(n_particles);
    std::vector<arma::vec> k2_v(n_particles);
    std::vector<arma::vec> k3_r(n_particles);
    std::vector<arma::vec> k3_v(n_particles);
    std::vector<arma::vec> k4_r(n_particles);
    std::vector<arma::vec> k4_v(n_particles);

    
    //Using forward Euler to move all particles one time-step
    //evolve_forward_Euler(dt); // *** I added this, this is similar to what i think we should do

    for (int i=0; i<particles.size(); i++)
    {   
        arma::vec v = particles.at(i).v;
        arma::vec F = total_force(i);
        
        k1_v.at(i) = (F/m)*dt;
        k1_r.at(i) = v*dt;
        
    }

    for (int i=0; i<particles.size(); i++)
    {   
        arma::vec v = particles_copy.at(i).v;
        arma::vec r = particles_copy.at(i).r;

        particles.at(i).v = v + 0.5*k1_v.at(i);  //*** deleted dt
        particles.at(i).r = r + 0.5*k1_r.at(i); //*** deleted dt
        
    }

    for (int i=0; i<particles.size(); i++)
    {   
        arma::vec F = total_force(i);
        arma::vec v = particles.at(i).v;

        k2_v.at(i) = dt * F/m;
        k2_r.at(i) = dt * v;
       
    }

    for (int i=0; i<particles.size(); i++)
    {   
        arma::vec v = particles_copy.at(i).v;
        arma::vec r = particles_copy.at(i).r; 
        particles.at(i).v = v + 0.5*k2_v.at(i);  //*** deleted dt
        particles.at(i).r = r + 0.5*k2_r.at(i); //*** deleted dt
        
    }

    for (int i=0; i<particles.size(); i++)
    {   
        arma::vec F = total_force(i);
        arma::vec v = particles.at(i).v;
        k3_v.at(i) = dt * F/m;
        k3_r.at(i) = dt * v;
        
    }

    for (int i=0; i<particles.size(); i++)
    {   
        arma::vec v = particles_copy.at(i).v;
        arma::vec r = particles_copy.at(i).r;
        particles.at(i).v = v + k3_v.at(i);  //*** deleted dt
        particles.at(i).r = r + k3_r.at(i); //*** deleted dt
        
    }

    for (int i=0; i<particles.size(); i++)
    {   
        arma::vec F = total_force(i);
        arma::vec v = particles.at(i).v;  
        k4_v.at(i) = dt * F/m;
        k4_r.at(i) = dt * v; 
    }
    for (int i=0; i<particles.size(); i++)
    { 
        arma::vec v = particles_copy.at(i).v;
        arma::vec r = particles_copy.at(i).r;

        particles.at(i).v = v + 1/6.*(k1_v.at(i) + 2.*k2_v.at(i) + 2.*k3_v.at(i) + k4_v.at(i) );
        particles.at(i).r = r + 1/6.*(k1_r.at(i) + 2.*k2_r.at(i) + 2.*k3_r.at(i) + k4_r.at(i) );

        //if norm of r for particle i is smaller than d, set v to zero
        if (arma::norm(particles.at(i).r) > d && particles.at(i).q!=0)
        {
            particles.at(i).v = arma::vec(std::vector<double> {0,0,0});
            particles.at(i).q=0;

        }
        
        
    }
    

}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
    std::vector<Particle> particles_new = particles;

    for (int i = 0; i < particles.size(); i++)
    {
        double q=particles.at(i).q;
        double m=particles.at(i).m;
        

        arma::vec r = particles.at(i).r;
        arma::vec v = particles.at(i).v;
        arma::vec F = total_force(i);
        arma::vec a = F/m;


        //Evolve using foward euler
        particles_new.at(i).r=r+v*dt;
        particles_new.at(i).v= v + a*dt;
    }

    // update all global particles *** updates how?
    // ***  This is what Even added
    for (int i=0; i < particles.size(); i++)
    {
    particles.at(i).r = particles_new.at(i).r;
    particles.at(i).v = particles_new.at(i).v;
    }
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j)
{

    //define empty 3 dimentional vector for our force
    arma::vec force_ij = arma::vec(3, arma::fill::zeros);

    //define some constants and other variables that will be used 
    double k_e = 1.38935333*pow(10.,5.);
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
    //If the perturbation is on we use the perturbed field
    if (perturbation == true)
    {
    force_external = q*external_E_field_perturbed(r, time,omega,f) + q*arma::cross(v,external_B_field(r));
    std::cout << "Perturbation is on" << std::endl;
    std::cout << "The perturbed E-field is: " << external_E_field_perturbed(r, time,omega,f) << std::endl;
    std::cout << "The unperturbed E-field is: " << external_E_field(r) << std::endl;
    
    std::cout << "At time: " << time << std::endl;
    
    }
    
    //If the perturbation is off we use the unperturbed field
    else
    
    {
    force_external = q*external_E_field(r) + q*arma::cross(v,external_B_field(r) );
    }

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
    if (n_particles>1){
    //calculate the force on particle_i from all the other particles
    for (int j = 0; j < n_particles; j++){
        if (j != i){
            force_particles += force_particle(i,j);
        }
    }
    }
    return force_particles;

}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i)
{   
    
    arma::vec force_total = arma::vec(3, arma::fill::zeros);

    if (interaction==false){
        force_total = total_force_external(i);

    }
    else if (interaction==true){

        force_total = total_force_external(i) + total_force_particles(i);

    }
    else{
        std::cout << "Error: interaction must be 0 or 1" << std::endl;
    }
    
    return force_total;

}


