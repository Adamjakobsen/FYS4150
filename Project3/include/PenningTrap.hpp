#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <armadillo>
#include <iostream>
#include "Particle.hpp"


class PenningTrap
{
    public:
        // these are attributes aka variables (or constants) of a class
        double B0, V0, d;
        bool interaction;
        bool perturbation;
        std::vector<Particle> particles;

        PenningTrap(double B0_in, double V0_in, double d_in);

        // Add a particle to the trap   
        void add_particle(Particle p_in);

        //method for evaluating ext e-field 
        arma::vec external_E_field(arma::vec r);

        //method for evaluating ext b-field
        arma::vec external_B_field(arma::vec r);  
        
        // Evolve the system one time step (dt) using Runge-Kutta 4th order
        void evolve_RK4(double dt);

        // Evolve the system one time step (dt) using Forward Euler
        void evolve_forward_Euler(double dt);

        // Force on particle_i from particle_j
        arma::vec force_particle(int i, int j);

        // The total force on particle_i from the external fields
        arma::vec total_force_external(int i);

        // The total force on particle_i from the other particles
        arma::vec total_force_particles(int i);

        // The total force on particle_i from both external fields and other particles
        arma::vec total_force(int i);

        // Option to turn on/off the interaction between particles
        void toggle_interaction(bool interaction);

        // External E_field with perturbation
        arma::vec external_E_field_perturbed(arma::vec r, double t, double omega,double f);

        // Option to turn on/off the perturbation
        void toggle_perturbation(bool perturbation);




};

#endif