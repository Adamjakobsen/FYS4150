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
        std::vector<Particle> particles;

        PenningTrap(double B0_in, double V0_in, double d_in);
        //method for evaluating ext e-field 
        arma::vec external_E_field(arma::vec r);
        //method for evaluating ext b-field
        arma::vec external_B_field(arma::vec r);  
        
        void add_particle(Particle p_in);
        void evolve_RK4(double dt);
        void evolve_forward_Euler(double dt);

};

#endif