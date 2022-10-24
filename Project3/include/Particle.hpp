#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

class Particle
{
    public:
        // these are attributes aka variables (or constants) of a class
        double q, m;
        arma::vec r, v;

        double time;

        Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in);
        void print_attributes();
};

#endif
