#include <iostream>
#include <cmath>

#include <armadillo>
#include "Particle.hpp"
#include "PenningTrap.hpp"

int main()
{
    arma::vec v = {0.,0.,0.};
    arma::vec r = {20.,0.,20.};

    //remember to define correct elementary charge value here!!
    Particle particle = Particle(1.6, 2.5, r, v);
    PenningTrap PT = PenningTrap(1., 25., 500.);

    particle.print_attributes();
    arma::vec electric_field = PT.external_E_field(particle.r);

    std::cout << 'External E-field is' << electric_field << std::endl;
    return 0;
}
