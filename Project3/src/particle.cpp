#include "../Include/Particle.hpp"

Particle::Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in)
{
    q = q_in;
    m = m_in;
    r = r_in;
    v = v_in;
}

void Particle::print_attributes()
{
    std::cout << "Particle's charge is: " << q << std::endl;
    std::cout << "Particle's mass   is: " << m << std::endl;
}
