#ifndef __Schrodinger_hpp__
#define __Schrodinger_hpp__
#include <armadillo>
#define ARMA_USE_SUPERLU
#include <complex>
using namespace std;

class Schrodinger
{
public:
    double v0; // potential
    double h;  // step size spacial
    double dt; // step size time
    int M;     // number of steps

    arma::cx_mat V;    // potential matrix
    arma::sp_cx_mat A; // A matrix
    arma::sp_cx_mat B; // B matrix
    arma::superlu_opts opts;

    // constructor
    Schrodinger(double v0_in, double h_in, double dt_in, int M_in);

    // Mapping 2D indices to 1D
    int ij_k(int i, int j);

    // position to index
    int pos_to_idx(double pos);

    // set potential
    void set_potential(string slits);

    // set A and B
    void set_A_B();

    // Solve for next timestep
    arma::cx_vec evolve(arma::cx_vec u);

    // Gaussian wavepacket function
    arma::cx_double gaussian(double x, double y, double sigma_x, double sigma_y, double x0, double y0, double px, double py);

    // Initialise state
    arma::cx_mat initialise_state(double sigma_x, double sigma_y, double x0, double y0, double px, double py);

    // convert matrix to vector
    arma::cx_vec get_u_vec(arma::cx_mat U);
};

#endif