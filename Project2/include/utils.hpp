#include <iostream>
#include <armadillo>
#include <cmath>
#include <iomanip>

const double PI = 4. * atan(1.);

arma::mat tridiagonal(int N);
double max_offdiag_symmetric(arma::mat A, int &k, int &l, int N);
arma::vec check_analytical_eval(int N, arma::mat A);
arma::mat check_analytical_evec(int N, arma::mat evecs_mat);
arma::mat jacobi(
    arma::mat A,
    int k,
    int l,
    int N,
    std::string ask,
    int &iterations,
    double epsilon = std::pow(10, -8));
void output_iterations_to_file(
    int width,
    int prec,
    std::ofstream &ofile,
    int N,
    int iterations);
