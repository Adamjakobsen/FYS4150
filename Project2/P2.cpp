#include <iostream>
#include <armadillo>

int main()
{
    int N = 6;     // dim of square matrix
    int n = N + 1; // steps
    double h = 1. / n;
    // initialising v that I will be solving for
    // giving it a dummy first index
    // it lacks the last boundary condition index
    // fill in vector with 1's just to know what's up
    arma::vec v = arma::vec(N).fill(1.);
    v(0) = 0.;
    // std::cout << v_star;
    double a = -1 / (h * h);
    double d = 2 / (h * h);

    arma::mat A = arma::mat(N, N).fill(0.);
    // sets up the tridiagonal
    A.diag() = arma::vec(N).fill(d);
    A.diag(1) = arma::vec(N - 1).fill(a);
    A.diag(-1) = arma::vec(N - 1).fill(a);

    arma::vec eigval; // vector of eigenvalues
    arma::mat eigvec; // matrix of eigenvectors
    eig_sym(eigval, eigvec, A);

    // not important yet, only needed for presentation at the end
    arma::mat norm_eigvec = normalise(eigvec); // normalized matrix of eigenvectors

    // realize we can check all at once
    for (int i = 0; i < N; i++)
    {
        arma::mat left_side = A * eigvec.col(i);
        arma::mat rigt_side = eigval(i) * eigvec.col(i);
        std::cout << left_side - rigt_side << '\n';
    }
}
