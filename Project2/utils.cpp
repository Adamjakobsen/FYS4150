#include "include/utils.hpp"

arma::mat tridiagonal(int N)
{
    // returns the tridiagonal matrix of the required dimensions
    int n = N + 1; // steps

    double h = 1. / n;       // stepsize
    double a = -1 / (h * h); // diagonal elements
    double d = 2 / (h * h);  // sub and super-diagonal elements

    arma::mat A = arma::mat(N, N).fill(0.);

    // sets up the proper tridiagonal
    A.diag() = arma::vec(N).fill(d);
    A.diag(1) = arma::vec(N - 1).fill(a);
    A.diag(-1) = arma::vec(N - 1).fill(a);
    return A;
}

double max_offdiag_symmetric(arma::mat A, int &k, int &l, int N)
{
    // Will return the largest offdiag element of A.
    // Also updates k and l value in memory to be used later.
    k = 0;
    l = 1;
    for (int i = 0; i <= N - 1; i++)
    {
        for (int j = i + 1; j <= N - 1; j++)
        {
            if (abs(A(i, j)) >= abs(A(k, l)))
            {
                k = i;
                l = j;
            }
        }
    }
    return A(k, l);
}

arma::vec check_analytical_eval(int N, arma::mat A)
{
    double a = A(0, 0);
    double d = A(0, 1);
    arma::vec evals_vec = arma::vec(N).fill(1.);
    for (int i = 0; i < N; i++)
    {
        evals_vec(i) = d + 2 * a * cos((i + 1) * PI / (N + 1));
    }
    return evals_vec;
}

arma::mat check_analytical_evec(int N, arma::mat evecs_mat)
{
    arma::vec evec_i;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            evecs_mat(i, j) = sin((i + 1.) * (j + 1.) * PI / (N + 1.));
        }
    }
    return evecs_mat;
}

arma::mat jacobi(arma::mat A, int k, int l, int N, std::string ask, int &iterations, double epsilon)
{
    iterations = 0;
    // step 1 - initialization
    // we have defined epsilon and passed A. Now we need R
    arma::mat R = arma::eye(N, N);

    // step 2 - call the max_offdiag_symmetric
    double max = max_offdiag_symmetric(A, k, l, N);

    // step 3 while loop
    double tau, t, c, s;
    double a_k_k, a_i_k, r_i_k; // temporary variables

    while (std::abs(max) > epsilon)
    {
        iterations += 1;
        // 3.1
        tau = (A(l, l) - A(k, k)) / (2 * max);

        // 3.2
        t = std::min(1 / (tau + sqrt(1 + tau * tau)), -1 / (-tau + sqrt(1 + tau * tau)));
        c = 1 / sqrt(1 + t * t);
        s = c * t;

        // 3.3
        a_k_k = A(k, k); // save before update
        A(k, k) = A(k, k) * c * c - 2 * A(k, l) * c * s + A(l, l) * s * s;
        A(l, l) = A(l, l) * c * c + 2 * A(k, l) * c * s + a_k_k * s * s; ///////////////// EXTRA CAREFUL HERE
        A(k, l) = 0;
        A(l, k) = 0;

        for (int i = 0; i < N; i++)
        {
            // skip the k,l elements
            if (i == k || i == l)
            {
                continue;
            }

            a_i_k = A(i, k); // save before update
            A(i, k) = A(k, i) = A(i, k) * c - A(i, l) * s;
            // a(k,i) = a(i,k);
            A(i, l) = A(l, i) = A(i, l) * c + a_i_k * s; // EXTRA CAREFUL HERE
            // a(l,i) = a(i,l)
        }

        // 3.4 update rotation matrix
        for (int i = 0; i < N; i++)
        {
            r_i_k = R(i, k); // save before update
            R(i, k) = R(i, k) * c - R(i, l) * s;
            R(i, l) = R(i, l) * c + r_i_k * s; // EXTRA CAREFUL HERE
        }

        // 3.5 find the big chungus again
        max = max_offdiag_symmetric(A, k, l, N);
    }

    if (ask == "val")
    {
        std::cout << "N: " << N << ",  "
                  << "iterations: " << iterations << "\n"
                  << std::endl;
        return A; // this is after the loop so this matrix is ready
    }
    else
    {
        return R;
    }
}

void output_iterations_to_file(int width, int prec, std::ofstream &ofile, int N, int iterations)
{
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << N
          << std::setw(width) << std::setprecision(prec) << std::scientific << iterations << std::endl;
}