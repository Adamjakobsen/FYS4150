#include <iostream>
#include <cmath>

#include <armadillo>

const double PI = 4. * atan(1.);

arma::vec check_analytical_eval(int N, double a, double d)
{
    arma::vec evals_vec = arma::vec(N).fill(1.);

    for (int i = 0; i < N; i++)
    {
        evals_vec(i) = d + 2 * a * cos((i + 1) * PI / (N + 1));
    }

    return evals_vec;
}

arma::mat check_analytical_evec(int N, arma::mat evecs_mat, double a)
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

double max_offdiag_symmetric(arma::mat &A, int N, int &k, int &l)
{
    // we will call this funtion multiple times so it is important to always reinitialize since we passed the memory address
    k = 0;
    l = 1;
    // this is O(n^2) unfortunately because we need to return the indexes.
    for (int i = 0; i <= N - 1; i++)
    {
        for (int j = i + 1; j <= N - 1; j++)
        {
            if (std::abs(A(i, j)) >= std::abs(A(k, l)))
            {
                k = i;
                l = j;
            }
        }
    }

    return A(k, l);
}

arma::mat jacobi(arma::mat A, int k, int l, int N, std::string ask, double epsilon = std::pow(10, -8))
{
    // step 1 - initialization
    // we have defined epsilon and passed A. Now we need R
    arma::mat R = arma::eye(N, N);

    // step 2 - call the max_offdiag_symmetric
    double max = max_offdiag_symmetric(A, N, k, l);

    // step 3 while loop
    double tau, t, c, s;
    double a_k_k, a_i_k, r_i_k; // temporary variables

    while (std::abs(max) > epsilon)
    {
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
            A(i, l) = A(l, i) = A(i, l) * c + a_i_k * s; ///////////////// EXTRA CAREFUL HERE
            // a(l,i) = a(i,l)
        }

        // 3.4 update rotation matrix
        for (int i = 0; i < N; i++)
        {
            r_i_k = R(i, k); // save before update
            R(i, k) = R(i, k) * c - R(i, l) * s;
            R(i, l) = R(i, l) * c + r_i_k * s; ///////////////// EXTRA CAREFUL HERE
        }

        // 3.5 find the big chungus again
        max = max_offdiag_symmetric(A, N, k, l);
    }

    if (ask == "val")
    {
        return A; // this is after the loop so this matrix is ready to eat
    }
    else
    {
        return R;
    }
}

int main()
{
    // Below is the test matrix
    int N = 6;     // dim of sqr matrix
    int n = N + 1; // steps
    double h = 1. / n;

    double a = -1 / (h * h);
    double d = 2 / (h * h);

    arma::mat A = arma::mat(N, N).fill(0.);
    // sets up the tridiagonal
    A.diag() = arma::vec(N).fill(d);
    A.diag(1) = arma::vec(N - 1).fill(a);
    A.diag(-1) = arma::vec(N - 1).fill(a);

    // testing the method
    int k = 0;
    int l = 1;

    arma::mat diag_A = jacobi(A, k, l, N, "val");
    arma::mat R_vec = jacobi(A, k, l, N, "vec");
    std::cout << "evectors_jacobi: \n " << R_vec << " | " << std::endl;

    arma::vec evals_jacobi_vec = arma::vec(N).fill(1.);

    arma::vec evals_vec = arma::vec(N);
    arma::mat evecs_mat = arma::mat(N, N).fill(1.);
    evals_vec = check_analytical_eval(N, a, d);
    evecs_mat = check_analytical_evec(N, evecs_mat, a);

    for (int i = 0; i < N; i++)
    {
        evals_jacobi_vec(i) = diag_A(i, i);
    }

    std::cout << "evals_vec analytical: \n " << evals_vec << " | " << std::endl;
    std::cout << "evals_jacobi_vec: \n " << evals_jacobi_vec << " | " << std::endl;

    // ALL THE EIGENVALUES ARE THERE, THEY ARE JUST ARANGED IN A DIFERENT ORDER
    // IF WE WERE TO REARANGE THEM, WE WOULD DO
    // 1 -> 1
    // 2 -> 3
    // 3 -> 5
    // 4 -> 6
    // 5 -> 4
    // 6 -> 2
    // SO THIS IS WHAT HAPPENED TO THE COLUMNS OF R AS WELL, WITH RESPECT TO THE COLUMNS OF THE MATRIX CONTAINING THE EIGENVECTORS

    std::cout << "evectors analytical : \n " << evecs_mat << " | " << std::endl;
    std::cout << "evectors_jacobi: \n " << R_vec << " | " << std::endl;

    return 0;
}