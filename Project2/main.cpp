#include <string>
#include "include/utils.hpp"

int main(int argc, char *argv[])
{
    // Which problem to solve?
    std::string problem = argv[1];

    if (problem == "P2")
    {
        std::cout << ">>> PROBLEM 2 COMPARISSON OF SOLUTIONS" << std ::endl;
        ;
        int N = 6; // dim of square matrix

        arma::mat A = tridiagonal(N);
        arma::vec eigval; // vector of eigenvalues
        arma::mat eigvec; // matrix of eigenvectors

        // Armadillo's solution for the eigen-problem is stored in eigval and eigvec
        eig_sym(eigval, eigvec, A);

        // Normalize the eigenvector as it will be useful in the future
        arma::mat norm_eigvec = normalise(eigvec); // normalized matrix of eigenvectors

        // We can check all eigenvectors and eigenvalues at once
        arma::mat mat_of_eivals = arma::diagmat(eigval);

        // Left side of the eigenvalue problem
        arma::mat left_side = A * norm_eigvec;

        // Right side of the eigenvalue problem
        arma::mat rigt_side = norm_eigvec * mat_of_eivals;

        // Subtract both sides to compare if we are solving the eigen-problem
        std::cout << left_side - rigt_side << '\n';
    }
    else if (problem == "P3")
    {
        // Below is the test matrix
        int N = 4; // dim of square matrix
        arma::mat A = arma::mat(N, N).fill(0.);
        // sets up the test matrix
        A.diag() = arma::vec(N).fill(1);
        A(1, 2) = -0.7;
        A(2, 1) = -0.7;
        A(0, 3) = 0.5;
        A(3, 0) = 0.5;

        int k = 0;
        int l = 1;

        double max = max_offdiag_symmetric(A, k, l, N);

        std::cout << ">>> PROBLEM 3\n"
                  << ">>> The largest (abs val) offdiag element in \n"
                  << A << " is: \n"
                  << ">>> " << max << std ::endl;
    }
    else if (problem == "P4")
    {
        // Below is the test matrix
        int N = 6; // dim of square matrix
        arma::mat A = tridiagonal(N);

        // testing the method
        int k = 0;
        int l = 1;
        int iterations;

        arma::mat diag_A = jacobi(A, k, l, N, "val", iterations);
        arma::mat R_vec = jacobi(A, k, l, N, "vec", iterations);
        arma::vec evals_jacobi_vec = arma::vec(N).fill(1.);

        arma::vec evals_vec = arma::vec(N);
        arma::mat evecs_mat = arma::mat(N, N).fill(1.);
        arma::mat evecs_mat_norm = arma::mat(N, N).fill(1.);
        evals_vec = check_analytical_eval(N, A);
        evecs_mat = check_analytical_evec(N, evecs_mat);
        evecs_mat_norm = normalise(evecs_mat);

        for (int i = 0; i < N; i++)
        {
            evals_jacobi_vec(i) = diag_A(i, i);
        }

        arma::uvec eval_index = arma::sort_index(evals_jacobi_vec);
        arma::vec evals_jacobi_vec_sorted = evals_jacobi_vec(eval_index);

        arma::mat evecs_jacobi_sorted = R_vec.cols(eval_index);
        std::cout << ">>> PROBLEM 4\n"
                  << ">>> evals_vec_analytical: \n " << evals_vec << " | \n"
                  << ">>> evals_jacobi_vec_sorted: \n " << evals_jacobi_vec_sorted << " | \n"
                  << ">>> evectors analytical normalised : \n " << evecs_mat_norm << " | \n"
                  << ">>> evec_jacobi_sorted: \n " << evecs_jacobi_sorted << " | " << std::endl;

        // Here, we expect to output to be zero as we are subtracting the analytical matrix containing our eigenvectors from our jacobi
        // our jacobi matrix has random negative vectors which are still consistent with our results when we find the absolute difference below
        // this is a result of
        std::cout << "checking if outputs match: \n " << arma::abs(evecs_mat_norm) - arma::abs(evecs_jacobi_sorted) << "|" << std::endl;
    }
    else if (problem == "P5")
    {
        std::cout << ">>> PROBLEM 5\n";
        std::vector<std::string>
            A_profiles_vec{"Tridiag", "Dense"};

        // Format parameters
        int width = 10;
        int prec = 5;
        // HERE WE CAN START A LOOP WITH VARIOUS N VALUES and dense or not
        for (std::string A_profile : A_profiles_vec)
        {
            std::string filename = "NOperationsP5" + A_profile + ".txt";
            // Create output file stream
            std::ofstream ofile;
            ofile.open(filename);

            for (int N = 2; N <= 100; N++)
            {
                int k = 0, l = 1;
                // for counting the iterations
                int iterations;

                if (A_profile != "Dense")
                {
                    // sets up the tridiagonal
                    arma::mat A = tridiagonal(N);
                    // Diagonalize
                    arma::mat diag_A = jacobi(A, k, l, N, "val", iterations);
                    output_iterations_to_file(width, prec, ofile, N, iterations);
                }
                else
                {
                    // Generate random N*N matrix
                    arma::mat A = arma::mat(N, N).randn();
                    // Symmetrize the matrix by reflecting the upper triangle to lower triangle
                    A = arma::symmatu(A);
                    arma::mat diag_A = jacobi(A, k, l, N, "val", iterations);
                    output_iterations_to_file(width, prec, ofile, N, iterations);
                }
            }
            ofile.close();
        }
    }
    else if (problem == "P6")
    {
        int N = atoi(argv[2]); // dim of sqr matrix

        arma::mat A = tridiagonal(N);

        int k = 0;
        int l = 1;
        int iterations;

        // Diagonalize and get the rotation matrix
        arma::mat diag_A = jacobi(A, k, l, N, "val", iterations, std::pow(10, -3));
        arma::mat R_vec = jacobi(A, k, l, N, "vec", iterations, std::pow(10, -3));

        // Will be the vetor of eigenvalues
        arma::vec evals_jacobi_vec = arma::vec(N).fill(1.);
        arma::vec evals_vec = arma::vec(N);

        // Will be the matrix of eigenvectors
        arma::mat evecs_mat = arma::mat(N, N).fill(1.);
        arma::mat evecs_mat_norm = arma::mat(N, N).fill(1.);

        evals_vec = check_analytical_eval(N, A);
        evecs_mat = check_analytical_evec(N, evecs_mat);

        evecs_mat_norm = normalise(evecs_mat);

        for (int i = 0; i < N; i++)
        {
            evals_jacobi_vec(i) = diag_A(i, i);
        }

        arma::uvec eval_index = arma::sort_index(evals_jacobi_vec);
        arma::vec evals_jacobi_vec_sorted = evals_jacobi_vec(eval_index);

        arma::mat evecs_jacobi_sorted = R_vec.cols(eval_index);

        std::cout << "evals_jacobi_vec_sorted: \n " << evals_jacobi_vec_sorted << " | " << std::endl;
        std::cout << "evec_jacobi_sorted: \n " << evecs_jacobi_sorted << " | " << std::endl;

        // Eigenboys refer to eigenvectors found using Jacobi's rotation method
        std::string filename = "P6eigenboys.txt";

        // Create output file stream
        std::ofstream ofile;
        ofile.open(filename);
        int width = 30;
        int prec = 10;
        arma::vec X = arma::linspace(0, 1, N);

        // output the eigenboys
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                ofile << std::setw(width) << std::setprecision(prec) << evecs_jacobi_sorted(i, j);
            }
            ofile << std::endl;
        }

        ofile.close();
    }
    return 0;
}