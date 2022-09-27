#include "include/utils.hpp"

int main(int argc, char *argv[])
{
    int N = atoi(argv[1]); // dim of sqr matrix

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

    return 0;
}
