#include <iostream>
#include <armadillo>

double max_offdiag_symmetric(arma::mat A, int &k, int &l, int N)
{
    // this is O(n^2) unfortunately because we need to return the indexes.
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

int main()
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
    std::cout << max << std ::endl;
}
