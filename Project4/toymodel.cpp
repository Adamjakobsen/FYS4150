#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <vector>

arma::mat init_random_config(int L);
arma::mat find_interacting_pairs(int index, arma::vec raveled_config, int L);
arma::mat evolve(arma::mat config, double beta);

// define your physical system's variables here:

double kb = 1.;
double T = 1.;
double beta = 1. / (kb * T);

int main(int argc, char *argv[])
{
    // main will only evolve the system and output the results to a file in matrix form
    // we will leave the find the total energy of config and other quantities of the system to python
    // here we just evolve and output the results
    int L = atoi(argv[1]);
    int epochs = atoi(argv[2]);
    arma::mat config = init_random_config(L);
    std::ofstream file;
    file.open("config.txt");
    for (int i = 0; i < epochs; i++)
    {
        config = evolve(config, beta);
        file << config << std::endl;
    }
}

arma::mat init_random_config(int L)
{
    arma::mat config = arma::zeros<arma::mat>(L, L);
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            config(i, j) = 2 * (rand() % 2) - 1; // the %2 makes it be a number between 0 and 1, and the 2* makes it be a number between 0 and 2, and the -1 makes it be a number between -1 and 1
        }
    }
    return config;
}

arma::mat find_interacting_pairs(int index, arma::vec raveled_config, int L)
{
    // create a vector of vectors to store the interacting pairs

    arma::mat interacting_pairs = arma::zeros<arma::mat>(4, 2); // rows are the interacting pairs, columns are the coordinates of the interacting pairs
    arma::mat side_neighbours = arma::zeros<arma::vec>(2);
    arma::mat vert_neighbours = arma::zeros<arma::vec>(2);

    // notice the intentional use of integer division
    side_neighbours(0) = (index - 1) % L + L * (index / L); // same as index-1, but if index is 0, it will wrap around to the end of the row
    side_neighbours(1) = (index + 1) % L + L * (index / L); // same as index+1, but if index is L-1, it will wrap around to the start of the row

    vert_neighbours(0) = (index - L) % (L * L); // same as index-L, but if index is 0, it will wrap around to the end of the column
    vert_neighbours(1) = (index + L) % (L * L); // same as index+L, but if index is L*L-1, it will wrap around to the start of the column

    int index_value = raveled_config(index);

    interacting_pairs.col(0) = arma::vec(4).fill(index_value);

    // there must be a better way to do this
    interacting_pairs(0, 1) = raveled_config(side_neighbours(0));
    interacting_pairs(1, 1) = raveled_config(side_neighbours(1));
    interacting_pairs(2, 1) = raveled_config(vert_neighbours(0));
    interacting_pairs(3, 1) = raveled_config(vert_neighbours(1));

    return interacting_pairs;
}

arma::mat evolve(arma::mat config, double beta)
{
    int L = config.n_rows;
    arma::mat raveled_config = arma::reshape(config, L * L, 1);
    arma::mat interacting_pairs = arma::zeros<arma::mat>(4, 2);
    double delta_E = 0;
    double energy = 0;
    for (int i = 0; i < L * L; i++)
    {
        interacting_pairs = find_interacting_pairs(i, raveled_config, L);
        for (int j = 0; j < 4; j++) // always has 4 neghbors so this is general
        {
            energy += -(interacting_pairs(j, 0) * interacting_pairs(j, 1));
        }
        delta_E = -2 * energy; // because because delta_E = E_new - E_old = -((-a)*b) - (-(a*b)) = 2*(a*b) = -2*energy
        if (delta_E <= 0)
        {
            raveled_config(i) *= -1;
        }
        else
        {
            double p = exp(-beta * delta_E);
            if ((rand() % 100) / 100 < p)
            {
                raveled_config(i) *= -1;
            }
        }
    }
    return arma::reshape(raveled_config, L, L);
}