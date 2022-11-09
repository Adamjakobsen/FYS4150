#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <vector>
#include <stdlib.h> /* srand, rand */

arma::mat init_random_config(int L, double &exp_E, double &exp_M);
arma::mat find_interacting_pairs(int index, arma::vec raveled_config, int L);
arma::mat evolve(arma::mat config, double beta, double &exp_E, double &exp_M);
void monte_carlo(int L, int mc_cycles);

// define your physical system's variables here:

double kb = 1.;
double T = 1.;
double beta = 1. / (kb * T);

int main(int argc, char *argv[])
{
    // main will only evolve the system and output the results to a file in matrix form
    // we will leave the find the total energy of config and other quantities of the system to python
    // here we just evolve and output the results
    srand(time(NULL));

    int L = atoi(argv[1]);
    int mc_cycles = atoi(argv[2]); // 1 = L^2 runs, 2 = 2*L^2 runs, etc
    monte_carlo(L, mc_cycles);
}

arma::mat init_random_config(int L, double &exp_E, double &exp_M)
{
    int N = L * L;
    arma::mat config = arma::zeros<arma::mat>(L, L);
    arma::mat raveled_config = arma::reshape(config, N, 1);
    arma::mat interacting_pairs = arma::zeros<arma::mat>(4, 2);

    for (int i = 0; i < N; i++)
    {
        raveled_config(i) = 2 * (rand() % 2) - 1; // the %2 makes it be a number between 0 and 1, and the 2* makes it be a number between 0 and 2, and the -1 makes it be a number between -1 and 1
        exp_M += raveled_config(i);
    }
    // notice that to find the initial energy, we need another loop
    for (int i = 0; i < N; i++)
    {
        arma::mat interacting_pairs = find_interacting_pairs(i, raveled_config, L);
        for (int j = 0; j < 4; j++) // always has 4 neghbors so this is general
        {
            exp_E += -(interacting_pairs(j, 0) * interacting_pairs(j, 1));
        }
    }
    return arma::reshape(raveled_config, L, L);
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

arma::mat evolve(arma::mat config, double beta, double &exp_E, double &exp_M)
{
    int L = config.n_rows;
    arma::mat raveled_config = arma::reshape(config, L * L, 1);
    arma::mat interacting_pairs = arma::zeros<arma::mat>(4, 2);
    double delta_E = 0;
    double energy = 0;
    for (int iter = 0; iter < L * L; iter++)
    {
        int i = rand() % (L * L);
        interacting_pairs = find_interacting_pairs(i, raveled_config, L);
        for (int j = 0; j < 4; j++) // always has 4 neghbors so this is general
        {
            energy += -(interacting_pairs(j, 0) * interacting_pairs(j, 1));
        }
        delta_E = -2 * energy; // because because delta_E = E_new - E_old = -((-a)*b) - (-(a*b)) = 2*(a*b) = -2*energy
        if (delta_E <= 0)
        {
            raveled_config(i) *= -1;
            exp_E += delta_E;               // should it not be exp_e += exp(-beta*delta_E) ?
            exp_M += 2 * raveled_config(i); // factor of 2 because we are only changing one spin
        }
        else
        {
            double p = exp(-beta * delta_E);
            if ((rand() % 100) / 100 < p)
            {
                raveled_config(i) *= -1;
                exp_E += delta_E;               // should it not be exp_E += exp(-beta*delta_E) ?
                exp_M += 2 * raveled_config(i); // factor of 2 because we are only changing one spin
            }
        }
    }
    return arma::reshape(raveled_config, L, L);
}

void monte_carlo(int L, int mc_cycles)
{
    std::ofstream file1;
    std::ofstream file2;
    std::string filename = "config_L_" + std::to_string(L) + "_mc_" + std::to_string(mc_cycles) + ".txt";
    std::string filename2 = "quantities_L_" + std::to_string(L) + "_mc_" + std::to_string(mc_cycles) + ".txt";

    file1.open(filename);
    file2.open(filename2);

    // notice these will be normalized per spin in the output file
    double exp_E = 0;
    double exp_E2 = 0;
    double exp_M = 0;
    double exp_M2 = 0;
    double exp_M_abs = 0;
    double Cv = 0;
    double chi = 0;
    arma::mat config = init_random_config(L, exp_E, exp_M);

    double N = L * L;

    int epochs = mc_cycles * N;
    for (int i = 0; i < epochs; i++)
    {
        config = evolve(config, beta, exp_E, exp_M);
        exp_E2 += exp_E * exp_E;
        exp_M2 += exp_M * exp_M;
        exp_M_abs += abs(exp_M);
        Cv = (exp_E2 - exp_E * exp_E) / (N * kb * T * T);
        chi = (exp_M2 - exp_M_abs * exp_M_abs) / (N * kb * T);
        file1 << config << std::endl;
        file2 << exp_E / N << " " << exp_E2 / N << " " << exp_M / N << " " << exp_M2 / N << " " << exp_M_abs / N << " " << Cv << " " << chi << std::endl;
    }
    file1.close();
    file2.close();
}