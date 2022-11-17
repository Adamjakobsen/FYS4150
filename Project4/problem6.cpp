#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <random> // for std::mt19937
#include <omp.h>
#include <armadillo>
#include <vector>
#include <algorithm> // for openmp ?
#include <stdlib.h>  /* srand, rand */

arma::mat init_random_config(int L, double &exp_E, double &exp_M, int align);
arma::mat find_interacting_pairs(int index, arma::vec raveled_config, int L);
arma::mat evolve(arma::mat config, double beta, double &exp_E, double &exp_M, arma::vec p_vec);
int mt_random_int(int low, int high);
double mt_random_float(int low, int high);

void monte_carlo(int L, int mc_cycles, int burn_pct, double T, int align);
std::random_device rd;
std::mt19937 gen(rd()); // gen(2) for a seed of 2
// define your physical system's variables here:

double kb = 1.;
int main(int argc, char *argv[])
{
    // main will only evolve the system and output the results to a file in matrix form
    // we will leave the find the total energy of config and other quantities of the system to python
    // here we just evolve and output the results

    int L = atoi(argv[1]);

    int mc_cycles = atoi(argv[2]); // 1 = L^2 runs, 2 = 2*L^2 runs, etc
    int burn_pct = atoi(argv[3]);  // 10 = 10% burn in, 20 = 20% burn in, etc
    double lower_temp = atof(argv[4]);
    double upper_temp = atof(argv[5]);
    double temp_step = atof(argv[6]);
    int align = atoi(argv[7]); // 0 = random, 1 = aligned up (1), 2 = aligned down (-1)
    double original_temp_step = temp_step;

#pragma omp parallel for
    for (double T = lower_temp; T <= upper_temp; T += temp_step)
    {
        if (T > 2.2 && T <= 2.4)
        {
            temp_step = 0.1 * original_temp_step;
        }
        if (T > 2.4)
        {
            temp_step = original_temp_step;
        }
        monte_carlo(L, mc_cycles, burn_pct, T, align);
    }
}

arma::mat init_random_config(int L, double &E, double &M, int align)
{
    int N = L * L;
    arma::mat config = arma::zeros<arma::mat>(L, L);
    arma::mat raveled_config = arma::reshape(config, N, 1);
    arma::mat interacting_pairs = arma::zeros<arma::mat>(4, 2);

    for (int i = 0; i < N; i++)
    {
        if (align == 0)
        {
            raveled_config(i) = mt_random_int(0, 1) * 2 - 1; // makes -1 or 1 with equal probability
        }
        else if (align == 1)
        {
            raveled_config(i) = 1;
        }
        else if (align == 2)
        {
            raveled_config(i) = -1;
        }

        M += raveled_config(i);
    }
    // notice that to find the initial energy, we need another loop
    for (int i = 0; i < N; i++)
    {
        arma::mat interacting_pairs = find_interacting_pairs(i, raveled_config, L);
        for (int j = 0; j < 4; j++) // always has 4 neghbors so this is general
        {
            E += -(interacting_pairs(j, 0) * interacting_pairs(j, 1));
        }
    }
    E = E / 2; // since we double count each interaction

    // std::cout << "Initial energy: " << E << std::endl;
    return arma::reshape(raveled_config, L, L);
}

arma::mat find_interacting_pairs(int index, arma::vec raveled_config, int L)
{
    // create a vector of vectors to store the interacting pairs

    arma::mat interacting_pairs = arma::zeros<arma::mat>(4, 2); // rows are the interacting pairs, columns are the coordinates of the interacting pairs
    arma::mat side_neighbours = arma::zeros<arma::vec>(2);
    arma::mat vert_neighbours = arma::zeros<arma::vec>(2);

    // notice the intentional use of integer division
    side_neighbours(0) = (L + index - 1) % L + L * (index / L); // same as index-1, but if index is 0, it will wrap around to the end of the row
    side_neighbours(1) = (L + index + 1) % L + L * (index / L); // same as index+1, but if index is L-1, it will wrap around to the start of the row

    vert_neighbours(0) = (L * L + index - L) % (L * L); // same as index-L, but if index is 0, it will wrap around to the end of the column
    vert_neighbours(1) = (L * L + index + L) % (L * L); // same as index+L, but if index is L*L-1, it will wrap around to the start of the column

    int index_value = raveled_config(index);

    interacting_pairs.col(0) = arma::vec(4).fill(index_value);
    //  there must be a better way to do this
    interacting_pairs(0, 1) = raveled_config(side_neighbours(0));
    interacting_pairs(1, 1) = raveled_config(side_neighbours(1));
    interacting_pairs(2, 1) = raveled_config(vert_neighbours(0));
    interacting_pairs(3, 1) = raveled_config(vert_neighbours(1));

    return interacting_pairs;
}

arma::mat evolve(arma::mat config, double beta, double &E, double &M, arma::vec p_vec)
{
    int L = config.n_rows;
    arma::mat raveled_config = arma::reshape(config, L * L, 1);
    arma::mat interacting_pairs = arma::zeros<arma::mat>(4, 2);
    for (int iter = 0; iter < L * L; iter++)
    {
        double neigbours_energy = 0; // energy of the neighbours of the randomly chosen spin that is why it has to be reset
        double delta_E = 0;

        int i = mt_random_int(0, L * L - 1); // choose a random spin index to flip
        // std::cout << i << std::endl;
        interacting_pairs = find_interacting_pairs(i, raveled_config, L);
        for (int j = 0; j < 4; j++) // always has 4 neghbors so this is general
        {
            neigbours_energy += -(interacting_pairs(j, 0) * interacting_pairs(j, 1));
        }

        delta_E = -2 * neigbours_energy; // because because delta_E = E_new - E_old = -((-a)*b) - (-(a*b)) = 2*(a*b) = -2*neigbours_energy

        if (delta_E <= 0)
        {
            raveled_config(i) *= -1;
            E += delta_E;
            M += 2 * raveled_config(i); // factor of 2 because we are only changing one spin
        }
        else
        {
            double p = std::exp(-beta * delta_E);
            // double p = p_vec[(delta_E + 8) / 4]; // this did not make it faster

            if (mt_random_float(0, 1) < p)
            {
                raveled_config(i) *= -1;
                E += delta_E;
                M += 2 * raveled_config(i); // factor of 2 because we are only changing one spin
            }
        }
    }
    return arma::reshape(raveled_config, L, L);
}

void monte_carlo(int L, int mc_cycles, int burn_pct, double T, int align)
{
    double beta = 1. / (kb * T);
    arma::vec p_vec = arma::zeros<arma::vec>(5);
    arma::vec delta_e_vec = arma::vec(std::vector<double>{-8, -4, 0, 4, 8});
    for (int i = 0; i < 5; i++)
    {
        p_vec(i) = std::exp(-beta * delta_e_vec(i));
        // std::cout << "delta_e_vec(i)" << delta_e_vec(i) << std::endl;
    }

    std::ofstream file2;

    std::string rounded_T = std::to_string(T).substr(0, std::to_string(T).find(".") + 3 + 1);

    std::cout
        << "T: " << T << std::endl;
    std::string filename2 = "problem6/" + std::to_string(L) + "/epsilon_L" + std::to_string(L) + "_A" + std::to_string(align) + "_mc" + std::to_string(mc_cycles) + "_burn" + std::to_string(burn_pct) + "_t" + rounded_T + ".txt";

    file2.open(filename2);

    // avg will always be wrt to the number of mc cycles
    double E = 0;
    double M = 0;

    double cumul_E = 0;
    double cumul_M = 0;

    double cumul_e = 0;
    double cumul_m = 0;

    double cumul_E2 = 0;
    double avg_e = 0;
    double avg_E = 0;
    double avg_E2 = 0;

    double avg_e2 = 0;
    double var_E = 0;

    double var_M = 0;
    double avg_M = 0;
    double avg_M2 = 0;

    double avg_mabs = 0;

    double cumul_M2 = 0;
    double cumul_m2 = 0;
    double cumul_e2 = 0;
    double avg_Mabs = 0;

    double avg_m2 = 0;

    double avg_Mabs2 = 0;
    double cumul_Mabs = 0;

    double cumul_mabs = 0;

    double Cv = 0;
    double chi = 0;
    arma::mat config = init_random_config(L, E, M, align);

    double N = L * L;
    int time_steps = mc_cycles * N;
    int burn_in = int((burn_pct / 100.) * time_steps);
    // std::cout << "burned mc_cycles: " << burn_in << std::endl;

    // std::cout << "time_steps = " << time_steps << std::endl;
    for (int i = 0; i < time_steps; i++)
    {
        config = evolve(config, beta, E, M, p_vec);
        // the following qtd will be calculated in a dumb way just to be didatic but uses lots of flops

        cumul_E += E;
        std::cout << "E = " << E << std::endl;
        //  if (i > burn_in)
        //{
        // if (i % int(N) == 0 && i != 0)
        //{

        file2 << std::setw(25) << std::setprecision(15)
              << E / N << std::endl;
        //}
        //}
    }
    file2.close();
}

int mt_random_int(int low, int high)
{
    std::uniform_int_distribution<> dist(low, high);
    return dist(gen);
}

double mt_random_float(int low, int high)
{
    std::uniform_real_distribution<double> dist2(low, high);
    return dist2(gen);
}