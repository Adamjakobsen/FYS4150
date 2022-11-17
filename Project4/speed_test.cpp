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
arma::mat evolve(arma::mat config, double beta, double &exp_E, double &exp_M, arma::vec p_vec);
int mt_random_int(int low, int high);
double mt_random_float(int low, int high);

void monte_carlo(int L, int mc_cycles, int burn_pct, double T, int align);
std::random_device rd;
std::mt19937 gen(2); // gen(2) for a seed of 2
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
    arma::mat interacting_pairs; //= arma::zeros<arma::mat>(4, 2);

    for (int i = 0; i < N; i++)
    {
        int k = i / L;
        int l = i % L;
        if (align == 0)
        {
            config(k, l) = mt_random_int(0, 1) * 2 - 1; // makes -1 or 1 with equal probability
        }
        else if (align == 1)
        {
            config(k, l) = 1;
        }
        else if (align == 2)
        {
            config(k, l) = -1;
        }

        M += config(k, l);
    }
    // notice that to find the initial energy, we need another loop
    for (int i = 0; i < N; i++)
    {
        int k = i / L;
        int l = i % L;

        E += -config(k, l) * config((k + 1) % L, l);
        E += -config(k, l) * config((k - 1 + L) % L, l);
        E += -config(k, l) * config(k, (l + 1) % L);
        E += -config(k, l) * config(k, (l - 1 + L) % L);
    }
    E = E / 2; // since we double count each interaction

    // std::cout << "Initial energy: " << E << std::endl;
    return config;
}

arma::mat find_interacting_pairs(int k, int l, arma::mat config, int L)
{
    // create a vector of vectors to store the interacting pairs
    arma::mat interacting_pairs = arma::zeros<arma::mat>(4, 2); // rows are the interacting pairs, columns are the coordinates of the interacting pairs

    // notice the intentional use of integer division

    int index_value = config(k, l);
    interacting_pairs.col(0) = arma::vec(4).fill(index_value);
    //  there must be a better way to do this
    interacting_pairs(0, 1) = config(k, (l + 1) % L);
    interacting_pairs(1, 1) = config(k, (l - 1 + L) % L);
    interacting_pairs(2, 1) = config((k + 1) % L, l);
    interacting_pairs(3, 1) = config((k - 1 + L) % L, l);

    return interacting_pairs;
}

arma::mat evolve(arma::mat config, double beta, double &E, double &M, arma::vec p_vec)
{
    int L = config.n_rows;
    arma::mat interacting_pairs = arma::zeros<arma::mat>(4, 2);
    for (int iter = 0; iter < L * L; iter++)
    {
        double neigbours_energy = 0; // energy of the neighbours of the randomly chosen spin that is why it has to be reset
        double delta_E = 0;

        int i = mt_random_int(0, L * L - 1); // choose a random spin index to flip

        // figure out k and l from i
        int k = i / L;
        int l = i % L;

        // std::cout << i << std::endl;

        neigbours_energy += -config(k, l) * config((k + 1) % L, l);
        neigbours_energy += -config(k, l) * config((k - 1 + L) % L, l);
        neigbours_energy += -config(k, l) * config(k, (l + 1) % L);
        neigbours_energy += -config(k, l) * config(k, (l - 1 + L) % L);

        delta_E = -2 * neigbours_energy; // because because delta_E = E_new - E_old = -((-a)*b) - (-(a*b)) = 2*(a*b) = -2*neigbours_energy

        if (delta_E <= 0)
        {
            config(k, l) *= -1;
            E += delta_E;
            M += 2 * config(k, l); // factor of 2 because we are only changing one spin
        }
        else
        {
            // double p = std::exp(-beta * delta_E);
            double p = p_vec[(delta_E + 8) / 4]; // this did not make it faster

            if (mt_random_float(0, 1) < p)
            {
                config(k, l) *= -1;
                E += delta_E;
                M += 2 * config(k, l); // factor of 2 because we are only changing one spin
            }
        }
    }
    return config;
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

    // std::ofstream file1;
    std::ofstream file2;

    std::string rounded_T = std::to_string(T).substr(0, std::to_string(T).find(".") + 3 + 1);

    std::cout
        << "T: " << T << std::endl;
    // std::string filename = std::to_string(L) + "/cfg_L" + std::to_string(L) + "_A" + std::to_string(align) + "_mc" + std::to_string(mc_cycles) + "_burn" + std::to_string(burn_pct) + "_t" + rounded_T + ".txt";
    std::string filename2 = std::to_string(L) + "/_speed_test_fast_qt_L" + std::to_string(L) + "_A" + std::to_string(align) + "_mc" + std::to_string(mc_cycles) + "_burn" + std::to_string(burn_pct) + "_t" + rounded_T + ".txt";

    // file1.open(filename);
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
        // std::cout << "E = " << E << std::endl;

        cumul_e = cumul_E / N;
        cumul_E2 += E * E;
        cumul_e2 = cumul_E2 / (N * N);
        avg_e = cumul_e / (i + 1);
        avg_e2 = cumul_e2 / (i + 1);
        avg_E2 = cumul_E2 / (i + 1);
        avg_E = cumul_E / (i + 1);
        var_E = avg_E2 - avg_E * avg_E;

        cumul_M += M;
        cumul_m = cumul_M / N;
        cumul_M2 += M * M;
        cumul_m2 = cumul_M2 / (N * N);
        cumul_Mabs += std::abs(M);
        avg_Mabs = cumul_Mabs / (i + 1);

        cumul_mabs = cumul_Mabs / N;

        avg_mabs = cumul_mabs / (i + 1);
        avg_m2 = cumul_m2 / (i + 1);
        avg_M = cumul_M / (i + 1);
        avg_M2 = cumul_M2 / (i + 1);
        avg_Mabs2 = avg_Mabs * avg_Mabs;

        var_M = avg_M2 - avg_Mabs2;
        Cv = (var_E) / (N * kb * T * T);
        chi = (var_M) / (N * kb * T);
        if (i > burn_in)
        {
            if (i % int(N) == 0 && i != 0) // output only at the end of a mc_cycle
            {
                // file1
                //     << config << std::endl;
                file2 << std::setw(25) << std::setprecision(15)
                      << avg_e << " " << avg_mabs << " " << Cv << " " << chi
                      << std::endl;
            }
        }
    }
    // file1.close();
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
