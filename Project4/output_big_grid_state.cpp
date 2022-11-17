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

void monte_carlo(int L, int mc_cycles, int burn_pct, double T, int align, int thread_id);
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
    int temp_n = int((upper_temp - lower_temp) / temp_step);

#pragma omp parallel
    {
        int i;
        double T;
        double temp_step = atof(argv[6]);
#pragma omp for
        for (i = 0; i < temp_n + 1; i++)
        {
            T = lower_temp + i * temp_step;
            if (T > 2.2 && T <= 2.4)
            {
                temp_step = 0.1 * original_temp_step;
            }
            if (T > 2.4)
            {
                temp_step = original_temp_step;
            }
            monte_carlo(L, mc_cycles, burn_pct, T, align, omp_get_thread_num());
            std::cout << "Thread " << omp_get_thread_num() << " finished T = " << T << std::endl;
        }
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

        E += -config(k, l) * (config((k + 1) % L, l) + config((k - 1 + L) % L, l) + config(k, (l + 1) % L) + config(k, (l - 1 + L) % L));
    }
    E = E / 2; // since we double count each interaction

    // std::cout << "Initial energy: " << E << std::endl;
    return config;
}

// arma::mat find_interacting_pairs(int k, int l, arma::mat config, int L)
//{
//     // create a vector of vectors to store the interacting pairs
//     arma::mat interacting_pairs = arma::zeros<arma::mat>(4, 2); // rows are the interacting pairs, columns are the coordinates of the interacting pairs
//
//     // notice the intentional use of integer division
//
//     int index_value = config(k, l);
//     interacting_pairs.col(0) = arma::vec(4).fill(index_value);
//     //  there must be a better way to do this
//     interacting_pairs(0, 1) = config(k, (l + 1) % L);
//     interacting_pairs(1, 1) = config(k, (l - 1 + L) % L);
//     interacting_pairs(2, 1) = config((k + 1) % L, l);
//     interacting_pairs(3, 1) = config((k - 1 + L) % L, l);
//
//     return interacting_pairs;
// }

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

        neigbours_energy += -config(k, l) * (config((k + 1) % L, l) + config((k - 1 + L) % L, l) + config(k, (l + 1) % L) + config(k, (l - 1 + L) % L));

        delta_E = -2 * neigbours_energy; // because because delta_E = E_new - E_old = -((-a)*b) - (-(a*b)) = 2*(a*b) = -2*neigbours_energy

        if (delta_E <= 0)
        {
            config(k, l) *= -1;
            E += delta_E;
            M += 2 * config(k, l); // factor of 2 because we are only changing one spin
        }
        else
        {
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

void monte_carlo(int L, int mc_cycles, int burn_pct, double T, int align, int thread_id)
{
    double beta = 1. / (kb * T);
    arma::vec p_vec = arma::zeros<arma::vec>(5);
    arma::vec delta_e_vec = arma::vec(std::vector<double>{-8, -4, 0, 4, 8});
    for (int i = 0; i < 5; i++)
    {
        p_vec(i) = std::exp(-beta * delta_e_vec(i));
    }

    std::ofstream file2;

    std::string rounded_T = std::to_string(T).substr(0, std::to_string(T).find(".") + 3 + 1);

    std::cout
        << "T: " << T << std::endl;
    std::string filename1 = std::to_string(L) + "/s1_cfg_L" + std::to_string(L) + "_A" + std::to_string(align) + "_mc" + std::to_string(mc_cycles) + "_burn" + std::to_string(burn_pct) + "_t" + rounded_T + ".csv";
    std::string filename2 = std::to_string(L) + "/s2_cfg_L" + std::to_string(L) + "_A" + std::to_string(align) + "_mc" + std::to_string(mc_cycles) + "_burn" + std::to_string(burn_pct) + "_t" + rounded_T + ".csv";
    std::string filename3 = std::to_string(L) + "/s3_cfg_L" + std::to_string(L) + "_A" + std::to_string(align) + "_mc" + std::to_string(mc_cycles) + "_burn" + std::to_string(burn_pct) + "_t" + rounded_T + ".csv";
    //std::string filename2 = std::to_string(L) + "/_TEST_qt_L" + std::to_string(L) + "_A" + std::to_string(align) + "_mc" + std::to_string(mc_cycles) + "_burn" + std::to_string(burn_pct) + "_t" + rounded_T + ".csv";

    // file1.open(filename);

    // avg will always be wrt to the number of mc cycles
    double E = 0;
    double M = 0;

    double cumul_E = 0;

    double cumul_E2 = 0;
    double avg_e = 0;
    double avg_E = 0;
    double avg_E2 = 0;

    double var_E = 0;

    double var_M = 0;
    double avg_M2 = 0;

    double avg_mabs = 0;

    double cumul_M2 = 0;
    double avg_Mabs = 0;

    double avg_Mabs2 = 0;
    double cumul_Mabs = 0;

    double cumul_mabs = 0;

    double Cv = 0;
    double chi = 0;
    arma::mat config = init_random_config(L, E, M, align);

    std::ofstream file_s1;
    std::ofstream file_s2;
    std::ofstream file_s3;

    file_s1.open(filename1);
    file_s2.open(filename2);
    file_s3.open(filename3);
    config.save(file_s1, arma::csv_ascii);
    file_s1.close();

    double N = L * L;
    int time_steps = mc_cycles * N;
    int burn_in = int((burn_pct / 100.) * time_steps);
    int mc_counter = 0;

    arma::mat output_data = arma::zeros<arma::mat>(mc_cycles - 1, 4);
    arma::mat output_state = arma::zeros<arma::mat>(L, L);
    // std::cout << "burned mc_cycles: " << burn_in << std::endl;

    // std::cout << "time_steps = " << time_steps << std::endl;
    for (int i = 0; i < time_steps; i++)
    {
        config = evolve(config, beta, E, M, p_vec);
        // the following qtd will be calculated in a dumb way just to be didatic but uses lots of flops

        cumul_E += E;
        cumul_E2 += E * E;
        cumul_Mabs += std::abs(M);
        cumul_M2 += M * M;
        // std::cout << "E = " << E << std::endl;

        if (i > burn_in)
        {
            if (i % int(N) == 0 && i != 0) // output only at the end of a mc_cycle
            {
                
                
                

                mc_counter += 1;

                if (mc_counter==int(mc_cycles/2))
            {
                config.save(file_s2, arma::csv_ascii);
                file_s2.close();
            }
            }
        }
    }
    config.save(file_s3, arma::csv_ascii);
    file_s3.close();
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
