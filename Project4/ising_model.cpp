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
        int iter_;
        double T;
        double temp_step = atof(argv[6]);
#pragma omp for
        for (iter_ = 0; iter_ < temp_n + 1; iter_++)
        {
            T = lower_temp + iter_ * temp_step;
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

    return config;
}

arma::mat evolve(arma::mat config, double beta, double &E, double &M, arma::vec p_vec)
{
    int L = config.n_rows;
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

    std::ofstream file1;

    std::string rounded_T = std::to_string(T).substr(0, std::to_string(T).find(".") + 3 + 1);

    std::cout << "T: " << T << std::endl;
    std::string filename1 = std::to_string(L) + "/qt_L" + std::to_string(L) + "_A" + std::to_string(align) + "_mc" + std::to_string(mc_cycles) + "_burn" + std::to_string(burn_pct) + "_t" + rounded_T + ".csv";

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

    double N = L * L;
    int time_steps = mc_cycles * N;
    int burn_in = int((burn_pct / 100.) * mc_cycles);
    int mc_counter = 0;

    arma::mat output_data = arma::zeros<arma::mat>(mc_cycles, 4);
    for (int i = 0; i < mc_cycles; i++)
    {
        config = evolve(config, beta, E, M, p_vec);

        cumul_E += E;
        cumul_E2 += E * E;
        cumul_Mabs += std::abs(M);
        cumul_M2 += M * M;

        if (i >= burn_in)
        {

            avg_E2 = cumul_E2 / (i + 1);
            avg_E = cumul_E / (i + 1);
            var_E = avg_E2 - avg_E * avg_E;

            avg_Mabs = cumul_Mabs / (i + 1);
            cumul_mabs = cumul_Mabs / N;

            avg_M2 = cumul_M2 / (i + 1);
            avg_Mabs2 = avg_Mabs * avg_Mabs;

            var_M = avg_M2 - avg_Mabs2;

            avg_e = avg_E / N;
            avg_mabs = avg_Mabs / N;

            Cv = (var_E) / (N * kb * T * T);
            chi = (var_M) / (N * kb * T);

            output_data(i, 0) = avg_e;
            output_data(i, 1) = avg_mabs;
            output_data(i, 2) = Cv;
            output_data(i, 3) = chi;
        }
    }
    file1.open(filename1);
    output_data.save(file1, arma::csv_ascii);
    file1.close();
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