#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>

int main(int argc, char *argv[])
{
    // Number of steps
    int n_steps = atoi(argv[1]);

    // Define number of points we want to solve for
    int n = n_steps + 1;
    // Number of points in complete solution
    int m = n + 2;
    // Set parameters and boundary values
    double x_min = 0.0;
    double x_max = 1.0;
    double h = (x_max - x_min) / n_steps;

    double v0 = 0.0;
    double vm = 0.0;

    // Define vector for complete solution and fill in boudary values
    std::vector<double> v(n + 1);
    v[0] = v0;
    v[n - 1] = vm;

    // Main diagonal vector a of length n
    std::vector<double> a(n, -1.0);
    // sub/superdiagonals with length n-1
    std::vector<double> b(n, 2.0);
    std::vector<double> c(n, -1.0);
    // Defining vector g
    std::vector<double> g(n);
    std::vector<double> g_tilde(n);
    // Define vector x
    std::vector<double> x(n);
    x[0] = x_min;

    double exp(double x);
    for (int i = 0; i < n_steps + 1; i++)
    {
        x[i] = i * h;
        g[i] = h * h * 100 * exp(-10 * x[i]);
    }
    double a_b;

    g_tilde[1] = g[1];                    /////////// here starts the difference between problem 7 and problem 9
    std::vector<double> pre_comp_b_tilde; // pre computed b tilde vectors
    for (int i = 1; i < n + 1, i++)       // notice it has the same dimensions as b
    {
        pre_comp_b_tilde.push_back((i + i) / i)
    }
    std::cout << "ASSERT" << b.size() << "/" << pre_comp_b_tilde.size() << std::endl;

    for (int i = 2; i < n_steps; i++)
    {
        g_tilde[i] = g[i] + (g_tilde[i - 1] / pre_comp_b_tilde[i - 1]);
    }
    v[n_steps - 1] = g_tilde[n - 1] / pre_comp_b_tilde[n - 1];

    for (int i = n - 2; i > 0; i--)
    {
        v[i] = (g_tilde[i] - c[i] * v[i + 1]) / pre_comp_b_tilde[i];
    }

    // Set filename
    std::string filename = "output_problem9_" + std::to_string(n_steps) + ".txt";

    // Create output file stream
    std::ofstream ofile;
    ofile.open(filename);

    // Format parameters
    int width = 18;
    int prec = 10;

    // Loop and write to file
    for (int i = 0; i < n_steps + 1; i++)
    {
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x[i]
              << std::setw(width) << std::setprecision(prec) << std::scientific << v[i] << std::endl;
    }
    // Close the output file
    ofile.close();
    return 0;
}
