#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

// Some width and precision parameters we will use to format the output
int width = 12;
int prec = 4;

string base_filename = "output_problem8_";

double solution(double x)
{

    return 1 - (1 - exp(-10)) * x - exp(-10 * x);
}

vector<double> ThomasAlgorithm(double n, vector<double> g, vector<double> a, vector<double> b, vector<double> c) // sub, main, super diagonal respectively
{
    double m = n - 2; // m is the number of solution points the algorithm gives us

    // vector<double> g_tilda, b_tilda;
    vector<double> v(m + 1);

    double g_tilda_i, b_tilda_i, v_i;

    // parameters limit cases
    // g_tilda.push_back(g[0]);
    // b_tilda.push_back(b[0]);

    // for loop of forward substitution
    for (int i = 1; i < m + 1; i++) // remember g has only m = n - 2 elements and we already have g_0
    {
        double omega = a[i] / b[i - 1];
        g_tilda_i = g[i] - (omega * g[i - 1]);
        b_tilda_i = b[i] - (omega * c[i - 1]);

        g[i] = g_tilda_i;
        b[i] = b_tilda_i;
    }

    // solutions limit cases
    v[m] = (g_tilda_i / b_tilda_i); // since loop is over, push last value instead vector's last entry

    // for loop of back substitution
    for (int i = m - 1; i >= 0; i--) // m - 1 because we are starting at 0
    {
        v_i = (g[i] - (c[i] * v[i + 1])) / b[i];
        v[i] = v_i;
    }

    return v;
}

vector<double> source_function(double n, vector<double> x_vec)
{
    vector<double> f;
    for (int i = 0; i < n; i++) // notice there will be n iterations
    {
        f.push_back(100 * exp(-10 * x_vec[i]));
    }
    return f;
}

int main()
{
    vector<double> steps_vector{10, 100, 1000, 10000, 100000, 1000000, 10000000};
    for (double n : steps_vector) // n is for number of steps
    {

        double x = 0;
        double x_min = 0;
        double x_max = 1;

        vector<double> u_vec;
        vector<double> x_vec;

        double step = (x_max - x_min) / n;

        for (int i = 0; i < n + 1; i++) // notice there will be n + 1 steps
        {
            x_vec.push_back(x); // notice the order matters for the boundary
            x = x + step;
        }

        vector<double> f = source_function(n, x_vec);

        double m = n - 2;
        vector<double> a(m, -1.); // # of elements is less than that of the main diag
        vector<double> c(m, -1.);
        vector<double> b(m + 1, 2.);
        vector<double> g(m + 1);

        // limit cases for g
        g[0] = step * step * f[x_vec[1]];     // + v[0] but recall v[0] = u[0] = 0
        g[m] = step * step * f[x_vec[m + 1]]; // + v[n-1] but recall v[n-1] = u[n-1] = 0
        for (int i = 1; i < m; i++)
        {
            g[i] = step * step * f[i + 1];
        }

        vector<double> v = ThomasAlgorithm(n, g, a, b, c); // notice that this vector has m entries but is it a good idea to fill in the boundary terms

        // fill in vector v
        v.insert(v.begin(), 0); // pushes 0 to the front of vector
        v.push_back(0);         // pushes 0 to the end of vector

        // OUTPUT TO FILE
        ofstream ofile;
        string filename = base_filename + to_string(int(n)) + ".txt";
        ofile.open(filename);
        for (int i = 1; i < n; i++) // notice there will be m iterations
        {
            double diff = v.at(i) - solution(x_vec.at(i));
            ofile
                << setw(width) << setprecision(prec) << scientific << x_vec.at(i)
                << ", "
                << setw(width) << setprecision(prec) << scientific << v.at(i)
                << ", "
                << setw(width) << setprecision(prec) << scientific << solution(x_vec.at(i))
                << ", "
                << setw(width) << setprecision(prec) << scientific << diff
                << endl;
        }
        ofile.close();
    }
    // all is well, return 0
    return 0;
}
