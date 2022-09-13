

#include <iostream>
#include <vector>
#include <cmath>

// double plotter(double x); // Declaration of plotter func

// the section creating a textfile was inspired by a code on stackoverflow https://stackoverflow.com/questions/478075/creating-files-in-c
#include <fstream>

double calc(double x);

void text()
{

    std::ofstream outfile("data.txt");

    for (double i = 0.0; i <= 1000.0; ++i)
    {
        double var = 0.0 + i / 1000;

        outfile << var << " " << calc(var) << std::endl;
    }
    outfile.close();
}

int main()
{

    text();

    return 0;
}

double calc(double x)
{
    double u;
    u = 1 - (1 - exp(-10)) * x - exp(-10 * x);
    return u;
}
