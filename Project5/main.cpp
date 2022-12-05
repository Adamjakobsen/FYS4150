#include "./include/Schrodinger.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

int main(int argc, char *argv[])
{
    // Declaring variables
    arma::cx_mat U0;
    double h, dt, T, x_0, sigma_x, px, y_0, sigma_y, py;
    std::string slits;
    double v0 = std::pow(10, 10) * 1.; // potential

    // Input filename from command line
    std::string input_filename = argv[1];
    std::string path_input = "./input/";

    // create filestream to read in inputs from file
    std::ifstream infile;
    infile.open(path_input + input_filename + ".txt");

    if (!infile.is_open())
    {
        std::cout << "Error opening file" << std::endl;
        return 1;
    }

    // read in inputs from file

    std::string line;

    while (std::getline(infile, line))
    {
        // loop through lines and sett corect variables. name and value separated by "="
        std::string name = line.substr(0, line.find("="));
        std::string value = line.substr(line.find("=") + 1);
        if (name == "h")
        {
            h = std::stod(value);
        }
        else if (name == "dt")
        {
            dt = std::stod(value);
        }
        else if (name == "T")
        {
            T = std::stod(value);
        }
        else if (name == "x_0")
        {
            x_0 = std::stod(value);
        }
        else if (name == "sigma_x")
        {
            sigma_x = std::stod(value);
        }
        else if (name == "px")
        {
            px = std::stod(value);
        }
        else if (name == "y_0")
        {
            y_0 = std::stod(value);
        }
        else if (name == "sigma_y")
        {
            sigma_y = std::stod(value);
        }
        else if (name == "py")
        {
            py = std::stod(value);
        }
        else if (name == "slits")
        {
            slits = value;
        }
        else
        {
            std::cout << "Error reading file" << std::endl;
            return 1;
        }
    }
    // print variables
    std::cout << "Running configurat from " << input_filename + " :" << std::endl;
    std::cout << "h=" << h << std::endl;
    std::cout << "dt=" << dt << std::endl;
    std::cout << "T=" << T << std::endl;
    std::cout << "x_0=" << x_0 << std::endl;
    std::cout << "sigma_x=" << sigma_x << std::endl;
    std::cout << "px=" << px << std::endl;
    std::cout << "y_0=" << y_0 << std::endl;
    std::cout << "sigma_y=" << sigma_y << std::endl;
    std::cout << "py=" << py << std::endl;
    std::cout << "slits=" << slits << std::endl;
    std::cout << "" << std::endl;

    int M = int(1.0 / h) + 1; // number of points
    int N = M - 2;            // number of inner points

    Schrodinger sc = Schrodinger(v0, h, dt, M);

    sc.set_potential(slits);

    sc.set_A_B();

    U0 = sc.initialise_state(sigma_x, sigma_y, x_0, y_0, px, py);
    std::cout << "State initialised" << std::endl;
    std::cout << "" << std::endl;

    double T_dt = T / dt;

    int n_timesteps = (int)(T / dt);

    // Define matrix: holds the state as vector at each timestep
    arma::cx_mat U = arma::cx_mat(N * N, n_timesteps);

    U.col(0) = sc.get_u_vec(U0);

    for (int i = 1; i < n_timesteps; i++)
    {
        U.col(i) = sc.evolve(U.col(i - 1));
    }
    // Save U to file
    string path = "./data/";
    string filename = "U_" + input_filename + ".bin";
    U.save(path + filename);
    std::cout << "data saved to " << path + filename << std::endl;
    std::cout << "" << std::endl;
    std::cout << "======================================================================" << std::endl;
    return 0;
}
