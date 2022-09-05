#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>


int main(int argc, char* argv[])
{
//Number of steps
int n_steps= atoi(argv[1]);

//Define number of points we want to solve for
int n= n_steps + 1;
//Number of points in complete solution
int m=n+2;
//Set parameters and boundary values
double x_min = 0.0;
double x_max = 1.0;
double h = (x_max - x_min) / n_steps;

double v0=0.0;
double vm=0.0;

//Define vector for complete solution and fill in boudary values
std::vector<double> v(m);
v[0]=v0;
v[m-1]=vm;


//Main diagonal vector a of length n
std::vector<double> a(n , 2.0);
//sub/superdiagonals with length n-1
std::vector<double> b(n-1 , -1.0);
std::vector<double> b_tilde(n-1);
std::vector<double> c(n-1 , -1.0);
//Defining vector g
std::vector<double> g(n);
std::vector<double> g_tilde(n);
//Define vector x
std::vector<double> x(m);
x[0]=x_min;



for (int i=0;i<=n_steps;i++)
{
	x[i]= i*h;	
	g[i]= h*h*100*exp(-10*x[i]);
}
double a_b;
b_tilde[1]=b[1];
g_tilde[1]=g[1];
for (int i=2; i<=n;i++)
{
	a_b =  a[i]/b_tilde[i-1];
	b_tilde[i]=b[i] - a_b*c[i-1];
	g_tilde[i]=g[i] - a_b*g_tilde[i-1];
}
v[n-1]=g_tilde[n-1]/b_tilde[n-1];

for (int i=n-1;i>=1;i--)
{
v[i]=(g_tilde[i] - c[i]*v[i+1])/b_tilde[i];
}

//Set filename
std::string filename = "General_solution_"+std::to_string(n_steps)+".txt";

//Create output file stream
std::ofstream ofile;
ofile.open(filename);

//Format parameters
int width = 14;
int prec = 4;


//Loop and write to file
for (int i = 0; i<=n_steps;i++)
{
ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x[i]
	<<std::setw(width) << std::setprecision(prec) << std::scientific << v[i] <<std::endl;

}
//Close the output file
ofile.close();
return 0;
}
