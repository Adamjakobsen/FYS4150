#include "./include/Schrodinger.hpp"
#include <iostream>
#include <iomanip>
#include <armadillo>

double v0=std::pow(10,10)*1.;// potential
double h=0.005; // step size
double dt=std::pow(10,-5)*2.5; // step size time


double T=0.008; // total time
double x_0=0.25; // initial position x
double sigma_x=0.05; // initial width x
double px=200; // initial momentum x
double y_0=0.5; // initial position y
double sigma_y=0.1; // initial width y
double py=0; // initial momentum y

int M = int(1.0/h);// number of steps
int N = M-2; // number of inner steps


Schrodinger sc=Schrodinger(v0,h,dt,M);

int main(){

arma::cx_mat U0;

std::string slits="double";
std::cout<<"dt:"<<dt<<std::endl;
sc.set_potential(slits);
std::cout << "V set" << std::endl;
sc.set_A_B();
std::cout << "AB set" << std::endl;
U0=sc.initialise_state(sigma_x,sigma_y,x_0,y_0,px,py);
std::cout << "State initialised" << std::endl;

double T_dt=T/dt;

int n_timesteps = (int)(T/dt);

arma::cx_mat U = arma::cx_mat(N*N,n_timesteps);
//print shape of U
std::cout << "U shape: " << U.n_rows << " " << U.n_cols << std::endl;



U.col(0)=sc.get_u_vec(U0);
std::cout << "U(0) set" << std::endl;

for (int i=1;i<n_timesteps;i++){
    U.col(i)=sc.evolve(U.col(i-1));
}
string path="./data/";
string filename="U_"+slits+".bin";
U.save(path+filename);

return 0;
}

