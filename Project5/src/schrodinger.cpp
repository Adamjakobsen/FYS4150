#include "../include/Schrodinger.hpp"
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <complex>
using namespace std;





Schrodinger::Schrodinger(double v0_in, double h_in,double dt_in,int M_in)
{
    v0 = v0_in;
    h = h_in;
    dt = dt_in;
    M = M_in;
    
    arma::cx_mat V ;
    arma::mat A ;
    arma::mat B ;
}

int Schrodinger::ij_k(int i, int j) {
    // converts 2D indices to 1D
    
    return i*(M-2) + j;
}
int Schrodinger::pos_to_idx(double pos) {
    //positon to index
    int  index=int(pos/h);
    
    
    return index;
}

void Schrodinger::set_potential(string slits)
{
    
    // set potential
    
    V = arma::zeros<arma::cx_mat>(M,M);
    
    V.col(0)=arma::cx_vec(M).fill(v0);
    
    V.col(M-1).fill(v0);
    V.row(0).fill(v0);
    V.row(M-1).fill(v0);
    cout<<"v0:"<<v0 <<endl;
    V.row(pos_to_idx(0.5)).fill(v0);
    
    if (slits=="double")
    
    {
        cout<<"double slits"<<endl;
        for (int i=0;i<M;i++){
        if ( (i>=pos_to_idx(0.5-0.075)) && (i<=pos_to_idx(0.5-0.05/2))){

            for (int j=pos_to_idx(0.5-0.01);j<=pos_to_idx(0.5+0.01);j++){
                V(j,i)=0;
            }
        }
        if ( (i>=pos_to_idx(0.5+0.05/2)) && (i<=pos_to_idx(0.5+0.05+0.05/2)) ){
            for (int j=pos_to_idx(0.5-0.01);j<=pos_to_idx(0.5+0.01);j++){
                V(j,i)=0;
            }
            
        }
    }
    }
    //make slit in middle wall
    if (slits=="single")
    {
        for (int i=0;i<M;i++){
        if (i>=pos_to_idx(0.5-0.05/2) && i<=pos_to_idx(0.5+0.05/2)){
            
            for (int j=pos_to_idx(0.5-0.01);j<=pos_to_idx(0.5+0.01);j++){
                V(j,i)=0;
            }
        }
    }
    }
    V.save("./data/V.bin");
    
        
    

}

void Schrodinger::set_A_B()
{   
    
    cout<<"set A_B"<< (M-2)<<"M:"<<M <<endl;
    A = arma::sp_cx_mat((M-2)*(M-2),(M-2)*(M-2));
    B = arma::sp_cx_mat((M-2)*(M-2),(M-2)*(M-2));

    
    complex<double> ak;
    complex<double> bk;
    complex<double> r=1i*dt/(2*h*h);
    
    int i;
    int j;
    for (int k=0; k<(M-2)*(M-2); k++){
        
        i=k/(M-2);
        j=k%(M-2);
        ak = 1. + 4.*r -1i*dt*V(i,j)/2.;
        bk = 1. - 4.*r -1i*dt*V(i,j)/2.;
        A(k,k)=ak;
        B(k,k)=bk;
         
        if (i!=0){
            A(k , ij_k(i-1,j))=-r;
            A(ij_k(i-1,j) , k)=-r;

            B(k , ij_k(i-1,j))=r;
            B(ij_k(i-1,j) , k)=r;
        }
        if (j!=0){
            A(k , ij_k(i,j-1))=-r;
            A(ij_k(i,j-1) , k)=-r;

            B(k , ij_k(i,j-1))=r;
            B(ij_k(i,j-1) , k)=r;

        }
        
    }        
}

arma::cx_vec Schrodinger::evolve(arma::cx_vec u_vec){
    arma::superlu_opts opts;
    opts.symmetric = true;
    // evolve wavefunction
    arma::cx_vec u_vec_new;
    u_vec_new = arma::zeros<arma::cx_vec>(u_vec.n_elem);
    u_vec_new = arma::spsolve(A,B*u_vec, "superlu",opts);
    return u_vec_new;
}

complex<double> Schrodinger::gaussian(double x, double y, double sigma_x, double sigma_y, double x0, double y0,double px,double py){
    // gaussian wavepacke
    complex<double> psi;
    psi = exp(-pow(x-x0,2)/(2*pow(sigma_x,2)) - pow(y-y0,2)/(2*pow(sigma_y,2)))*exp(1i*(px*(x-x0) + py*(y-y0)));
    return psi;
}

arma::cx_mat Schrodinger::initialise_state(double sigma_x, double sigma_y, double x0, double y0,double px,double py){
    // initialize wavefunction
    arma::cx_mat U;
    U = arma::zeros<arma::cx_mat>(M-2,M-2);
    for (int i=0;i<M-2;i++){
        for (int j=0;j<M-2;j++){
           
            U(i,j) = gaussian(i*h,j*h,sigma_x,sigma_y,x0,y0,px,py);
        }
    }

    /* //Boundary conditions
    U.row(0).fill(0);
    U.row(M-1).fill(0);
    U.col(0).fill(0);
    cout<<"0 col set"<<endl;
    U.col(M-1).fill(0);
    cout<<"M-1 col set"<<endl;
     */
    //Normalize
    U = U/sqrt(arma::accu(arma::conj(U)%U));
    cout << "U set and normalised" << endl;
    return U;
}

arma::cx_vec Schrodinger::get_u_vec(arma::cx_mat U){
    // convert 2D matrix to 1D vector
    arma::cx_vec u_vec;
    u_vec = arma::zeros<arma::cx_vec>((M-2)*(M-2));
    for (int i=0;i<M-2;i++){
        for (int j=0;j<M-2;j++){
            u_vec(ij_k(i,j)) = U(i,j);
        }
    }
    return u_vec;
}
