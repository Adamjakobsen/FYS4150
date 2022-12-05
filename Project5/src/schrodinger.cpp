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
    int width_walls=pos_to_idx(0.02);
    V.cols(0,width_walls).fill(v0);
    
    V.cols(M-1-width_walls ,M-1)= arma::ones<arma::cx_mat>(M,width_walls+1)*v0;
    V.rows(0,width_walls) = arma::ones<arma::cx_mat>(width_walls+1,M)*v0;
    V.rows(M-1-width_walls ,M-1)= arma::ones<arma::cx_mat>(width_walls+1,M)*v0;
    cout<<"v0:"<<v0 <<endl;
    V.rows(pos_to_idx(0.5-0.01),pos_to_idx(0.5+0.01))=arma::ones<arma::cx_mat>(width_walls+1,M)*v0;


    int start_x=pos_to_idx(0.5-0.01);
    int end_x=pos_to_idx(0.5+0.01);
    
    if (slits=="double")
    
    {   
        int left_start_y=pos_to_idx(0.5-0.075);
        int left_end_y=pos_to_idx(0.5-0.025);

        int right_start_y=pos_to_idx(0.5+0.025);
        int right_end_y=pos_to_idx(0.5+0.075)+1;

        V.submat(arma::span(start_x,end_x), arma::span(left_start_y,left_end_y))    =arma::zeros<arma::cx_mat>(end_x-start_x+1,left_end_y-left_start_y+1); 
        V.submat(arma::span(start_x,end_x), arma::span(right_start_y,right_end_y))  =arma::zeros<arma::cx_mat>(end_x-start_x+1,right_end_y-right_start_y+1); 
        
        
    }
    
    //make slit in middle wall
    if (slits=="single")
    {   
        int middle_start_y=pos_to_idx(0.5-0.025);
        int middle_end_y=pos_to_idx(0.5+0.025);
        V.submat(arma::span(start_x,end_x), arma::span(middle_start_y,middle_end_y))=arma::zeros<arma::cx_mat>(end_x-start_x+1,middle_end_y-middle_start_y+1); 

        
        

    
    }

    
    //Triple slit
    if (slits=="triple")
    {
    

    int middle_start_y=pos_to_idx(0.5-0.025);
    int middle_end_y=pos_to_idx(0.5+0.025);

    int left_start_y=pos_to_idx(0.5-0.075-0.05);
    int left_end_y=pos_to_idx(0.5-0.075)+1;

    int right_start_y=pos_to_idx(0.5+0.075);
    int right_end_y=pos_to_idx(0.5+0.075+0.05);

    V.submat(arma::span(start_x,end_x), arma::span(middle_start_y,middle_end_y))=arma::zeros<arma::cx_mat>(end_x-start_x+1,middle_end_y-middle_start_y+1); 
    V.submat(arma::span(start_x,end_x), arma::span(left_start_y,left_end_y))    =arma::zeros<arma::cx_mat>(end_x-start_x+1,left_end_y-left_start_y+1); 
    V.submat(arma::span(start_x,end_x), arma::span(right_start_y,right_end_y))  =arma::zeros<arma::cx_mat>(end_x-start_x+1,right_end_y-right_start_y+1); 
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
        ak = 1. + 4.*r +1i*dt*V(i,j)/2.;
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
    A.save("./data/A.bin");
    cout<<"A set"<<endl;
    
             
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
