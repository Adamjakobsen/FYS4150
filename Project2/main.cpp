#include <iostream>
#include <armadillo>

int main(){
int N = 6;
int n = N + 1;
double h = 1./n;
//initialising v that I will be solving for
//giving it a dummy first index
//it lacks the last boundary condition index
//fill in vector with 1's just to know what's up 
arma::vec v = arma::vec(N).fill(1.);
v(0) = 0.;
// std::cout << v_star;
double a = -1/(h*h);
double b = 2/(h*h);
arma::mat A = arma::mat(N,N).fill(0.);
A.diag() = arma::vec(N).fill(b);
A.diag(1) = arma::vec(N-1).fill(a); 
A.diag(-1) = arma::vec(N-1).fill(a);
arma::vec eigval;
arma::mat eigvec;
eig_sym(eigval, eigvec, A);
//scaling of eigenvectors is not fixed yet; gives some
//not important yet, only needed for presentation at the end
arma::vec norm_eigvec = normalise(eigvec);
std::cout << norm_eigvec << '\n' ;


int check_col = 3;
arma::vec check1 = A*eigvec.col(check_col);
arma::vec check2 = eigval(check_col)*eigvec.col(check_col);
arma::vec final = check1 - check2;
// std::cout << final << '\n';


}


