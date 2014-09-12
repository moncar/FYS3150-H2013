#include <armadillo>

void forwardEulerDiffusion(double dt, double T, double dx, arma::colvec ic, arma::colvec bc, std::string filepath);
void backwardEulerDiffusion(double dt, double T, double dx, arma::colvec ic, arma::colvec bc, std::string filepath);
void crankNicolsonDiffusion(double dt, double T, double dx, arma::colvec ic, arma::colvec bc, std::string filepath);
