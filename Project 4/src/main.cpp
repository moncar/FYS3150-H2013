#include <iostream>
#include <armadillo>

#include "diffusionIntegration.hpp"

using namespace std;
using namespace arma;

// Time constants
const double T = 1;

int main(int argc, char* argv[]) {
  if (argc != 5) {
    cout << "Usage: <exe> 1/2/3 dx dt savefile" << endl << "1: FE, 2: BE, 3: CN" << endl;
    return 1;
  }

  int method = atoi(argv[1]);
  double dx = atof(argv[2]);
  double dt = atof(argv[3]);
  string savefile = argv[4];

  // Initial condition
  int nx = (1.0 / dx) + 1;
  colvec ic = zeros(nx);

  // Boundary conditions
  colvec bc = zeros(2); bc(0) = 1; bc(1) = 0;

  switch (method) {
    case 1:
      forwardEulerDiffusion(dt,T,dx,ic,bc,savefile);
      break;
    case 2:
      backwardEulerDiffusion(dt,T,dx,ic,bc,savefile);
      break;
    case 3:
      crankNicolsonDiffusion(dt,T,dx,ic,bc,savefile);
      break;
  }

  cout << "Reached end of main." << endl;
  return 0;
}
