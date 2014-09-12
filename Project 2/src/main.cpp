#include <stdlib.h>
#include <armadillo>

#include "JacobiRotationProblem.hpp"
#include "potentials.hpp"

using namespace std;
using namespace arma;

const double K = 1; // Boltzmann constant?
const double M_ELECTRON = 9.10938291e-31; // [kg]
const double HBAR = 1.05457173e-34; // [m^2 kg / s]
const double ALPHA = pow( (HBAR * HBAR) / (M_ELECTRON * K), 0.25);
const double C_LAMBDA = (2.0 * M_ELECTRON * ALPHA * ALPHA ) / (HBAR * HBAR);

void radialSchrodingerOneElectron(double,int,string);
void radialSchrodingerTwoElectrons(double,int,double,string);
void testCaseSymmetricMatrix();

int main(int argc, char *argv[]) {
  // Print usage if lacking certain number of electrons
  string usage = "Usage: ./<exe>.x -#e <#electrons> -n <Size matrix> -wr <omega r> -rhoMax <max radius> -o <savefilename>";
  if (argc < 2) {
    cout << usage << endl;
    return 1;
  }

  // For storing commandline arguments
  int electrons = 0;
  int n = 0;
  double omegaR = 0;
  string savefile = " ";
  double rhoMax = 0;

  // Traverse commandline and store
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i],"-n") == 0) {
      n = atoi(argv[i+1]);
      i++;
    } else if (strcmp(argv[i],"-#e") == 0) {
      electrons = atoi(argv[i+1]);
      i++;
    } else if (strcmp(argv[i],"-wr") == 0) {
      omegaR = atof(argv[i+1]);
      i++;
    } else if (strcmp(argv[i],"-rhoMax") == 0) {
      rhoMax = atof(argv[i+1]);
      i++;
    } else if (strcmp(argv[i],"-o") == 0) {
      savefile = argv[i+1];
      i++;
    }
  }

  // Check that necessary are present
  if (n == 0) {
    cout << "Lacking size of matrix. Not optional arg." << endl;
    cout << usage << endl;
    return 1;
  }
  if (electrons == 0) {
    cout << "Lacking number of electrons. Not optional arg." << endl;
    return 1;
  }
  if (rhoMax == 0) {
    cout << "Lacking rhoMax. Not optional arg." << endl;
    return 1;
  }
  if (strcmp(savefile.c_str()," ") == 0) {
    cout << "Lacking output filename. Not optional arg." << endl;
    return 1;
  }

  // Run correct simulation
  if (electrons == 1) {
    radialSchrodingerOneElectron(rhoMax,n,savefile);
  } else if (electrons == 2) {
    if (omegaR == 0) {
      cout << "Can't run for two electrons without strength of oscillator." << endl;
      return 1;
    }
    radialSchrodingerTwoElectrons(rhoMax,n,omegaR,savefile);
  }

  return 0;
}

/*
 * Uses the Jacobi eigenvalue algorithm to find the eigenvalues for the radial
 * Schrodinger equation for one electron in a three-dimensional harmonic
 * oscillator well. Calculation is done using dimensionless variables.
 *
 * @param rhoMax Dimensionless max radius.
 * @param n Size of matrix.
 * @param savefile File to save eigenvalues to as sorted list.
 */
void radialSchrodingerOneElectron(double rhoMax, int n, string savefile) {
  double rhoMin = 0;                            // Minimum dimensionless radius
  double h = (rhoMax - rhoMin) / (double) (n+1);  // Step size
  double E = -1.0 / (h*h);                    // Non-diagonal elements value

  vec rho = zeros<vec>(n);                    // Discrete radii
  vec V = zeros<vec>(n);                    // Discrete potential
  vec d = zeros<vec>(n);                    // Diagonal matrix elements
  vec e = zeros<vec>(n-1);                      // Non-diagonal matrix elements

  for (int i = 0; i < n; i++) {
    rho(i) = rhoMin + (i+1)*h;
    V(i) = harOscV(rho(i));
    d(i) = 2.0 / (h*h) + V(i);
    if (i < n-1) { e(i) = E; }
  }

  // Create the matrix
  mat problem = zeros<mat>(n,n);
  problem.diag() = d;
  problem.diag(1) = e;
  problem.diag(-1) = e;

  // Solve with Jacobis method
  JacobiRotationProblem prob = JacobiRotationProblem(problem);
  prob.solve(1e-8);
  prob.giveParameters(rhoMax,1);
  prob.saveResult(savefile);
}

/*
 * Same as above problem only now with two electrons and including the Coulomb
 * potential.
 *
 * @param rhoMax Dimensionless max radius.
 * @param n Size of matrix.
 * @param omegaR Parameter decing strength of harmonic oscillator potential.
 * @param savefile File to save eigenvalues to as sorted list.
 */
void radialSchrodingerTwoElectrons(double rhoMax, int n, double omegaR, string savefile) {
  double rhoMin = 0;                            // Minimum dimensionless radius
  double h = (rhoMax - rhoMin) / (double) (n+1);  // Step size
  double E = -1.0 / (h*h);                    // Non-diagonal elements value

  vec rho = zeros<vec>(n);                    // Discrete radii
  vec V = zeros<vec>(n);                    // Discrete potential
  vec d = zeros<vec>(n);                    // Diagonal matrix elements
  vec e = zeros<vec>(n-1);                      // Non-diagonal matrix elements

  for (int i = 0; i < n; i++) {
    rho(i) = rhoMin + (i+1)*h;
    V(i) = harOscTwoElectronsV(rho(i),omegaR);
    d(i) = 2.0 / (h*h) + V(i);
    if (i < n-1) { e(i) = E; }
  }

  // Create the matrix
  mat problem = zeros<mat>(n,n);
  problem.diag() = d;
  problem.diag(1) = e;
  problem.diag(-1) = e;

  // Solve with Jacobis method
  JacobiRotationProblem prob = JacobiRotationProblem(problem);
  prob.solve(1e-8);
  prob.giveParameters(rhoMax,2,omegaR);
  prob.saveResult(savefile);
}

/*
 * Contains a simple test case. Eigenvectors should be 8, 3 and 6.
 */
void testCaseSymmetricMatrix() {
  mat mymat = randu<mat>(3,3);

  mymat(0,0) = 6;
  mymat(0,1) = -2;
  mymat(0,2) = -1;
  mymat(1,0) = -2;
  mymat(1,1) = 6;
  mymat(1,2) = -1;
  mymat(2,0) = -1;
  mymat(2,1) = -1;
  mymat(2,2) = 5;

  cout << "Matrix:" << endl << mymat << endl;

  JacobiRotationProblem prob = JacobiRotationProblem(mymat);
  prob.solve(1e-8);
  prob.printResultMatrix();
}
