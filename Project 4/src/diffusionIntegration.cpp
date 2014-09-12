#include "diffusionIntegration.hpp"

using namespace arma;
using namespace std;

/*
 * @param dt Time step.
 * @param T Runtime.
 * @param dx Step size.
 * @param ic Initial condition.
 * @param bc Boundary conditions.
 * @param filepath Path to savefile.
 */
void forwardEulerDiffusion(double dt, double T, double dx, colvec ic, colvec bc, string filepath) {
  // Spatial parameters
  double a = 0;
  double b = 1;
  int nx = b / dx;
  colvec x = linspace<colvec>(a, b, nx);

  // Time parameters
  double t = 0;

  // Physical constants
  double D = 1;

  // Create central difference matrix
  mat ddx2 = zeros(x.size(), x.size());
  ddx2.diag(0).fill(-2);
  ddx2.diag(1).fill(1);
  ddx2.diag(-1).fill(1);

  // Check on alpha status
  double alpha = dt / (dx*dx);
  if (alpha > 0.5) {
    cout << "Bad choice of steps sizes." << endl;
  }

  // Open savefile
  ofstream outfile;
  outfile.open(filepath.c_str());

  // Write header
  outfile << "#HEADER#" << endl;
  outfile << "Method: ForwardEuler" << endl;
  outfile << "T: " << T << endl << "dt: " << dt << endl;
  outfile << "dx: " << dx << endl;
  outfile << "x: ";
  for (int i = 0; i < nx; i++) { outfile << x(i) << " "; }
  outfile << endl;
  outfile << "#/HEADER#" << endl;

  // Initial state
  colvec u = D*alpha*ddx2*ic;

  // Add boundary conditions
  u(0) = bc(0);
  u(u.size()-1) = bc(1);

  // Save initial state
  for (int i = 0; i < nx; i++) { outfile << u(i) << " "; }
  outfile << endl;

  while (t < T) {
    // Find change
    u = u + alpha*ddx2*u;

    // Add boundary conditions
    u(0) = bc(0);
    u(u.size()-1) = bc(1);

    // Save current state to file
    for (int i = 0; i < nx; i++) { outfile << u(i) << " "; }
    outfile << endl;

    t += dt;
  }

  // Close savefile
  outfile.close();
}

/*
 * @param dt Time step.
 * @param T Runtime.
 * @param dx Step size.
 * @param ic Initial condition.
 * @param bc Boundary conditions.
 * @param filepath Path to savefile.
 */
void backwardEulerDiffusion(double dt, double T, double dx, colvec ic, colvec bc, string filepath) {
  // Spatial parameters
  double a = 0;
  double b = 1;
  int nx = (b / dx) + 1;
  colvec x = linspace<colvec>(a, b, nx);

  // Time parameters
  double t = 0;

  // Physical constants
  double D = 1;

  // Create central difference matrix
  mat ddx2 = zeros(x.size(), x.size());
  ddx2.diag(0).fill(-2);
  ddx2.diag(1).fill(1);
  ddx2.diag(-1).fill(1);

  // Check on alpha status
  double alpha = dt / (dx*dx);
  if (alpha > 0.5) {
    cout << "Bad choice of steps sizes." << endl;
  }

  // Open savefile
  ofstream outfile;
  outfile.open(filepath.c_str());

  // Write header
  outfile << "#HEADER#" << endl;
  outfile << "Method: BackwardEuler" << endl;
  outfile << "T: " << T << endl << "dt: " << dt << endl;
  outfile << "dx: " << dx << endl;
  outfile << "x: ";
  for (int i = 0; i < nx; i++) { outfile << x(i) << " "; }
  outfile << endl;
  outfile << "#/HEADER#" << endl;

  // Initial state
  colvec u = zeros(nx);
  u = D*alpha*ddx2*ic;

  // Add boundary conditions
  u(0) = bc(0);
  u(u.size()-1) = bc(1);

  // Save initial state
  for (int i = 0; i < nx; i++) { outfile << u(i) << " "; }
  outfile << endl;

  // Need identity matrix
  mat I = eye(nx,nx);

  mat A; colvec B; colvec sol;
  while (t < T) {
    // Find change by solving linear set of equations
    A = I - alpha*ddx2;
    B = u;
    sol = solve(A,B);
    u(span(1,nx-2)) = sol(span(1,nx-2));

    // Save current state to file
    for (int i = 0; i < nx; i++) { outfile << u(i) << " "; }
    outfile << endl;

    t += dt;
  }

  // Close savefile
  outfile.close();
}

/*
 * @param dt Time step.
 * @param T Runtime.
 * @param dx Step size.
 * @param ic Initial condition.
 * @param bc Boundary conditions.
 * @param filepath Path to savefile.
 */
void crankNicolsonDiffusion(double dt, double T, double dx, colvec ic, colvec bc, string filepath) {
  // Spatial parameters
  double a = 0;
  double b = 1;
  double nx = b / dx;
  colvec x = linspace<colvec>(a, b, nx);

  // Time parameters
  double t = 0;

  // Physical constants
  double D = 1;

  // Create central difference matrix
  mat ddx2 = zeros(x.size(), x.size());
  ddx2.diag(0).fill(-2);
  ddx2.diag(1).fill(1);
  ddx2.diag(-1).fill(1);

  double alpha = dt / (dx*dx);

  // Open savefile
  ofstream outfile;
  outfile.open(filepath.c_str());

  // Write header
  outfile << "#HEADER#" << endl;
  outfile << "Method: CrankNicolson" << endl;
  outfile << "T: " << T << endl << "dt: " << dt << endl;
  outfile << "dx: " << dx << endl;
  outfile << "x: ";
  for (int i = 0; i < nx; i++) { outfile << x(i) << " "; }
  outfile << endl;
  outfile << "#/HEADER#" << endl;

  // Initial state
  colvec u = D*alpha*ddx2*ic;

  // Add boundary conditions
  u(0) = bc(0);
  u(u.size()-1) = bc(1);

  // Save initial state
  for (int i = 0; i < nx; i++) { outfile << u(i) << " "; }
  outfile << endl;

  // Need identity matrix
  mat I = eye(nx,nx);

  mat A; colvec B;
  while (t < T) {
    // Solve equation to find change
    A = I - 0.5*alpha*ddx2;
    B = (I + 0.5*alpha*ddx2)*u;
    u = solve(A,B);

    // Add boundary conditions
    u(0) = bc(0);
    u(u.size()-1) = bc(1);

    // Save current state to file
    for (int i = 0; i < nx; i++) { outfile << u(i) << " "; }
    outfile << endl;

    t += dt;
  }

  // Close savefile
  outfile.close();
}
