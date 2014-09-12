#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <sstream>
#include <armadillo>

#include "ThomasAlgorithm.hpp"
#include "GaussianElimination.hpp"
#include "relativeError.hpp"

using namespace std;
using namespace arma;

void exerciseB(int);
void exerciseC(double,double,int);
void exerciseD(int);

void errorAnalysis(double*,int*,int);
double f(double);
double u(double);
double maxRelativeError(int,double,double*);
double duration(clock_t,clock_t);

int main(int argc, char* argv[])
{
  // Usage
  string usage = "Usage: ./<exe> -<exercise> -n <n> -logA <logA> -logB<logB>";
  if (argc == 1) {
    cout << usage << endl;
    return 1;
  }

  string exe = "none";
  int n = 0;
  double logA = 0;
  double logB = 0;

  // Traverse commandline
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i],"-b") == 0) {
      exe = "b";
    } else if (strcmp(argv[i],"-c") == 0) {
      exe = "c";
    } else if (strcmp(argv[i],"-d") == 0) {
      exe = "d";
    } else if (strcmp(argv[i],"-n") == 0) {
      n = atoi(argv[i+1]);
      i++;
    } else if (strcmp(argv[i],"-logA") == 0) {
      logA = atof(argv[i+1]);
      i++;
    } else if (strcmp(argv[i],"-logB") == 0) {
      logB = atof(argv[i+1]);
      i++;
    }
  }

  // Run correct exercise
  if (n == 0) {
    cout << usage << endl;
    return 1;
  }

  if (exe == "b") {
    exerciseB(n);
  } else if (exe == "c") {
    if (logA == 0 || logB == 0) {
      cout << usage << endl;
      return 1;
    }

    exerciseC(logA,logB,n);
  } else if (exe == "d") {
    exerciseD(n);
  }

  return 0;
}

// ***********************************************

/*
  The source term
*/
double f(double x) {
  return 100*exp(-10*x);
}

/*
 * The analytical solution. The closed form.
 */
double u(double x) {
  return 1 - (1 - exp(-10))*x - exp(-10*x);
}

/*
 * Returns the duration between two clock_t times (cpu-time) in microseconds
 */
double duration(clock_t start, clock_t finish) {
  return (finish-start)/((double)CLOCKS_PER_SEC/1000000.0);
}

/*
 * Finds the maximum relative error from two sets of solutions
 *
 * n: Length of num array. Must be n long.
 * h: Step size
 * num: The numerical solution
 */
double maxRelativeError(int n, double h, double* num) {
  double *errs = new double[n];
  double *ana = new double[n];
  for (int i = 0; i < n; i++) { ana[i] = u(i*h); }

  relativeError(errs,num,ana,n);

  delete [] ana;

  // Find max value from epsilon array
  double max = -10;
  for (int i = 1; i < n-1; i++) {
    if (errs[i] > max) {
      max = errs[i];
    }
  }

  delete [] errs;
  return max;
}

/*
 * Finds the relative error for a given spectrum of hs. The two arrays given
 * must be of same length.
 *
 * - eps: Empty array for error
 * - ns: Array of the number of steps
 * - n: Length of array h, how many steps to check for?
 */
void errorAnalysis(double* eps, int* ns, int n_errors) {
  double c1,h;
  double *d,*v;

  for (int i = 0; i < n_errors; i++) {
    // Solve the problem for given n

    // Ugly hack: Assume integration area is 1
    h = 1.0 / (ns[i] + 1);
    c1 = -h*h;

    d = new double[ns[i]+1]; // d right hand side coeffs
    v = new double[ns[i]+2]; // numerical solution

    for (int l = 0; l < ns[i]+2; l++) { v[l] = 0; }
    for (int l = 0; l < ns[i]+1; l++) { d[l] = c1*f(l*h); }

    GESolve(v,d,ns[i]);

    // The max error for this solution
    eps[i] = maxRelativeError(ns[i], h, v);

    delete [] d;
    delete [] v;

    cout << "Progress: i = " << i << " of " << n_errors << "\r" << flush;
  }
}

/*
  Exercise b
*/
void exerciseB(int n) {
  clock_t start,finish;

  double a = 0.0; double b = 1.0;
  double h = (b - a) / (n+1);
  double c1 = -h*h;
  double *v = new double[n+2];

  // Discrete grid xes
  double *x = new double[n+2];
  for (int i = 0; i < n+2; i++) { x[i] = i*h; v[i] = 0; }

  double *d = new double[n+1];
  for (int i = 0; i < n+2; i++) { d[i] = c1*f(x[i]); }

  start = clock();

  // Integration
  GESolve(v,d,n);

  finish = clock();

  // Write result to file
  char dat_time[80];
  time_t timer;
  struct tm * timeinfo;
  time(&timer);
  timeinfo = localtime(&timer);
  strftime(dat_time,80,"%a_%b_%d_%H_%M_%S.dat",timeinfo);

  ostringstream oss;
  oss << "data_n" << n << "_" << dat_time;
  string filename = oss.str();
  ofstream ofile;
  ofile.open(filename.c_str());

  // Print loop
  for (int i = 0; i < n+2; i++) {
    ofile << x[i] << " " << v[i] << endl;
  }

  ofile.close();

  cout << "Completed in: " << duration(start,finish) << " microseconds." << endl;

  // Delete dynamically allocated arrays
  delete [] d;
  delete [] x;
  delete [] v;
}

/*
 * Exercise c
 *
 * Performs an error analysis with ns from a to b.
 */
void exerciseC(double logA, double logB, int steps) {
  double increment = (logB - logA) / steps;
  double logn = logA;
  int *n_arr = new int[steps];
  double *eps_arr = new double[steps];
  double *logh_arr = new double[steps];

  int step_c = 0;
  while (logn < logB) {
    n_arr[step_c] = (int) pow(10,logn);
    logn += increment;
    logh_arr[step_c] = log10(1.0 / (n_arr[step_c] + 1));
    step_c++;
  }

  // Perform the analysis
  errorAnalysis(eps_arr,n_arr,steps);

  // Write eps and h to a file
  char dat_time[80];
  time_t timer;
  struct tm * timeinfo;
  time(&timer);
  timeinfo = localtime(&timer);
  strftime(dat_time,80,"%a_%b_%d_%H_%M_%S.dat",timeinfo);

  ostringstream oss;
  oss << "error_analysis_logA" << logA << "_logB" << logB << "_" << dat_time;
  string filename = oss.str();
  ofstream ofile;
  ofile.open(filename.c_str());

  // Print loop
  for (int i = 0; i < steps; i++) {
    ofile << logh_arr[i] << " " << eps_arr[i] << endl;
  }

  ofile.close();

  delete [] n_arr;
  delete [] eps_arr;
  delete [] logh_arr;
}

/*
 * Exercise d
 *
 * Comparing effectiveness with LU decomposition
 */
void exerciseD(int n) {
  clock_t start,finish;

  mat A = zeros(n,n);
  A(0,0) = 2; A(0,1) = -1;
  A(n-1,n-1) = 2; A(n-1,n-2) = -1;
  for (int i = 1; i < n-1; i++) {
    A(i,i) = 2;
    A(i,i-1) = -1;
    A(i,i+1) = -1;
  }

  double a = 0.0; double b = 1.0;
  double h = (b - a) / (n+1);
  double c1 = h*h;

  mat L,U;
  vec d = zeros(n);
  vec z = zeros(n);
  vec x = zeros(n);
  for (int i = 0; i < n; i++) { d(i) = c1*f(i*h); }

  start = clock();

  lu(L,U,A);
  solve(z,L,d);
  solve(x,U,z);

  finish = clock();

  vec solution = zeros(n+2);
  solution(0) = 0; solution(n+1) = 0;
  for (int i = 1; i < n+1; i++) {
    solution(i) = x(i-1);
  }

  // Write solution to file
  char dat_time[80];
  time_t timer;
  struct tm * timeinfo;
  time(&timer);
  timeinfo = localtime(&timer);
  strftime(dat_time,80,"%a_%b_%d_%H_%M_%S.dat",timeinfo);

  ostringstream oss;
  oss << "lu_decomp_sol_n" << n << "_" << dat_time;
  string filename = oss.str();
  ofstream ofile;
  ofile.open(filename.c_str());

  // Print loop
  for (int i = 0; i < n+2; i++) {
    ofile << i*h << " " << solution(i) << endl;
  }

  ofile.close();

  cout << "Completed in: " << duration(start,finish) << " microseconds." << endl;
}
