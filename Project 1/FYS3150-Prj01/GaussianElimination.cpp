#include "GaussianElimination.hpp"

/*
 * n: Number of steps used in conjugtion with h
 */
void GESolve(double *x, double *d, int n) {
  // Create the three diagonals
  double *a,*b,*c;
  double factor;
  a = new double[n+1]; b = new double[n+1]; c = new double[n+1];
  for (int i = 0; i < n+1; i++) { a[i] = 1; b[i] = -2; c[i] = 1; }
  // Remember last and first element
  a[1] = 0; c[n] = 0;

  // Forward substitution
  for (int i = 2; i < n+2; i++) {
    factor = a[i]/b[i-1];
    b[i] -= factor*c[i-1];
    d[i] -= factor*d[i-1];
  }

  // Backward substitution
  x[n] = d[n]/b[n];
  for (int i = n-1; i > 0; i--) {
    x[i] = (d[i] - c[i]*x[i+1])/b[i];
  }

  // Delete dynamically allocated arrays
  delete [] a;
  delete [] b;
  delete [] c;
}
