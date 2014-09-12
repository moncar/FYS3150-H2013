#include "ThomasAlgorithm.hpp"

/*
  Function solves a linear set of equations
  for a tridiagonal matrix (-1, 2, -1).

  Arguments:
  - *x pointer to array of x'es
  - *d pointer to array of d's
  - n number of elements
*/
void TASolve(double *x, double *d, int n) {
  // Create the three diagonals
  double *a,*b,*c;
  a = new double[n]; b = new double[n]; c = new double[n];
  for (int i = 0; i < n; i++) { a[i] = 1; b[i] = -2; c[i] = 1; }
  // Remember last and first element
  a[0] = 0; c[n-1] = 0;

  /*
    The solving algorithm
  */

  // Special first elements
  c[0] /= b[0];
  d[0] /= b[0];

  // Forward sweep
  for (int i = 1; i < n; i++) {
    c[i] /= (b[i] - a[i]*c[i-1]);
    d[i] = (d[i] - d[i-1]*a[i]) / (b[i] - a[i]*c[i-1]);
  }

  // Backward sweep
  x[n-1] = d[n-1];
  for (int i = n-1; i > 0; i--) {
    x[i] = d[i] - c[i]*x[i+1];
  }

  // Delete dynamically allocated arrays
  delete [] a;
  delete [] b;
  delete [] c;
}
