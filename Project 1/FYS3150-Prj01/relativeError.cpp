#include "relativeError.hpp"
#include <cmath>

/*
 * Calculates the relative error between the two given arrays. u is used as the
 * relative term.
 *
 * - v,u: Arrays containing the two solutions. Must be n long.
 * - eps: Array that will be filled with the errors. Must be n long.
 * - n: Length of v,u,eps
 */
void relativeError(double* error, double* numerical, double* analytical, int len) {
  // Make sure epsilon array is zeros
  for (int i = 0; i < len; i++) { error[i] = 0; }

  /*
   * Avoid the first and last element since they are the initial conditions and
   * the error is
   */
  for (int i = 1; i < len-1; i++) {
    error[i] = log10( fabs( (numerical[i]-analytical[i])/analytical[i] ) );
  }
}
