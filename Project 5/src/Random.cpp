#include <cmath>

#include "Random.hpp"
#include "constants.cpp"

using namespace arma;

/*
 * Initiates a new distribution.
 *
 * @param seed Seed to start at.
 */
Random :: Random(int seed) {
  this->seed = seed;
}

/*
 * @return A new random variable from distribution.
 */
double Random :: normal(double mean, double deviation) {
  return mean + deviation*gaussian_deviate(&seed);
}

/*
 * Relation found via integrating spherical coordinates under the assumption
 * that the volume element should be the same for all coordinate systems.
 *
 * @param R Radius of sphere to contain particles in.
 * TODO: Have as class constant?
 *
 * @return A vector representing a normally distributed point in a three
 * dimensional sphere. Note that the vector returned contains kartesian
 * coordinates.
 */
colvec Random :: randomSphericalPosition(double R) {
  double u = uniformStandard();
  double v = uniformStandard();
  double w = uniformStandard();

  double phi = 2.*PI*w;
  double theta = acos(1-2*v);
  double r = R*pow(u,0.333);

  colvec a = zeros(3);
  a[0] = r*sin(theta)*cos(phi);     // x
  a[1] = r*sin(theta)*sin(phi);     // y
  a[2] = r*cos(theta);              // z

  return a;
}

/*
 * @return A number from standard normal distribution (between 0 and 1).
 */
double Random :: uniformStandard() {
  return ran2(&seed);
}

/*
 * @return A random number from uniform distribution between a and b.
 */
double Random :: uniform(double a, double b) {
  return 0.1;
}

/*
 * @return A random number from standard normal distribution (mean 0 and
 * deviation 1).
 */
double Random :: normalStandard() {
  return gaussian_deviate(&seed);
}
