#include <armadillo>

#include "gaussiandeviate.hpp"

/*
 * Motivation for writing this container class for a random distribution
 *
 * --> We want to make it easier to choose standard deviation and mean value
 * --> We don't need to worry about seeds interfering with each other, just use
 *  a different one for each random generator.
 * --> Better compact code
 */
class Random {
  public:
    Random(int);

    double normal(double mean, double deviation);
    double normalStandard();
    double uniform(double a, double b);
    double uniformStandard();

    arma::colvec randomSphericalPosition(double R);

  private:
    long seed;

    double deviation,mean;
};
