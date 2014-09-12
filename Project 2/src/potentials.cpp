#include "potentials.hpp"

/*
 * Returns dimensionless potential for three-dimensional harmonic oscillator.
 *
 * @param rho Dimensionless radius.
 *
 * @return Potential as double.
 */
double harOscV(double rho) {
  return rho*rho;
}

/*
 * A dimensionless potential for three-dimensional harmonic oscillator with two
 * electrons including the coulomb interaction.
 *
 * @param rho Dimensionless radius.
 * @param omegaR Dimensionless parameter that reflects strength of oscillator
 * potential.
 *
 * @return Potential as double.
 */
double harOscTwoElectronsV(double rho, double omegaR) {
  return omegaR*omegaR*rho*rho + 1.0 / rho;
}
