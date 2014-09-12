#include "System.hpp"

using namespace arma;

/*
 * Advances system using Runge-Kutta 4 scheme.
 *
 * @param dt The time step to step forwards.
 */
void SystemRK4 :: advance(double dt) {
  // Original velocity and position
  mat pn = zeros<mat>(dimension,N);

  /*
   * Matrices used to store the Ks for both acceleration and velocity for all
   * objects in this system.
   */
  mat aK1 = zeros<mat>(dimension,N);
  mat aK2 = zeros<mat>(dimension,N);
  mat aK3 = zeros<mat>(dimension,N);
  mat aK4 = zeros<mat>(dimension,N);

  mat vK1 = zeros<mat>(dimension,N);
  mat vK2 = zeros<mat>(dimension,N);
  mat vK3 = zeros<mat>(dimension,N);
  mat vK4 = zeros<mat>(dimension,N);

  // Store first position for last step
  pn = getAllPositions();

  // K1
  aK1 = getAllAccelerations();
  vK1 = getAllVelocities();

  // Update system
  setAllVelocities( vK1 + 0.5*dt*aK1 );
  setAllPositions( pn + 0.5*dt*vK1 );

  // K2
  aK2 = getAllAccelerations();
  vK2 = getAllVelocities();

  // Update system
  setAllVelocities( vK1 + 0.5*dt*aK2 );
  setAllPositions( pn + 0.5*dt*vK2 );

  // K3
  aK3 = getAllAccelerations();
  vK3 = getAllVelocities();

  // Update system
  setAllVelocities( vK1 + dt*aK3 );
  setAllPositions( pn + dt*vK3 );

  // K4
  aK4 = getAllAccelerations();
  vK4 = getAllVelocities();

  // The real step
  setAllVelocities( vK1 + (1./6) * dt * (aK1 + 2*aK2 + 2*aK3 + aK4) );
  setAllPositions( pn + (1./6) * dt * (vK1 + 2*vK2 + 2*vK3 + vK4) );
  
  // After advance
  t += dt;

  // Update energies
  updateEnergies();

  // Update number of advances
  advanceCounter++;

  if (advanceCounter == saveEach) {
    // Save to file
    saveState();
    advanceCounter = 0;
  }
}
