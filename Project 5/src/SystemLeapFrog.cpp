#include "System.hpp"

using namespace arma;

/*
 * Advances the system using the leap frog scheme.
 *
 * @param dt The time step.
 */
void SystemLeapFrog :: advance(double dt) {
  mat v_half = getAllVelocities() + 0.5*dt*getAllAccelerations();
  
  setAllPositions( getAllPositions() + dt*v_half );
  setAllVelocities( v_half + 0.5*dt*getAllAccelerations() );

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
