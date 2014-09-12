#include <armadillo>
#include <iostream>
#include <ctime>

#include "System.hpp"
#include "constants.cpp"

using namespace arma;
using namespace std;

/*
 * @param N The number of random particles to contain in the system at startup.
 * @param dimension The dimension to run the simulation for. Usually 3.
 * @param R0 Size of cluster.
 * @param meanMass The mean mass of the cluster.
 * @param t0 The time when which the system starts. Usually 0.
 *
 * TODO: Get rid of all cluster related things in this (should be) generic
 * class. Create subclass Cluster.
 */
System :: System(int N, int dimension, double R0, double meanMass, double t0) {
  this->N = 0;
  this->dimension = dimension;
  this->t = t0;
  this->R0 = R0;
  this->meanMass = meanMass;

  this->rho0 = ( 3.0 * N * meanMass ) / ( 4 * PI * R0*R0*R0 );
  this->G = ( 3.0 * PI ) / ( 32.0 * rho0 );
  this->epsilon = 0;

  this->saveEach = 1;
  this->advanceCounter = 0; // For counting number of advances

  time_t seed = time(0);
  int seed1 = - seed;
  int seed2 = seed1 * 0.3;

  this->massGenerator = new Random(seed1);
  this->positionGenerator = new Random(seed2);

  // Initialize the system
  createInitialSystem(N);
}

/*
 * Creates the all particles to be contained in the system and adds them to
 * objects list.
 */
void System :: createInitialSystem(int N) {
  for (int i = 0; i < N; i++) {
    Particle newParticle = makeNewParticle(i);
    addParticle(newParticle);
  }
}

/*
 * Turn on the smoothing grav force.
 *
 * @param epsilon The value to use in smoothing function.
 */
void System :: setSmoothing(double epsilon) {
  this->epsilon = epsilon;
}

/*
 * Optional can save not every timestep. Good for small timesteps. Not that
 * much data.
 *
 * @param saveEach Int describing only to save each multiple of that int. (1 is
 * all, 10 is each tenth).
 */
void System :: setSaveEach(int saveEach) {
  this->saveEach = saveEach;
}

/*
 * Adds a particle to the list of objects.
 *
 * TODO: Error check for if particle already present.
 */
void System :: addParticle(Particle newParticle) {
  objects.push_back(newParticle);
  N++;
}

/*
 * @number An integer describing which number in creation this particle has.
 * Used as identifier for particle.
 *
 * @return A new particle.
 */
Particle System :: makeNewParticle(int number) {
  Particle newParticle = Particle( makeNewInitialPosition(),
                                   makeNewInitialVelocity(),
                                   makeNewMass(),
                                   number );
  return newParticle;
}

/*
 * @return A new mass for a particle from gaussian distribution using correct
 * units.
 */
double System :: makeNewMass() {
  return massGenerator->normal(meanMass,1);
}

/*
 * @return A vector containing a new initial position chosen from correct
 * gaussian distribution.
 */
colvec System :: makeNewInitialPosition() {
  return positionGenerator->randomSphericalPosition(R0);
}

/*
 * @return A vector containing a new initial velocity chosen from correct
 * gaussian distribution.
 */
colvec System :: makeNewInitialVelocity() {
  return zeros(dimension);
}

/* DATA STORAGE */

/*
 * Open file for saving.
 *
 * @param savefile Name of file to save to.
 */
void System :: open(string savefile) {
  datafile = new ofstream(savefile.c_str());

  // Write header with meta-information
  *datafile << "*********META********" << endl;
  *datafile << "R0=" << R0 << endl;
  *datafile << "N=" << N << endl;
  *datafile << "epsilon=" << epsilon << endl;
  *datafile << "saveEach=" << saveEach << endl;
  *datafile << "********/META********" << endl;
}

/*
 * Close datafile.
 */
void System :: close() {
  datafile->close();
}

/*
 * Saves the current state of the system (all positions) to a file using the
 * following syntax. For each call, a new line is added which includes a
 * parantheses with all information for that object.
 *
 * t=<current time> E=<tot eng> (x y z bound=False/True) (x y z bound=False/True) (x y z bound=False/True) ... for all particles
 *
 * Also general for all dimensions.
 */
void System :: saveState() {
  colvec position;
  double kineticEnergy,potentialEnergy,potentialEnergyBound;

  // Write current time
  *datafile << "t=" << t << " ";
  // Write current total energy
  *datafile << "E=" << getTotalEnergy(false) << " ";

  for (int i = 0; i < N; i++) {
    // Position
    position = objects[i].getPos();
    *datafile << "(x=" << position[0] << " y=" << position[1] << " z=" << position[2];

    // Potential and kinetic energy
    kineticEnergy = objects[i].getKineticEnergy();
    potentialEnergy = objects[i].getPotentialEnergy(false);
    potentialEnergyBound = objects[i].getPotentialEnergy(true);
    *datafile << " kinetic=" << kineticEnergy;
    *datafile << " potential=" << potentialEnergy;
    *datafile << " potential_bound=" << potentialEnergyBound;

    // Is particle bound ?
    *datafile << " bound=";
    if (isBound(objects[i])) {
      *datafile << "True";
    } else {
      *datafile << "False";
    }

    *datafile << ")";
  }

  *datafile << endl;
}

void System :: updateEnergies() {
  for (int i = 0; i < N; i++) {
    objects[i].setKineticEnergy(calculateKineticEnergy(objects[i]));
    objects[i].setPotentialEnergy(calculatePotentialEnergy(objects[i],false),
                                  calculatePotentialEnergy(objects[i],true));
  }
}

/* GET METHODS */

/*
 * @paramo object Object to check if is bound.
 *
 * @return Bool indicating whether or not given particle is bound to the
 * system.
 */
bool System :: isBound(Particle object) {
  if ( object.getKineticEnergy() + object.getPotentialEnergy(false) < 0) {
    return true;
  }
  return false;
}

/*
 * @param bound Bool indicating wether or not only count bound objects in
 * potential.
 *
 * @return The total energy of the system in this state.
 */
double System :: getTotalEnergy(bool bound) {
  double kinetic,potential;

  for (int i = 0; i < N; i++) {
    kinetic += objects[i].getKineticEnergy();
    potential += objects[i].getPotentialEnergy(bound);
  }

  return kinetic + 0.5*potential;
}

/*
 * @return The total kinetic energy of a given particle.
 */
double System :: calculateKineticEnergy(Particle object) {
  return 0.5*object.getMass()*pow(norm(object.getVel(),2),2);
}

/*
 * @param object The object to find potential energy on.
 * @param bound Switch to only count from bound bodies
 *
 * @return The potential energy from all other particles on the given particle.
 */
double System :: calculatePotentialEnergy(Particle object, bool bound) {
  if (bound && !isBound(object)) {
    return 0;
  }

  double potential = 0;
  colvec r;

  for (int i = 0; i < N; i++) {
    if (objects[i].getId() == object.getId()) {
      continue;
    }
    if (bound && !isBound(objects[i])) {
      continue;
    }

    r = object.getDistTo(objects[i]);
    potential += - G * objects[i].getMass()*object.getMass() / norm(r,2);
  }

  return potential;
}

/*
 * @return All velocities as arma matrix.
 */
mat System :: getAllVelocities() {
  mat v = zeros<mat>(dimension,N);

  for (int i = 0; i < N; i++) {
    v.col(i) = objects[i].getVel();
  }

  return v;
}

/*
 * @return All positions as arma matrix.
 */
mat System :: getAllPositions() {
  mat p = zeros<mat>(dimension,N);

  for (int i = 0; i < N; i++) {
    p.col(i) = objects[i].getPos();
  }

  return p;
}

/*
 * @return All current accelerations as arma matrix.
 */
mat System :: getAllAccelerations() {
  mat a = zeros<mat>(dimension,N);

  for (int i = 0; i < N; i++) {
    a.col(i) = getForces(objects[i]) / objects[i].getMass();
  }

  return a;
}

/*
 * @return A vector representing the total force working on this object from
 * all other objects in the system.
 */
colvec System :: getForces(Particle object) {
  vec force = zeros(dimension); vec r;

  for (int i = 0; i < N; i++) {
    if (object.getId() == objects[i].getId()) {
      continue;
    }

    // Gravitational force
    r = object.getDistTo(objects[i]);
    force += G * r * objects[i].getMass()*object.getMass() / ((pow(norm(r,2),2) + pow(epsilon,2))*norm(r,2));
  }

  return force;
}

/*
 * @return Number of objects currently in system.
 */
int System :: getNoOfObjects() {
  return N;
}

/* SET METHODS */

/*
 * Set velocity of all objects
 *
 * @param velocities Arma matrix containing all velocities as columns.
 */
void System :: setAllVelocities(mat velocities) {
  for (int i = 0; i < N; i++) {
    objects[i].setVel(velocities.col(i));
  }
}

/*
 * Set position of all objects
 *
 * @param positions Arma matrix containing all positions as columns.
 */
void System :: setAllPositions(mat positions) {
  for (int i = 0; i < N; i++) {
    objects[i].setPos(positions.col(i));
  }
}
