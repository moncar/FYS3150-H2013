#include <iostream>
#include <armadillo>
#include <sstream>

#include "SolarSystem.hpp"

using namespace std;
using namespace arma;

/*
 * Constants
 */
const int DIMENSIONALITY = 2;
const string OBJECTS_DATA_PATH = "../data/objects";

/*
 * Constructor that takes datafile. An already setup system can be read from
 * this file. Assumes file is in a certain syntax. Examples can be found in
 * `data` directory.
 *
 * This function is dependent upon that the dimensionality is set to 2. This is
 * bad.
 * TODO: Refine and make reading independent of dimensions. Maybe take as int
 * arg.
 *
 * @param systemfile The file to read the system data from.
 */
SolarSystem :: SolarSystem(string systemfile) {
  // Values needed to temporary store information from file in
  ifstream datafile;
  datafile.open(systemfile.c_str());
  string line,col,id;
  bool record = false; // Switch for knowing when past meta part in systemfile
  double x0,y0,v0x,v0y,m; // For containing the values

  // Procedure for reading the columns line by line
  while (datafile.good()) {
    // Wait till past meta header
    if (!record) {
      datafile >> line;
      if (strcmp(line.c_str(),"#DATA#") == 0) { record = true; }
    } else {
      // Get the data
      datafile >> id >> x0 >> y0 >> v0x >> v0y >> m;

      // Check for EOF
      if (!datafile.good()) { break; }

      // Create temporary vectors for the data (in order to feed to functions)
      vec position; position << x0 << y0;
      vec velocity; velocity << v0x << v0y;

      // Create the dynamic file name
      ostringstream oss;
      oss << OBJECTS_DATA_PATH << "/" << id << ".dat";
      string savefile = oss.str();

      // Create the object and store it
      CelestialObject newObject = CelestialObject(id,position,velocity,m,savefile);
      addObject(newObject);
    }
  }

  // Close the input systemfile
  datafile.close();
}

/*
 * Advances the solar system by a given timestep via the RungeKutta4 scheme.
 *
 * @param dt The timestep to advance.
 */
void SolarSystem :: advance(double dt) {
  // Original velocity and position
  mat pn = zeros<mat>(DIMENSIONALITY,getNoOfObjects());

  /*
   * Matrices used to store the Ks for both acceleration and velocity for all
   * objects in this system.
   */
  mat aK1 = zeros<mat>(DIMENSIONALITY,getNoOfObjects());
  mat aK2 = zeros<mat>(DIMENSIONALITY,getNoOfObjects());
  mat aK3 = zeros<mat>(DIMENSIONALITY,getNoOfObjects());
  mat aK4 = zeros<mat>(DIMENSIONALITY,getNoOfObjects());

  mat vK1 = zeros<mat>(DIMENSIONALITY,getNoOfObjects());
  mat vK2 = zeros<mat>(DIMENSIONALITY,getNoOfObjects());
  mat vK3 = zeros<mat>(DIMENSIONALITY,getNoOfObjects());
  mat vK4 = zeros<mat>(DIMENSIONALITY,getNoOfObjects());

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
  
  // Save positions to files
  saveAllPositions();

  // After advance
  t += dt;
}

/*
 * Closes all datafiles.
 */
void SolarSystem :: close() {
  for (int i = 0; i < getNoOfObjects(); i++) {
    objects[i].close();
  }
}

/*
 * Adds a new celestial object to the collection of objects.
 *
 * @param newObject The new object. Instance of `CelestialObject`.
 */
void SolarSystem :: addObject(CelestialObject newObject) {
  objects.push_back(newObject);
}

/*
 * Runs through list of celestial objects and makes all save their current
 * coordinate to their own datafile.
 */
void SolarSystem :: saveAllPositions() {
  for (int i = 0; i < getNoOfObjects(); i++) {
    objects[i].saveCurrentPos();
  }
}

/*
 * Get functions
 */

/*
 * @return All velocities as arma matrix.
 */
mat SolarSystem :: getAllVelocities() {
  mat v = zeros<mat>(DIMENSIONALITY,getNoOfObjects());

  for (int i = 0; i < getNoOfObjects(); i++) {
    v.col(i) = objects[i].getV();
  }

  return v;
}

/*
 * @return All positions as arma matrix.
 */
mat SolarSystem :: getAllPositions() {
  mat p = zeros<mat>(DIMENSIONALITY,getNoOfObjects());

  for (int i = 0; i < getNoOfObjects(); i++) {
    p.col(i) = objects[i].getPos();
  }

  return p;
}

/*
 * @return All current accelerations as arma matrix.
 */
mat SolarSystem :: getAllAccelerations() {
  mat a = zeros<mat>(DIMENSIONALITY,getNoOfObjects());

  for (int i = 0; i < getNoOfObjects(); i++) {
    a.col(i) = getForces(objects[i]) / objects[i].getM();
  }

  return a;
}

/*
 * Set methods
 */

/*
 * Set velocity of all objects
 *
 * @param velocities Arma matrix containing all velocities as columns.
 */
void SolarSystem :: setAllVelocities(mat velocities) {
  for (int i = 0; i < getNoOfObjects(); i++) {
    objects[i].setV(velocities.col(i));
  }
}

/*
 * Set position of all objects
 *
 * @param positions Arma matrix containing all positions as columns.
 */
void SolarSystem :: setAllPositions(mat positions) {
  for (int i = 0; i < getNoOfObjects(); i++) {
    objects[i].setPos(positions.col(i));
  }
}

/*
 * Returns a double representing the total force working on this object from
 * all other objects in the system.
 */
vec SolarSystem :: getForces(CelestialObject object) {
  vec force = zeros(DIMENSIONALITY); vec r;

  for (int i = 0; i < getNoOfObjects(); i++) {
    // Not find force from itself
    if (strcmp(objects[i].getId().c_str(),object.getId().c_str()) == 0) { continue; }

    // Gravitational force
    r = object.getDistTo(objects[i]);
    force += r * objects[i].getM()*object.getM() / pow(norm(r,2),3);
  }

  return force;
}

/*
 * @return The number of objects in the system.
 */
int SolarSystem :: getNoOfObjects() {
  return objects.size();
}

/*
 * @return The current time the system is in.
 */
double SolarSystem :: getT() {
  return t;
}
