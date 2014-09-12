#include <fstream>
#include <iostream>

#include "CelestialObject.hpp"

using namespace std;
using namespace arma;

/*
 * @param id Short string describing object.
 * @param coors Initial coordinates of celestial object.
 * @param vel Initial velocity of celestial object.
 * @param m Mass of object.
 * @param savefile Path to file where object can store points as calculation is
 * done.
 */
CelestialObject :: CelestialObject(string id, vec coors, vec vel, double m, string savefile) {
  this->m = m;
  this->pos = coors;
  this->v = vel;
  this->id = id;

  /*
   * Open and store datafile. Header is written at once.
   */
  datafile = new ofstream(savefile.c_str());
  *datafile << "Positions for: " << getId() << endl << "#syntax: x y" << endl;

  // Save first point
  saveCurrentPos();
}

/*
 * Close datafile
 */
void CelestialObject :: close() {
  datafile->close();
}

/*
 * Set functions
 */

/*
 * Saves the current force exerted on this object.
 *
 * @param newF The force exerted.
 */
void CelestialObject :: setF(vec newF) {
  F = newF;
}

/*
 * Set new positions.
 *
 * @param newPos New position of object.
 */
void CelestialObject :: setPos(vec newPos) {
  pos = newPos;
}

/*
 * Updates velocity of object.
 *
 * @param newVel New velocity of object.
 */
void CelestialObject :: setV(vec newVel) {
  v = newVel;
}

/*
 * Stores current position to its given savefile.
 */
void CelestialObject :: saveCurrentPos() {
  vec currentPosition = getPos();
  for (int i = 0; i < currentPosition.n_elem; i++) {
    *datafile << currentPosition[i] << " ";
  }
  *datafile << endl;
}

/*
 * Get functions
 */

/*
 * Returns the distance from this object to a different one.
 *
 * @param other The other celestial object to find distance to.
 */
arma::vec CelestialObject :: getDistTo(CelestialObject other) {
  vec dist = other.getPos() - getPos();
  return dist;
}

double CelestialObject :: getM() { return m; }
vec CelestialObject :: getPos() { return pos; }
vec CelestialObject :: getV() { return v; }
vec CelestialObject :: getForce() { return F; }
string CelestialObject :: getId() { return id; }
