#include "Particle.hpp"

using namespace arma;
using namespace std;

/*
 * @param pos Vector containing position.
 * @param vel Vector containing velocity.
 * @param mass Double describing mass.
 * @param id Id identifying particle (when created).
 */
Particle :: Particle(vec pos, vec vel, double mass, int id) {
  this->m = mass;
  this->pos = pos;
  this->vel = vel;
  this->id = id;
}

/*
 * @return Distance to another particle.
 */
vec Particle :: getDistTo(Particle other) {
  vec dist = other.getPos() - getPos();
  return dist;
}

/* GET METHODS */

double Particle :: getMass() {
  return m;
}

vec Particle :: getPos() {
  return pos;
}

vec Particle :: getVel() {
  return vel;
}

int Particle :: getId() {
  return id;
}

/*
 * @param bound Indicating whether or not return energy only from bound
 * particles.
 */
double Particle :: getPotentialEnergy(bool bound) {
  if (bound) {
    return potentialEnergyBound;
  }
  return potentialEnergy;
}

double Particle :: getKineticEnergy() {
  return kineticEnergy;
}

/* SET METHODS */

/*
 * @param potentialEnergy New potential energy.
 * @param potentialEnergy New potential energy from bound objects.
 */
void Particle :: setPotentialEnergy(double potentialEnergy, double potentialEnergyBound) {
  this->potentialEnergy = potentialEnergy;
  this->potentialEnergyBound = potentialEnergyBound;
}

/*
 * @param kineticEnergy New kinetic energy.
 */
void Particle :: setKineticEnergy(double kineticEnergy) {
  this->kineticEnergy = kineticEnergy;
}

/*
 * Set new positions.
 *
 * @param newPos New position of object.
 */
void Particle :: setPos(vec newPos) {
  pos = newPos;
}

/*
 * Updates velocity of object.
 *
 * @param newVel New velocity of object.
 */
void Particle :: setVel(vec newVel) {
  vel = newVel;
}
