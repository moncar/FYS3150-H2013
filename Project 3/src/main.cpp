#include <iostream>
#include <armadillo>

#include "SolarSystem.hpp"

using namespace std;
using namespace arma;

/*
 * Main method takes a systemfile containing the setup of the system being
 * simulated, a dt and a T (how long the simulation is to run). All times are
 * in years.
 *
 * Usage:
 *   ./<exe> <dt> <T> <systemfile>
 */
int main(int argc, char* argv[]) {
  if (argc != 4) {
    cout << "Need step size and number of steps pluss file to fetch system from." << endl;
    return 1;
  }

  double dt = atof(argv[1]);
  double T = atof(argv[2]);
  double t = 0;

  string infile = argv[3];

  SolarSystem mySystem = SolarSystem(infile);

  while (t < T) {
    mySystem.advance(dt);
    t += dt;
  }

  mySystem.close();

  cout << "Finished simulation." << endl;
  return 0;
}
