#include <vector>
#include <armadillo>

#include "CelestialObject.hpp"

class SolarSystem {
  public:
    SolarSystem(std::string);

    void addObject(CelestialObject);
    void advance(double);
    void close();

    int getNoOfObjects();
    double getT();

  private:
    double t;

    std::vector<CelestialObject> objects;

    arma::vec getForces(CelestialObject);

    arma::mat getAllVelocities();
    arma::mat getAllPositions();
    arma::mat getAllAccelerations();

    void setAllVelocities(arma::mat);
    void setAllPositions(arma::mat);
    void saveAllPositions();
};
