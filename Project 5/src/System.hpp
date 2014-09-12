#include <vector>
#include <armadillo>
#include <fstream>

#include "Particle.hpp"
#include "Random.hpp"

class System {
  public:
    /* METHODS */
    System(int N, int dimension, double R0, double meanMass, double t0);

    virtual void advance(double dt) { return; }

    void setSmoothing(double epsilon);
    void open(std::string savefilename);
    void close();
    void saveState();
    void setSaveEach(int);

    int getNoOfObjects();

  protected:
    /* METHODS */
    arma::mat getAllVelocities();
    arma::mat getAllPositions();
    arma::mat getAllAccelerations();

    double getTotalEnergy(bool);

    void setAllVelocities(arma::mat);
    void setAllPositions(arma::mat);
    void updateEnergies();

    /* FIELDS */
    int N, dimension, saveEach, advanceCounter;
    double t, rho0, G, meanMass, epsilon;

  private:
    /* METHODS */
    double makeNewMass();

    Particle makeNewParticle(int);

    arma::colvec getForces(Particle);
    arma::colvec makeNewInitialPosition();
    arma::colvec makeNewInitialVelocity();

    bool isBound(Particle);

    double calculateKineticEnergy(Particle);
    double calculatePotentialEnergy(Particle,bool);

    void createInitialSystem(int N);
    void addParticle(Particle);

    /* FIELDS */
    double R0;

    Random *massGenerator,*positionGenerator;

    std::vector<Particle> objects;
    std::ofstream* datafile;
};

class SystemRK4: public System {
  public:
    SystemRK4(int N, int dimension, double R0, double meanMass, double t0)
    : System(N,dimension,R0,meanMass,t0) { };

    void advance(double dt);
};

class SystemLeapFrog: public System {
  public:
    SystemLeapFrog(int N, int dimension, double R0, double meanMass, double t0)
    : System(N,dimension,R0,meanMass,t0) { };

    void advance(double dt);
};
