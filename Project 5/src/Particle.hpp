#include <armadillo>

class Particle {
  public:
    /* METHODS */
    Particle(arma::vec pos,arma::vec vel,double mass, int id);

    arma::vec getDistTo(Particle);
    arma::vec getPos(),getVel();

    int getId();

    double getMass();
    double getKineticEnergy();
    double getPotentialEnergy(bool);

    void setPos(arma::vec);
    void setVel(arma::vec);
    void setKineticEnergy(double);
    void setPotentialEnergy(double,double);

  private:
    /* FIELDS */
    int id;

    double m,potentialEnergy,kineticEnergy,potentialEnergyBound;

    arma::vec pos,vel;
};
