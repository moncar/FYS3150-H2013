#include <string>
#include <armadillo>
#include <fstream>

class CelestialObject {
  public:
    CelestialObject(std::string,arma::vec,arma::vec,double,std::string);

    arma::vec getDistTo(CelestialObject);
    arma::vec getPos(),getV(),getForce();

    std::string getId();

    double getM();

    void setF(arma::vec);
    void setPos(arma::vec);
    void setV(arma::vec);

    void saveCurrentPos();
    void close();

  private:
    double m;

    arma::vec pos,v,F;

    std::string id;
    std::ofstream* datafile;
};
