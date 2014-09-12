// Are these two lines needed in this file?
#include <armadillo>
#include <stdlib.h>
using namespace arma;
using namespace std;

/*
 * Represents a problem that uses the Jacobi rotation algorithm to find the
 * eigenvalues of a matrix.
 *
 * Functionality should include just giving the matrix and letting it solve for
 * itself without interference. Should be able to run solver several times for
 * different error tolerances.
 */
class JacobiRotationProblem {
  public:
    /*
     * Functions
     */
    JacobiRotationProblem(mat);         // Constructor only takes matrix to solve
    int getNumRotations();              // Number of rotations used to solve
    vec* getEigenvectors();             // Returns eigenvectors
    double* getEigenvalues();           // Returns eigenvalues
    void solve(double);                 // Start algorithm, arg is error tolerance
    void printResultMatrix();           // If finished, prints complete matrix
    void saveResult(string);            // If finished, saves eigenvalues to file

    /*
     * For giving physical parameters
     */
    void giveParameters(double,int);
    void giveParameters(double,int,double);

  private:
    /*
     * Fields
     */
    bool finished,hasParams;            // Switch telling if solver has run
    int rotations,k,l;                  // Int fields for coors and rotations
    mat matrix,originalMatrix;          // Matrix fields
    double t,c,s,maxNonDiagElm;         // For storing rotation values and maxelm

    /*
     * Physical parameters
     */
    double rhoMax,omegaR;
    int numElectrons;

    /*
     * Functions
     */
    void updateRotationValues();        // Finds rotation values
    void updateMaxNonDiagElm();         // Updates max non diag elm (and coors)
    bool isFinished();                  // If rotations is run
};
