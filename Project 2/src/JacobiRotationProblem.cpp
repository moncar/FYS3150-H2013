#include <cmath>
#include <fstream>

#include "JacobiRotationProblem.hpp"

/*
 * Constructor
 *
 * Saves matrix in the two pointers.
 *
 * @param newProblem The matrix to find eigenvalues for.
 */
JacobiRotationProblem :: JacobiRotationProblem(mat newProblem) {
  matrix = newProblem;
  originalMatrix = newProblem;

  // Set special fields
  finished = false;
  rotations = 0;
}

/*
 * Uses the algorithm to rearange the matrix into a diagonal matrix which will
 * contain the eigenvalues.
 *
 * @param errorTolerance The error tolerance to make all non-diagonal elements
 * smaller than.
 */
void JacobiRotationProblem :: solve(double errorTolerance) {
  double a_ik,a_il,a_ll,a_kk,a_kl;

  // Decide interval of progress update from matrix size
  int progUpdate = (int) ( matrix.n_rows / 2.0 ) + 1;
  int counter = 0;  // For counting progress updates

  do {
    // Update max non diag elm and coors
    updateMaxNonDiagElm();

    // Fetch certain elements
    a_ll = matrix(l,l);
    a_kk = matrix(k,k);
    a_kl = matrix(k,l);

    // Find new rotation values
    updateRotationValues();

    matrix(k,l) = 0.0;
    matrix(l,k) = 0.0;

    // Perform rotation
    matrix(k,k) = a_kk*c*c - 2*a_kl*c*s + a_ll*s*s;
    matrix(l,l) = a_ll*c*c + 2*a_kl*c*s + a_kk*s*s;

    for (int i = 0; i < matrix.n_rows; i++) {
      if (i == k || i == l) {
        // Skip diagonal elements
        continue;
      }
      a_ik = matrix(i,k);
      a_il = matrix(i,l);

      matrix(i,k) = a_ik*c - a_il*s;
      matrix(k,i) = matrix(i,k);

      matrix(i,l) = a_il*c + a_ik*s;
      matrix(l,i) = matrix(i,l);
    }

    // Upp the number of rotations
    rotations++;

    /*
     * Display progress. Progess is shown in current max non diag element and
     * difference to the error tolerance.
     */
    if (counter == progUpdate) {
      counter = 0;
      cout << "\rCur: " << abs(maxNonDiagElm) << " Diff: " <<
        (abs(maxNonDiagElm) - errorTolerance) << "             ";
      cout.flush();
    }
    counter++;

  } while (abs(maxNonDiagElm) > errorTolerance);

  // Progress must end the line
  cout << endl;

  finished = true;
}

/*
 * Find the rotation matrix from the current state of the matrix.
 */
void JacobiRotationProblem :: updateRotationValues() {
  // Values needed for rotation
  double tau = ( matrix(l,l) - matrix(k,k) ) / ( 2.0 * matrix(k,l) );

  if (tau > 0) {
    t = 1.0 / (tau + sqrt(1 + tau*tau));
  } else {
    t = 1.0 / (tau - sqrt(1 + tau*tau));
  }

  c = 1.0 / sqrt(1 + t*t);
  s = c * t;
}

/*
 * Updates the current max non diag element and also gives coordinates.
 */
void JacobiRotationProblem :: updateMaxNonDiagElm() {
  /*
   * This uses the assumption that given matrix at least has to rows. This is a
   * reasonable assumption.
   */
  maxNonDiagElm = matrix(1,0);
  k = 1; l = 0;

  for (int j = 0; j < matrix.n_rows; j++) {
    for (int i = j+1; i < matrix.n_cols; i++) {
      if (j == i) {
        // Skip elements along diagonal
        continue;
      }
      if (abs(matrix(j,i)) > abs(maxNonDiagElm)) {
        // Store coordinates
        k = j;
        l = i;

        // Store max element
        maxNonDiagElm = matrix(j,i);
      }
    }
  }
}

/*
 * Prints the complete matrix.
 */
void JacobiRotationProblem :: printResultMatrix() {
  if (!isFinished()) {
    return;
  }

  cout << "Printing resulting matrix." << endl << "Found after " << rotations
    << " rotations." << endl;
  cout << matrix << endl;
}

/*
 * Prints the eigenvalues to a file named by given argument.
 *
 * @param filename The name of the outputfile, full name.
 */
void JacobiRotationProblem :: saveResult(string filename) {
  if (!isFinished()) {
    cout << "Tried to print result when solver hasn't run." << endl;
    return;
  }
  
  ofstream outfile;
  outfile.open(filename.c_str());
  vec eigenvalues = matrix.diag();
  eigenvalues = sort(eigenvalues);

  // Write header with metadata
  if (hasParams) {
      outfile << "******[META]******" << endl;
      outfile << "#e: " << numElectrons << endl;
      outfile << "rhoMax: " << rhoMax << endl;
      outfile << "n: " << matrix.n_rows << endl;
      if (numElectrons > 1) {
        outfile << "omegaR: " << omegaR << endl;
      }
      outfile << "***[SORTED EIGENVALUES]***" << endl;
  }

  for (int i = 0; i < matrix.n_rows; i++) {
    outfile << eigenvalues(i) << endl;
  }
  outfile.close();
  cout << "Saved eigenvalues to file: " << filename << endl;
}

/*
 * Gives physical parameters for the problem to the solver. Used for printing
 * parameters in the datafile.
 *
 * NOTE: Bad programming convention. Only relevant for Schrodinger equation
 * problem.
 *
 * @param rMax rhoMax, max dimensionless radius.
 * @param numElec The number of electrons.
 */
void JacobiRotationProblem :: giveParameters(double rMax, int numElec) {
  rhoMax = rMax;
  numElectrons = numElec;
  hasParams = true;
}

/*
 * Overloading of above function for two electrons.
 *
 * @param wr OmegaR the strength of the harmonic oscillator potential.
 */
void JacobiRotationProblem :: giveParameters(double rhoMax, int numElec, double wr) {
  giveParameters(rhoMax,numElec);
  omegaR = wr;
}

/** Get methods **/

/*
 * Returns the number of rotations needed to solve the problem. If solution is
 * not found yet it returns 0.
 */
int JacobiRotationProblem :: getNumRotations() {
  if (!isFinished()) {
    return 0;
  }
  return rotations;
}

/*
 * Returns status of problem. True if solution is found and false if not.
 */
bool JacobiRotationProblem :: isFinished() {
  return finished;
}
