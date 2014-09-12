#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>

#include "System.hpp"

using namespace std;
namespace po = boost::program_options;

/* CONSTANTS */

const double T0 = 0;              // Time to start simulation at

/* METHOD DECLARATIONS */

void rk4Simulation(po::variables_map);
void lfSimulation(po::variables_map);

void progress(int,int,int);

/* MAIN */

int main(int argc, char* argv[]) {
  /*
   * Initialize arg parser and options
   */
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("N", po::value<int>()->default_value(100), "number of particles in cluster")
    ("T", po::value<double>()->default_value(2.0), "how long to run simulation")
    ("dt", po::value<double>()->default_value(0.01), "the time step")
    ("saveEach", po::value<int>()->default_value(1), "only save each timestep")
    ("epsilon", po::value<double>()->default_value(0.15), "value for smoothing function")
    ("noSmooth", po::value<int>()->default_value(0), "use smoothing in forces")
    ("R0", po::value<double>()->default_value(20.0), "radius of the cluster")
    ("meanMass", po::value<double>()->default_value(10.0), "mean mass of objects")
    ("dim", po::value<int>()->default_value(3), "dimension of simulation")
    ("method", po::value<string>()->default_value("rk4"), "scheme to integrate by")
    ("savefile", po::value<string>(), "name of save file")
  ;
  // Create object for storing options and parse command line
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  // Print help message
  if (vm.count("help")) {
    cout << desc << endl;
    return 1;
  }

  if (vm["noSmooth"].as<int>() != 1 && vm["noSmooth"].as<int>() != 0) {
    cout << "noSmooth must be 0 or 1 (True / False)." << endl;
    return 1;
  }

  // Must give savefile
  if (!vm.count("savefile")) {
    cout << "Must give savefile.\n--savefile <filename>" << endl;
    return 1;
  }

  cout << "\nStarting simulation." << endl;

  // Call correct method
  if (strcmp(vm["method"].as<string>().c_str(),"rk4") == 0) {
    rk4Simulation(vm);
  } else if (strcmp(vm["method"].as<string>().c_str(),"lf") == 0) {
    lfSimulation(vm);
  }

  cout << "\nSimulation complete.\n" << endl;

  return 0;
}

/* METHODS */

/*
 * Runs simulation with RK4 scheme
 */
void rk4Simulation(po::variables_map vm) {
  SystemRK4 system = SystemRK4(vm["N"].as<int>(),vm["dim"].as<int>(),vm["R0"].as<double>(),vm["meanMass"].as<double>(),T0);

  if (!vm["noSmooth"].as<int>()) {
    system.setSmoothing(vm["epsilon"].as<double>());
  }
  system.setSaveEach(vm["saveEach"].as<int>());

  system.open(vm["savefile"].as<string>());

  int n = vm["T"].as<double>() / vm["dt"].as<double>();

  system.saveState();

  for (int i = 0; i*vm["dt"].as<double>() <= vm["T"].as<double>(); i++) {
    system.advance(vm["dt"].as<double>());
    progress(i, n, 40);
  }

  system.close();
}

/*
 * Runs simulation with leap frog scheme
 */
void lfSimulation(po::variables_map vm) {
  SystemLeapFrog system = SystemLeapFrog(vm["N"].as<int>(),vm["dim"].as<int>(),vm["R0"].as<double>(),vm["meanMass"].as<double>(),T0);

  if (!vm["noSmooth"].as<int>()) {
    system.setSmoothing(vm["epsilon"].as<double>());
  }
  system.setSaveEach(vm["saveEach"].as<int>());

  system.open(vm["savefile"].as<string>());

  int n = vm["T"].as<double>() / vm["dt"].as<double>();

  system.saveState();

  for (int i = 0; i*vm["dt"].as<double>() <= vm["T"].as<double>(); i++) {
    system.advance(vm["dt"].as<double>());
    progress(i, n, 40);
  }

  system.close();
}

/*
 * Function prints progress as nice bar, also skips to beginning of line so
 * function can be used to print again showing a dynamically updating progress
 * bar.
 *
 * @param currentStep The current step number.
 * @param totalSteps The total number of step sizes.
 * @param length Total number of # in progress bar.
 */
void progress(int currentStep, int totalSteps, int length) {
  double fraction = currentStep / ( (double) totalSteps );
  double percent = fraction*100;
  double numberOfSquares = fraction*length;
  double counter = 0;

  cout << "\r[";
  while (counter <= numberOfSquares) {
    cout << "#";
    counter++;
  }
  while (counter <= length) {
    cout << " ";
    counter++;
  }
  cout << "] " << setprecision(1) << fixed << percent << " \% completed ";
  cout.flush();
}
