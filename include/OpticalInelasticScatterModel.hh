#ifndef OpticalInelasticScatterModel_h
#define OpticalInelasticScatterModel_h 1

// Headers
#include "globals.hh"
#include <complex>

using namespace std;

// Class
class OpticalInelasticScatterModel {

public:
  OpticalInelasticScatterModel();
  ~OpticalInelasticScatterModel();

  double computeMFP (double E);
  double sampleELF (double E);
  double calcELFpoint (double E, double omega);
  double calcFirstBorn(double E, double omega);
  double calcSecondBorn(double E, double omega);
  double calcExchange(double E, double omega);
  double calcELF (double k, double omega);
  complex<double> calcDielectricFunction(double k, double omega);

private:
    
};

#endif

