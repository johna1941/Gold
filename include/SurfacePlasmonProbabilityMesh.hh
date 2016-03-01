#ifndef SurfacePlasmonProbabilityMesh_h
#define SurfacePlasmonProbabilityMesh_h 1

// Headers
#include <cmath>
#include <array>
#include <complex>

#include "globals.hh"
#include "G4VDiscreteProcess.hh"

using namespace std;

class SurfacePlasmonProbabilityMesh {
public:
  struct outputs {
    G4double omegaOut;
    G4double thetaOut;
    G4double phiOut;
  };
  
public:
  SurfacePlasmonProbabilityMesh(double E_eV, double s_m, double alpha);
  ~SurfacePlasmonProbabilityMesh();

  complex<double> calcRc (double omega, double v, double s,
			complex<double> epsilon, double k_s);
  complex<double> calcDielectricFunction(double k, double omega);
  double calcPs (double k_s, double v, double alpha, double k,
		 complex<double> epsilon, complex<double> Rc,
		 array<double,3> v_vec, array<double,3> k_vec,
		 double omega);

  G4double getOmega ();
  G4double getTheta ();
  G4double getPhi ();
};

#endif
