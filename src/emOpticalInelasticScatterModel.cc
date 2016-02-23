#include "emOpticalInelasticScatterModel.hh"
#include "OpticalInelasticScatterModel.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "Randomize.hh"

#include <iostream>
#include <iomanip>

#include <stdlib.h> 
#include <fstream>

#include <cmath>
#include <array>
#include <complex>

#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;
const int nMFP = 50;
const int nELF = 100;
//const int nMFP = 100;
//const int nELF = 200;

G4ThreadLocal gsl_spline *MFPspline= gsl_spline_alloc(gsl_interp_cspline,nMFP);
G4ThreadLocal gsl_interp_accel *MFPaccelerator = gsl_interp_accel_alloc();

const double lowElimit = 35;
const double highElimit = 5000;

G4ThreadLocal double ELFx[nELF+1];
G4ThreadLocal double ELFy[nELF+1];
G4ThreadLocal double ELF[nELF+1][nELF+1];

G4ThreadLocal G4ParticleChangeForLoss* opticalParticleChange = new G4ParticleChangeForLoss;

G4ThreadLocal G4bool isInitialisedOIS = false;

G4NistManager* manOpt = G4NistManager::Instance();

emOpticalInelasticScatterModel::emOpticalInelasticScatterModel(const G4ParticleDefinition*,
                                                     const G4String& nam)
  :G4VEmModel(nam)
{ 
}

emOpticalInelasticScatterModel::~emOpticalInelasticScatterModel()
{
  //delete opticalParticleChange;
}

void emOpticalInelasticScatterModel::Initialise(const G4ParticleDefinition* particle,
                                           const G4DataVector& /*cuts*/)
{
  OpticalInelasticScatterModel* model = new OpticalInelasticScatterModel();

  for (int i = 1 ; i < nELF ; ++i) {
    ELFx[i] = (double)0;
    ELFy[i] = (double)0;
     
    for (int j = 1 ; j < nELF ; ++j) {
      ELF[i][j] = (double)0;
    }
  }
  
  //////////////////////////////////////////////////
  // Build tables 1) Mean Free Path
  //////////////////////////////////////////////////

  G4cout << " *** Building Mean Free Path table *** " << endl;
  
  double Emax = 2100;//1486.71;
  double MFPx[nMFP];
  double MFPy[nMFP];

  MFPx[0] = 0;
  MFPy[0] = 0;
  for (int i = 1 ; i < nMFP ; ++i) {
    MFPx[i] = (i/(double)nMFP)*Emax;
    MFPy[i] = model->computeMFP(MFPx[i]);
    
    //G4cout << "MFP " << MFPx[i] << " " << MFPy[i]*m/nm << G4endl;
  }
  
  //MFPspline = gsl_spline_alloc(gsl_interp_cspline,nMFP);
  gsl_spline_init (MFPspline, MFPx, MFPy, nMFP);
  
  G4cout << " *** Built Mean Free Path table *** " << endl;

  //abort();
  
  //////////////////////////////////////////////////
  // Build tables 2) Energy Loss Function
  //////////////////////////////////////////////////

  G4cout << " *** Building Energy Loss Function table *** " << endl;

  //ELFx[0] = 0;
  //ELFy[0] = 0;
  for (int i = 0 ; i <= nELF ; ++i) {
    double E_loc = (i/(double)nELF)*Emax;
    ELFx[i] = E_loc;
    //G4cout << "ELF," << E_loc;
    for (int j = 1 ; j <= nELF ; ++j) {
      double omega_loc = (j/(double)nELF)*E_loc;
      if (i==(nELF-1)) ELFy[j] = omega_loc;
      
      ELF[i][j]=model->calcELFpoint(E_loc,omega_loc);
      //G4cout << "," << model->calcELFpoint(E_loc,omega_loc);
    }
    //G4cout << G4endl;
  }

  G4cout << " *** Built Energy Loss Function table *** " << endl;
  //abort();

  //////////////////////////////////////////////////
  // Clean up
  //////////////////////////////////////////////////

  delete model;
  
  if (isInitialisedOIS) { return; }
  isInitialisedOIS = true;

  G4cout << "emOpticalInelasticScatterModel::Initialise" << G4endl;
}

G4double emOpticalInelasticScatterModel::
CrossSectionPerVolume(const G4Material* material,
		      const G4ParticleDefinition* particleDefinition,
		      G4double ekin,
		      G4double,
		      G4double)
{
  //G4Material* gold = manOpt->FindOrBuildMaterial("G4_Au");
  G4double particleEnergy = ekin;
  G4double lowerLimit = DBL_MIN;
  
  //if (material != gold) {return lowerLimit;}
  //else if (particleEnergy/eV >= highElimit) {return lowerLimit;}
  if (particleEnergy/eV >= highElimit) {return lowerLimit;}
  else if (particleEnergy/eV <= lowElimit) {return lowerLimit;}
  G4double meanFreePath = gsl_spline_eval(MFPspline,
					  particleEnergy/eV,
					  MFPaccelerator);
  return (1/(meanFreePath*m));
}

void emOpticalInelasticScatterModel::
    SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
			const G4MaterialCutsCouple* couple,
			const G4DynamicParticle* particle, G4double, G4double) {
  //G4cout << "Optical inelastic scattering called" << G4endl;
    
  //G4ThreadLocal G4ParticleChangeForLoss opticalParticleChange;
  //opticalParticleChange = 0;
  opticalParticleChange = GetParticleChangeForLoss();
  
  G4double particleEnergy = particle->GetKineticEnergy()/eV;
  
  //////////////////////////////////////////////////
  // Obtain energy loss 
  //////////////////////////////////////////////////

  double ELFsample[nELF];
  
  // find adjacent x elements
  int index = 0;
  while ((particleEnergy>ELFx[index])&&(index<nELF)) ++index;
  double proportion = ((particleEnergy-ELFx[index-1])
		       /(ELFx[index]-ELFx[index-1]));
  
  // form interpolated ELF
  double ELFmax = 0;
  for (int i = 0 ; i < nELF ; ++i) {
    ELFsample[i] = ELF[index-1][i] + proportion*(ELF[index][i]-ELF[index-1][i]);
    ELFy[i] = ((double)i/nELF)*particleEnergy;
    if (ELFsample[i]>ELFmax) ELFmax = ELFsample[i];
  }

  gsl_interp_accel *ELFaccelerator = gsl_interp_accel_alloc();
  gsl_spline *energyLossFunction = gsl_spline_alloc(gsl_interp_cspline, nELF);

  gsl_spline_init (energyLossFunction, ELFy, ELFsample, nELF);
  
  double dOmega = ELFy[2]-ELFy[1];
  double yValue = 0;
  double omegaOut = 0;
  do {
    double R1 = G4UniformRand();
    omegaOut = R1*particleEnergy;
    yValue = gsl_spline_eval(energyLossFunction, omegaOut, ELFaccelerator);
  } while ((G4UniformRand()*ELFmax)>yValue);  
  double omega = omegaOut;
    
  //////////////////////////////////////////////////
  // Effect paticle changes
  //////////////////////////////////////////////////

  G4double newEnergy = particleEnergy-omega;

  opticalParticleChange->SetProposedKineticEnergy(newEnergy*eV);
  opticalParticleChange->ProposeLocalEnergyDeposit(omega*eV);
    
  /*
  if(!opticalParticleChange->CheckIt(*opticalParticleChange->GetCurrentTrack())){
    G4cout << "oh dear OIS "<< G4endl;
    //abort();
    }*/

  gsl_spline_free(energyLossFunction);
  gsl_interp_accel_free(ELFaccelerator);
}
