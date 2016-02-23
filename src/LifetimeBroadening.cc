#include "LifetimeBroadening.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Electron.hh"

#include "G4ParticleChangeForGamma.hh"
#include "G4ParticleChangeForLoss.hh"

#include <vector>
#include <array>
#include <random>
#include <cmath>

using namespace std;

/*const array <double,22> livermoreBE = {0.08096, // K
				       0.014336, // L1
				       0.013776, // L2
				       0.011922, // L3
				       0.0034014, // M1
				       0.0031443, // M2
				       0.002734, // M3
				       0.0023007, // M4
				       0.0022122, // M5
				       0.0007474, // N1
				       0.00063628, // N2
				       0.00053739, // N3
				       0.00035348, // N4
				       0.00033505, // N5
				       9.737e-05, // O1
				       9.339e-05, // N6
				       0.00011523, // N7
				       7.855e-05, // O2
				       6.137e-05, // O3
				       1.216e-05, // O4
				       1.046e-05, // O5
				       8.3e-06 // P1
};

const array <double,22> lineWidths = {52.1, // K
				      9.8, // L1
				      5.53, // L2
				      5.53, // L3
				      15, // M1
				      9.5, // M2
				      8.5, // M3
				      2.18, // M4
				      2.18, // M5
				      8.5, // N1
				      6.4, // N2
				      5.05, // N3
				      4.1, // N4
				      3.9, // N5
				      0, // O1
				      0.37, // N6
				      0.33, // N7
				      0, // O2
				      0, // O3
				      0, // O4
				      0, // O5
				      0 // P1
				      };*/

const array <double,19> livermoreBE = {0.080725, // K
				       0.014353, // L1
				       0.013734, // L2
				       0.011919, // L3
				       0.003425, // M1
				       0.003148, // M2
				       0.002743, // M3
				       0.002291, // M4
				       0.002206, // M5
				       0.007621, // N1
				       0.0006427, // N2
				       0.0005463, // N3
				       0.0003532, // N4
				       0.0003351, // N5
				       0.0001072, // O1
				       0.0000876, // N6
				       0.000084, // N7
				       0.0000742, // O2
				       0.0000572 // O3
};

// line widths taken from campbell and papp 2001 (0s where no data)
const array <double,19> lineWidths = {52.1, // K
				      9.8, // L1
				      5.53, // L2
				      5.53, // L3
				      15, // M1
				      9.5, // M2
				      8.5, // M3
				      2.18, // M4
				      2.18, // M5
				      8.5, // N1
				      6.4, // N2
				      5.05, // N3
				      4.1, // N4
				      3.9, // N5
				      0, // O1
				      0.37, // N6
				      0.33, // N7
				      0, // O2
				      0, // O3
};

//G4ThreadLocal G4ParticleChangeForLoss* lifetimeParticleChange =
//  new G4ParticleChangeForLoss;
G4ThreadLocal G4ParticleChangeForGamma* lifetimeParticleChange =
  new G4ParticleChangeForGamma;
//G4ThreadLocal G4ParticleChange* lifetimeParticleChange = new G4ParticleChange;

LifetimeBroadening::LifetimeBroadening(const G4String& processName)
  : G4VDiscreteProcess(processName)
{
  //pParticleChange = &fParticleChange;
}

LifetimeBroadening::~LifetimeBroadening() 
{
  //delete lifetimeParticleChange;
}


G4bool LifetimeBroadening::IsApplicable(const G4ParticleDefinition& particle)
{
  return (!strcmp(particle.GetParticleName(),"e-")); 
}

G4double 
LifetimeBroadening::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
							G4double, 
							G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;

  G4double particleEnergy = aTrack.GetKineticEnergy();
  G4String materialName = aTrack.GetMaterial()->GetName();

  //if (materialName != "G4_Au") return DBL_MAX;
  //else

  //cout << aTrack.GetCurrentStepNumber() << endl;
  //const G4String creatorProcessName
  //  = *aTrack.GetCreatorProcess()->GetProcessName();

  if(aTrack.GetCreatorProcess()!=0) {
    if (aTrack.GetCreatorProcess()->GetProcessName()=="phot" &&
	aTrack.GetCurrentStepNumber()==1) return 0.0;
    else return DBL_MAX;
  } else return DBL_MAX;
  //rack.GetCreatorProcess()->GetProcessName().compareTo("phot");
  
  /*  for (int i = 0 ; i < livermoreBE.size() ; ++i) {
    double error = (livermoreBE[i]-((1487.6*eV)-particleEnergy));
    if (abs(error)<(1*eV)) return 0.0;
    //cout << livermoreBE[i] << " " << ((1487.6*eV)-particleEnergy) << " " << error << endl;
    }*/
  //cout << "oli" << endl;
  return DBL_MAX;
}

G4VParticleChange* LifetimeBroadening::PostStepDoIt(const G4Track& aTrack, 
						   const G4Step&)
{
  //G4cout << "particle cut!" << G4endl;
  //G4ParticleChange* lifetimeParticleChange = new G4ParticleChange;
  lifetimeParticleChange->InitializeForPostStep(aTrack);
  //lifetimeParticleChange->Initialize(aTrack);
  
  double particleEnergy = aTrack.GetKineticEnergy();

  double lineWidth = 0.0;
  double bestError = DBL_MAX;
  bool found = false;
  int index = 0;
  for (int i = 0 ; i < livermoreBE.size() ; ++i) {
    //double error = (livermoreBE[i]-((1487.6*eV)-particleEnergy));
    //if (abs(error)<(1*eV)) lineWidth = lineWidths[i];

    double currentError = (livermoreBE[i]-((1487.6*eV)-particleEnergy));
    //cout << livermoreBE[i] << " " << particleEnergy << " " << (livermoreBE[i]-((1487.6*eV)-particleEnergy)) << endl;
    if (abs(currentError)<bestError && abs(currentError<10*eV)) {
      index = i;
      bestError = abs(currentError);
      found = true;
    }
  }
  if (found) lineWidth = lineWidths[index];

  //default_random_engine generator;
  //normal_distribution<double> distribution(5.0,2.0);
  //double newEnergy = distribution(generator);

  std::random_device rd;
  std::mt19937 e2(rd());
  std::cauchy_distribution<> dist(particleEnergy/eV, lineWidth/2); // 2*gamma=FWHM
  double newEnergy = dist(e2);

  //////////////////////////////////////////////////
  // shouldn't really need this...
  //////////////////////////////////////////////////
  while(newEnergy<0.0) newEnergy = dist(e2);
  //G4cout << "rand " << newEnergy << " " << lineWidth << G4endl;


  lifetimeParticleChange->SetProposedKineticEnergy(newEnergy*eV);
  //lifetimeParticleChange->ProposeEnergy(newEnergy*eV);
  //cout << "lifetime woop " << particleEnergy/eV << " " << newEnergy << " " << lineWidth<< endl;

  CLHEP::Hep3Vector currentMomentumDirection = aTrack.GetDynamicParticle()->GetMomentumDirection();
  CLHEP::Hep3Vector proposedMomentumDirection = currentMomentumDirection;
  //lifetimeParticleChange->
  //SetProposedMomentumDirection(proposedMomentumDirection.unit());

  /*if (!lifetimeParticleChange->CheckIt(*lifetimeParticleChange->GetCurrentTrack())) {
    G4cout << "oh dear lifetime " << particleEnergy/eV
	   << " " << lineWidth/2
	   << " " << newEnergy << G4endl;
	   };*/

  if(newEnergy<0.0) G4cout << "urgh lifetime "
			   << particleEnergy/eV << " "
			   << newEnergy << " "
			   << lineWidth/2 << " "
			   << dist(e2) 
			   << G4endl;
  
  return lifetimeParticleChange;
}

G4double LifetimeBroadening::GetMeanFreePath(const G4Track&,G4double,
					    G4ForceCondition*)
{
  return DBL_MAX;
}
