#ifndef emOpticalInelasticScatterModel_h
#define emOpticalInelasticScatterModel_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4LogLogInterpolation.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NistManager.hh"
#include "G4VAtomDeexcitation.hh"

#include "G4ParticleChangeForGamma.hh"
#include "G4ParticleChangeForLoss.hh"

class emOpticalInelasticScatterModel : public G4VEmModel
{

public:

  emOpticalInelasticScatterModel(const G4ParticleDefinition* p = 0, 
		          const G4String& nam = "emOpticalInelasticScatterModel");

  virtual ~emOpticalInelasticScatterModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double CrossSectionPerVolume(const G4Material* material,
					   const G4ParticleDefinition* p,
					   G4double ekin,
					   G4double emin,
					   G4double emax);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);


protected:

  //G4ParticleChangeForGamma* opticalParticleChange;
  //G4ThreadLocal G4ParticleChangeForLoss opticalParticleChange;
  
private:
  
  //G4VAtomDeexcitation*      fAtomDeexcitation;
  //G4Material* nistSi;
  //G4bool isInitialised;
  
  // emOpticalInelasticScatterModel & operator=(const emOpticalInelasticScatterModel &right);
  //emOpticalInelasticScatterModel(const emOpticalInelasticScatterModel&);

};

#endif
