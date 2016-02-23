#ifndef emNISTElasticScatterModel_h
#define emNISTElasticScatterModel_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4LogLogInterpolation.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NistManager.hh"
#include "G4VAtomDeexcitation.hh"

#include "G4ParticleChangeForGamma.hh"
#include "G4ParticleChangeForLoss.hh"

class emNISTElasticScatterModel : public G4VEmModel
{

public:

  emNISTElasticScatterModel(const G4ParticleDefinition* p = 0, 
		          const G4String& nam = "emNISTElasticScatterModel");

  virtual ~emNISTElasticScatterModel();

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

  G4double RandomTheta(G4double E);

protected:

  //G4ThreadLocal G4ParticleChangeForGamma NISTparticleChange;
  //G4ParticleChangeForLoss* NISTparticleChange;
  
private:
  
  //G4VAtomDeexcitation*      fAtomDeexcitation;
  //G4Material* nistSi;
  //G4bool isInitialised;
  
  //emNISTElasticScatterModel & operator=(const emNISTElasticScatterModel &right);
  //emNISTElasticScatterModel(const emNISTElasticScatterModel&);

};

#endif
