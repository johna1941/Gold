#ifndef ElectronCapture_h
#define ElectronCapture_h 1

#include "G4VDiscreteProcess.hh"
#include "globals.hh"
#include "G4ParticleChangeForGamma.hh"

class G4Region;

class ECut : public G4VDiscreteProcess
{
public:

  ECut(const G4String& regName, G4double ekinlimit);

  virtual ~ECut();

  void SetKinEnergyLimit(G4double);

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  virtual G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
							 G4double previousStepSize,
							 G4ForceCondition* condition);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:

  virtual G4double GetMeanFreePath(const G4Track&, G4double,G4ForceCondition*);

private:

  // hide assignment operator as private
  ECut(const ECut&);
  ECut& operator = (const ECut &right);

  G4double kinEnergyThreshold;
  G4String regionName;
  G4Region* region;
};

#endif
