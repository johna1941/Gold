#ifndef LifetimeBroadening_h
#define LifetimeBroadening_h 1

#include "G4VDiscreteProcess.hh"
#include "globals.hh"

class G4Region;

class LifetimeBroadening : public G4VDiscreteProcess
{
public:

  LifetimeBroadening(const G4String& processName = "lifeBroad");
~LifetimeBroadening();

  void SetKinEnergyLimit(G4double);

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  virtual G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
							 G4double previousStepSize,
							 G4ForceCondition* condition);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:

  virtual G4double GetMeanFreePath(const G4Track&, G4double,G4ForceCondition*);

private:

  // hide assignment operator as private
  LifetimeBroadening(const LifetimeBroadening&);
  LifetimeBroadening& operator = (const LifetimeBroadening &right);

};

#endif
