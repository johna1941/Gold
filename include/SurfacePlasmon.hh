#ifndef SurfacePlasmon_h
#define SurfacePlasmon_h 1

// Headers
#include "globals.hh"
#include "G4VDiscreteProcess.hh"

using namespace std;

class SurfacePlasmon : public G4VDiscreteProcess {
  public:
  
  SurfacePlasmon(const G4String& processName = "SurfacePlasmon");
  ~SurfacePlasmon();

  G4bool IsApplicable(const G4ParticleDefinition&);
  G4double GetMeanFreePath (const G4Track &aTrack,
  			    G4double previousStepSize,
  			    G4ForceCondition *condition);
  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
  						G4double   previousStepSize,
  						G4ForceCondition* condition);
  G4VParticleChange* PostStepDoIt(const G4Track & aTrack, const G4Step& step);
};

#endif
