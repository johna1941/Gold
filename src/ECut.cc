#include "ECut.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Electron.hh"

#include "G4ParticleChangeForGamma.hh"
#include "G4ParticleChangeForLoss.hh"

G4ThreadLocal G4ParticleChangeForGamma* eCutParticleChange = 0;

ECut::ECut(const G4String& regName, G4double ekinlim)
  : G4VDiscreteProcess("eCapture", fElectromagnetic), kinEnergyThreshold(ekinlim),
    regionName(regName), region(0)
{
  if(regName == "" || regName == "world") { 
    regionName = "DefaultRegionForTheWorld";
  }
  //pParticleChange = &fParticleChange;
  if (!eCutParticleChange) eCutParticleChange = new G4ParticleChangeForGamma;
}

ECut::~ECut() 
{
  delete eCutParticleChange;
}

void ECut::SetKinEnergyLimit(G4double val)
{
  kinEnergyThreshold = val;
  if(verboseLevel > 0) {
    G4cout << "### ECut: Tracking cut E(MeV) = " 
	   << kinEnergyThreshold/MeV << G4endl;
  }
}

void ECut::BuildPhysicsTable(const G4ParticleDefinition&)
{
  region = (G4RegionStore::GetInstance())->GetRegion(regionName);
  if(region && verboseLevel > 0) {
    G4cout << "### ECut: Tracking cut E(MeV) = " 
	   << kinEnergyThreshold/MeV << " is assigned to " << regionName 
	   << G4endl;
  }
}

G4bool ECut::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

G4double 
ECut::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
							G4double, 
							G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  
  G4double limit = DBL_MAX; 
  if (aTrack.GetKineticEnergy() < 50*eV) return 0.0;
  //if (aTrack.GetKineticEnergy() < 1370*eV) return 0.0;
  //else if (aTrack.GetKineticEnergy() > 1440*eV) return 0.0;
  return DBL_MAX;
}

G4VParticleChange* ECut::PostStepDoIt(const G4Track& aTrack, 
						   const G4Step&)
{
  //G4ThreadLocal G4ParticleChangeForGamma eCutParticleChange;//    = new G4ParticleChangeForGamma;
  eCutParticleChange->InitializeForPostStep(aTrack);
  eCutParticleChange->ProposeTrackStatus(fStopAndKill);
  eCutParticleChange->ProposeLocalEnergyDeposit(aTrack.GetKineticEnergy());
  eCutParticleChange->SetProposedKineticEnergy(0.0);

  return eCutParticleChange;
}

G4double ECut::GetMeanFreePath(const G4Track&,G4double,
					    G4ForceCondition*)
{
  return DBL_MAX;
}
