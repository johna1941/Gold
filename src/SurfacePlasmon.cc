#include "SurfacePlasmon.hh"
#include "SurfacePlasmonProbabilityMesh.hh"

//#include <string>
#include <iostream>
#include <fstream>
//#include <cstdio>
#include <math.h>

//#include "NISTElasticScatter.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

#include "G4DataSet.hh"
#include "G4DataVector.hh"

#include "G4CrossSectionDataSet.hh"
#include "G4LinInterpolation.hh"

#include "G4ParallelWorldProcess.hh"
#include "G4TransportationManager.hh"

#include "G4SystemOfUnits.hh"

#include "G4ParticleChangeForLoss.hh"

using namespace std;
  
//////////////////////////////////////////////////
// Constructor / Destructor 
//////////////////////////////////////////////////

SurfacePlasmon::SurfacePlasmon(const G4String& processName)
  : G4VDiscreteProcess(processName) {
}

SurfacePlasmon::~SurfacePlasmon() {
    G4cout << "SurfacePlasmon" << G4endl;
}

//////////////////////////////////////////////////
// General methods
//////////////////////////////////////////////////

G4bool SurfacePlasmon::IsApplicable(const G4ParticleDefinition& particle) {
  return (!strcmp(particle.GetParticleName(),"e-")); 
}

G4double SurfacePlasmon::GetMeanFreePath (const G4Track &aTrack,
					  G4double previousStepSize,
					  G4ForceCondition *condition) {
  //*condition = Forced;
  return DBL_MAX;
}

//////////////////////////////////////////////////
// Interaction length methods
/////////////////////////////////////////////////

G4double SurfacePlasmon::PostStepGetPhysicalInteractionLength (
			         const G4Track &aTrack,
				 G4double previousStepSize,
				 G4ForceCondition *condition) {
  //if (previousStepSize==0.0){
  //  *condition = InActivated;
  //}
  //else {
//    *condition = Forced;
    //}
//  return DBL_MAX;
  /*if (previousStepSize < 0.0) {
    G4cout << "checking1" << G4endl;
    return DBL_MAX;
  } else {
    G4cout << "checking2" << G4endl;
    return 0.0;
    }*/

  *condition = Forced;
  //  if (aTrack.GetStep().GetPostStepPoint()->GetStepStatus()==fGeomBoundary)
  //return DBL_MIN;
  return DBL_MIN;
}

//////////////////////////////////////////////////
// Do it methods
//////////////////////////////////////////////////

G4VParticleChange* SurfacePlasmon::PostStepDoIt(const G4Track & aTrack, const G4Step& step) {
  //  G4cout << aTrack.GetStep()->GetPostStepPoint()->GetStepStatus() << " "
  //	 << fGeomBoundary
  //	 << G4endl;
  G4ParticleChangeForLoss* SPparticleChange = new G4ParticleChangeForLoss;
  SPparticleChange->InitializeForPostStep(aTrack);
  
  if (step.GetPostStepPoint()->GetStepStatus()==fGeomBoundary) { // &&
    //      step.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="tgt_box") {
    //G4cout << "Surface Plasmon called" << G4endl;
    
    //SPparticleChange.ProposeTrueStepLength(0.1*nm);
    // obtain E
    G4double particleEnergy = aTrack.GetKineticEnergy();
    
    // obtain alpha
    G4bool valid;
    CLHEP::Hep3Vector postStepPoint = aTrack.GetStep()
      ->GetPostStepPoint()
      ->GetPosition();
    
    G4int hyperStepNavigatorId = G4ParallelWorldProcess::GetHypNavigatorID();
    
    vector<G4Navigator*>::iterator hyperStepNavigatorIdIterator =
      G4TransportationManager::GetTransportationManager()
      ->GetActiveNavigatorsIterator();
    
    CLHEP::Hep3Vector angleToSurfaceNormal =
      (hyperStepNavigatorIdIterator[hyperStepNavigatorId])
      ->GetGlobalExitNormal(postStepPoint,&valid); // three vector

    if (valid) angleToSurfaceNormal = -angleToSurfaceNormal;
    if (aTrack.GetDynamicParticle()->GetMomentumDirection() *
	angleToSurfaceNormal > 0.0) angleToSurfaceNormal = -angleToSurfaceNormal;
    
    // retreive mesh
    SurfacePlasmonProbabilityMesh* surfacePlasmonProbMesh = new SurfacePlasmonProbabilityMesh(particleEnergy/eV,1,angleToSurfaceNormal.theta());
    //SurfacePlasmonProbabilityMesh* surfacePlasmonProbMesh = new SurfacePlasmonProbabilityMesh(particleEnergy/eV,1,0);
    
    G4double omega = surfacePlasmonProbMesh->getOmega();
    G4double theta = surfacePlasmonProbMesh->getTheta();
    G4double phi = surfacePlasmonProbMesh->getPhi();
    
    CLHEP::Hep3Vector currentMomentumDirection = aTrack.GetDynamicParticle()->GetMomentumDirection();
    CLHEP::Hep3Vector proposedMomentumDirection = currentMomentumDirection;
    proposedMomentumDirection.setTheta(theta);
    proposedMomentumDirection.setPhi(phi);
    
    G4double newEnergy = particleEnergy-50*eV;//omega;
    //G4cout << "omega: " << omega << G4endl;
    if (newEnergy>0&&newEnergy<particleEnergy) {
    //      G4cout << "newEnergy: " << newEnergy << G4endl;
      SPparticleChange->SetProposedKineticEnergy(newEnergy);
    } else return SPparticleChange;
    //cout << "Eloss: " << eLoss/MeV << endl;
    SPparticleChange->SetProposedMomentumDirection(proposedMomentumDirection.unit());

    delete surfacePlasmonProbMesh;
    //G4cout << "Surface Plasmon generated" << G4endl;
    return SPparticleChange;
    } else return SPparticleChange;
}
