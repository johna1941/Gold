//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: Gold1SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file Gold1SteppingAction.cc
/// \brief Implementation of the Gold1SteppingAction class

#include "Gold1SteppingAction.hh"
#include "Gold1EventAction.hh"
#include "Gold1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4Electron.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Gold1SteppingAction::Gold1SteppingAction(Gold1EventAction* /*eventAction*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Gold1SteppingAction::~Gold1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::ofstream peOutFile("photo_electron_kinetic_energies.dat");
std::ofstream escapeOutFile("escape_electron_kinetic_energies.dat");
G4Mutex outFileMutex = G4MUTEX_INITIALIZER;

void Gold1SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4StepPoint* postStepPoint = step->GetPostStepPoint();

  G4String processName = postStepPoint->GetProcessDefinedStep()->GetProcessName();
  if (processName == "phot") {
    const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();
    for (const auto& track: *secondaries) {
      const G4ParticleDefinition* secondaryParticleDefinition = track->GetParticleDefinition();
      if (secondaryParticleDefinition == G4Electron::Electron()) {
        G4double kineticEnergy = track->GetKineticEnergy();
        G4MUTEXLOCK(&outFileMutex);
        peOutFile << kineticEnergy/eV << std::endl;
        G4MUTEXUNLOCK(&outFileMutex);
      }
    }
    return;
  }

  G4Track* track = step->GetTrack();
  if (track->GetDefinition() == G4Electron::Electron()) {
    G4String postStepPVName = postStepPoint->GetTouchableHandle()->GetVolume()->GetName();
    if (postStepPVName == "World") {
      G4double kineticEnergy = track->GetKineticEnergy();
      G4MUTEXLOCK(&outFileMutex);
      escapeOutFile << kineticEnergy/eV << std::endl;
      G4MUTEXUNLOCK(&outFileMutex);
    }
    return;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

