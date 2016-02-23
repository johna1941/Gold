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
// $Id: Gold1RunAction.cc 89630 2015-04-23 12:11:28Z gcosmo $
//
/// \file Gold1RunAction.cc
/// \brief Implementation of the Gold1RunAction class

#include "Gold1RunAction.hh"
#include "Gold1PrimaryGeneratorAction.hh"
#include "Gold1DetectorConstruction.hh"
#include "Gold1Run.hh"
#include "G4GeneralParticleSource.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

#include "Gold1Analysis.hh"


// I believe analysis manager is now redundent since I'm using <fstream>

Gold1RunAction::Gold1RunAction():G4UserRunAction()
{
    // Create analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetVerboseLevel(0);
    // Create ntuple
  analysisManager->CreateNtuple("Gold1", "NumberOfPhotons");
  analysisManager->CreateNtupleDColumn("NPhot");
  analysisManager->FinishNtuple();
  numberOfPhotons = 0;
}

Gold1RunAction::~Gold1RunAction()
{
    delete G4AnalysisManager::Instance();
}

G4Run* Gold1RunAction::GenerateRun()
{
  return new Gold1Run; 
}

void Gold1RunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    
  // Get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Open an output file
    analysisManager->OpenFile("Gold1");
}

void Gold1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  const Gold1Run* p1Run = static_cast<const Gold1Run*>(run);
  numberOfPhotons = p1Run->GetPhotons();

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const Gold1PrimaryGeneratorAction* generatorAction
  = static_cast<const Gold1PrimaryGeneratorAction*>
  (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4GeneralParticleSource* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }

  // Print
  //
  if (IsMaster()) {
    G4cout
    << "\n--------------------End of Global Run-----------------------";
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
     
      std::ofstream output("test.txt", std::ios::app); // std::ios::app = append mode: will add data onto the end of file
      output << numberOfPhotons << std::endl;
	output.close();
	      
  }
  else {
    G4cout
    << "\n--------------------End of Local Run------------------------";
  }

  G4cout
  << "\n The run consists of " << nofEvents << " " << runCondition
  << "\n Number of photons reaching sensitive detector: "
  << numberOfPhotons
  << "\n------------------------------------------------------------"
  << G4endl;
}
