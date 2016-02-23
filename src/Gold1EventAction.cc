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
// $Id: Gold1EventAction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file Gold1EventAction.cc
/// \brief Implementation of the Gold1EventAction class

#include "Gold1EventAction.hh"
#include "Gold1Run.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

Gold1EventAction::Gold1EventAction()
: fNumberOfPhotons(0)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Gold1EventAction::~Gold1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gold1EventAction::BeginOfEventAction(const G4Event*)
{    
  fNumberOfPhotons = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gold1EventAction::EndOfEventAction(const G4Event*)
{   
  // accumulate statistics in Gold1Run
  Gold1Run* run =
  static_cast<Gold1Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->AddPhotons(fNumberOfPhotons);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
