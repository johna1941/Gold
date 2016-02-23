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
// $Id: Gold1ActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file Gold1ActionInitialization.cc
/// \brief Implementation of the Gold1ActionInitialization class

#include "Gold1ActionInitialization.hh"
#include "Gold1PrimaryGeneratorAction.hh"
#include "Gold1RunAction.hh"
#include "Gold1EventAction.hh"
#include "Gold1SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Gold1ActionInitialization::Gold1ActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Gold1ActionInitialization::~Gold1ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gold1ActionInitialization::BuildForMaster() const
{
  SetUserAction(new Gold1RunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gold1ActionInitialization::Build() const
{
  SetUserAction(new Gold1PrimaryGeneratorAction);
  SetUserAction(new Gold1RunAction);

  Gold1EventAction* p1EventAction = new Gold1EventAction;
  SetUserAction(p1EventAction);
  
  SetUserAction(new Gold1SteppingAction(p1EventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
