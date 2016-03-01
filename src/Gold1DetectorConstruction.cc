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
// $Id: Gold1DetectorConstruction.cc 90623 2015-06-05 09:24:30Z gcosmo $
//
/// \file Gold1DetectorConstruction.cc
/// \brief Implementation of the Gold1DetectorConstruction class

#include "Gold1DetectorConstruction.hh"
#include "Gold1DetectorMessenger.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4TriangularFacet.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "Gold1SensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
// The above three are for the creation of ABS

#include <fstream>

Gold1DetectorConstruction::Gold1DetectorConstruction()
: fpDetectorMessenger(new Gold1DetectorMessenger(this))
{ }

Gold1DetectorConstruction::~Gold1DetectorConstruction()
{
 delete fpDetectorMessenger; 
}

G4VPhysicalVolume* Gold1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* gold = nist->FindOrBuildMaterial("G4_Au");

  G4String name;

  // World
  name = "World";
  G4Box* solidWorld =
    new G4Box(name,        //its name
       1.*m, 1.*m, 1.*m);     //its size
  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        name);               //its name
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      name,                  //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  // Gold block
  name = "Gold_block";
  G4VSolid* gold_block = new G4Box(name,5.*cm,5.*cm,5.*cm);
  G4LogicalVolume* gold_block_lv = new G4LogicalVolume(gold_block,gold,name);
  new G4PVPlacement(G4Translate3D(),gold_block_lv,name,logicWorld,0,false,checkOverlaps);

  return physWorld;
}

void Gold1DetectorConstruction::ConstructSDandField()
{
//  G4SDManager* pSDman = G4SDManager::GetSDMpointer();
//  G4VSensitiveDetector* fibreSD = new Gold1SensitiveDetector("Fibre");
//  pSDman->AddNewDetector(fibreSD);
//  fFibreLV->SetSensitiveDetector(fibreSD);
}

