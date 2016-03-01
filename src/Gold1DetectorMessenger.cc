#include "Gold1DetectorMessenger.hh"

#include "Gold1DetectorConstruction.hh"
#include "G4UIparameter.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADouble.hh"
#include <sstream>

Gold1DetectorMessenger::Gold1DetectorMessenger(Gold1DetectorConstruction * myDet)
:myDetector(myDet)
{
  fpGold1CommandDirectory = new G4UIcommand("/Gold/",this);
  fpGold1CommandDirectory->SetGuidance("Gold1 detector control.");

  fpGold1SetDirectory = new G4UIcommand("/Gold/set/",this);
  fpGold1SetDirectory->SetGuidance("Set commands.");

//  fpReflectivityCommand = new G4UIcmdWithADouble("/Gold/set/reflectivity",this);
//  fpReflectivityCommand->SetGuidance("Define reflectivity of chamber walls.");
//  fpReflectivityCommand->SetToBeBroadcasted(false);
//  fpReflectivityCommand->AvailableForStates(G4State_PreInit);
}

Gold1DetectorMessenger::~Gold1DetectorMessenger () {
//  delete fpReflectivityCommand;
  delete fpGold1SetDirectory;
  delete fpGold1CommandDirectory;
}

void Gold1DetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
//  if (command == fpReflectivityCommand)
//  {
//    myDetector->fReflectivity = fpReflectivityCommand->GetNewDoubleValue(newValues);
//  }
}
