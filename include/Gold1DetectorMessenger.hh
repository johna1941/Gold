#ifndef Gold1DetectorMessenger_h
#define Gold1DetectorMessenger_h 1

//#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIcommand;
class G4UIcmdWithADouble;

class Gold1DetectorConstruction;

class Gold1DetectorMessenger: public G4UImessenger
{
public:
  Gold1DetectorMessenger(Gold1DetectorConstruction * myDet);
  ~Gold1DetectorMessenger ();
  void SetNewValue(G4UIcommand * command,G4String newValues);
private:
  Gold1DetectorConstruction * myDetector;
  G4UIcommand* fpGold1CommandDirectory;
  G4UIcommand* fpGold1SetDirectory;
//  G4UIcmdWithADouble* fpReflectivityCommand;
};

#endif

