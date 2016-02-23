#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class StepMax;
class PhysicsListMessenger;
class NISTElasticScatter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
public:

  PhysicsList();
 ~PhysicsList();

  virtual void ConstructParticle();

  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);
        
  virtual void ConstructProcess();
    
  void AddSurfacePlasmons();
  void setupStandard();
  void setupEm();  

private:
};

#endif
