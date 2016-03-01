 #include "PhysicsList.hh"

//#include "G4EmLivermorePhysics.hh"
//#include "G4EmStandardPhysics_option3.hh"
//#include "G4EmStandardPhysics_option4.hh"
//#include "G4EmStandardPhysics_XPS.hh"
//#include "G4EmStandardPhysics_protonDetector.hh"
//#include "G4EmDNAPhysics.hh"

#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicsListHelper.hh"

#include "G4ProcessManager.hh"
#include "G4Decay.hh"

#include "StepMax.hh"
#include "SurfacePlasmon.hh"
#include "ECut.hh"
#include "LifetimeBroadening.hh"

//#include "OpticalInelasticScatter.hh"
//#include "NISTElasticScatter.hh"
#include "emOpticalInelasticScatter.hh"
#include "emNISTElasticScatter.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4VAtomDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4ProductionCutsTable.hh"

#include "G4DummyModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4UAtomicDeexcitation.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4eIonisation.hh"

#include <iostream>
#include <fstream>
using namespace std;

G4double fCutForGamma;
G4double fCutForElectron;
G4double fCutForPositron;

//G4String emName;
//G4VPhysicsConstructor* emPhysicsList;   

//emNISTElasticScatter* NISTElasticScatterPhysics;
//emOpticalInelasticScatter* OIS;
//ECut* energyCut;
//LifetimeBroadening* lifeBroad;
//SurfacePlasmon* surfacePlasmonPhysics;

G4PhotoElectricEffect* pe;
G4VEmModel* theLivermorePEModel;

PhysicsList::PhysicsList() : G4VModularPhysicsList() {
  //G4LossTableManager::Instance();
  
  // EM physics
  //emPhysicsList = new G4EmStandardPhysics_option3();
  //emPhysicsList = new G4EmStandardPhysics_XPS();
  //emPhysicsList = new G4EmStandardPhysics_protonDetector(); 
  //emPhysicsList = new G4EmLivermorePhysics();
  //emPhysicsList = new G4EmDNAPhysics();

  //  surfacePlasmonPhysics = new SurfacePlasmon("surfPlasmon");
}

G4bool isDel = false;

PhysicsList::~PhysicsList()
{
  G4cout << "PhysicsList" << G4endl;

  //  delete OIS;
  //delete NISTElasticScatterPhysics;
  //delete energyCut;
  //delete lifeBroad;
  //  delete surfacePlasmonPhysics;

  delete pe;
  delete theLivermorePEModel;
}

void PhysicsList::ConstructParticle()
{
  //emPhysicsList->ConstructParticle();
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

void PhysicsList::ConstructProcess()
{
  AddTransportation();

  pe = new G4PhotoElectricEffect();
  theLivermorePEModel = new G4LivermorePhotoElectricModel();
  
  setupStandard();
  //setupEm()
  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(50*eV, 400*MeV);
  //SetCutValue(0.1*nm, "gamma");
  //SetCutValue(0.1*nm, "e-");
  //SetCutValue(0.1*nm, "e+");
 
  //////////////////////////////////////////////////
  // Clean up
  //////////////////////////////////////////////////

  //delete OIS;
  
  // Deexcitation
  //G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  //G4LossTableManager::Instance()->SetAtomDeexcitation(de);
  //de->SetFluo(true);
  
  //defaultCutValue = 0.001*nm;
  //SetCutsWithDefault();
}

void PhysicsList::setupStandard()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();
    
    if (particleName == "gamma") {

      G4PhotoElectricEffect* phot = new G4PhotoElectricEffect();
      G4LivermorePhotoElectricModel* photModel
	= new G4LivermorePhotoElectricModel();
      phot->AddEmModel(0, photModel);

      ph->RegisterProcess(phot, particle);
      //ph->RegisterProcess(new G4eIonisation,         particle);

    } else if (particleName == "e-") {

      emOpticalInelasticScatter* OIS;
      OIS = new emOpticalInelasticScatter("inelasScat");
      OIS->SetProcessType(fElectromagnetic);
      OIS->SetProcessSubType(53);
      
      ph->RegisterProcess(OIS, particle);

      emNISTElasticScatter* NISTElasticScatterPhysics;
      NISTElasticScatterPhysics = new emNISTElasticScatter("elasScat");
      NISTElasticScatterPhysics->SetProcessType(fElectromagnetic);
      NISTElasticScatterPhysics->SetProcessSubType(51);
      
      ph->RegisterProcess(NISTElasticScatterPhysics, particle);
      
      LifetimeBroadening* lifeBroad;
      lifeBroad = new LifetimeBroadening("lifeBroad");
      lifeBroad->SetProcessType(fElectromagnetic);
      lifeBroad->SetProcessSubType(6);

      ph->RegisterProcess(lifeBroad, particle);
      
      // // SurfacePlasmon* surfacePlasmonPhysics;
      // // surfacePlasmonPhysics = new SurfacePlasmon("surfPlas");
      // // surfacePlasmonPhysics->SetProcessType(fElectromagnetic);
      // // surfacePlasmonPhysics->SetProcessSubType(13);

      // //ph->RegisterProcess(surfacePlasmonPhysics, particle);
      
      ECut* energyCut;
      energyCut = new ECut("tgt_box",0.0*eV);
      energyCut->SetProcessType(fElectromagnetic);
      energyCut->SetProcessSubType(1);

      ph->RegisterProcess(energyCut, particle);
    }
  }
  
  //delete energyCut;
  //delete lifeBroad;
  //delete surfacePlasmonPhysics;
}
