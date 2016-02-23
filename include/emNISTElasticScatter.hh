#ifndef emNISTElasticScatter_h
#define emNISTElasticScatter_h 1

#include "G4VEmProcess.hh"
#include "G4Electron.hh"

// Available models
#include "emNISTElasticScatterModel.hh"

class emNISTElasticScatter : public G4VEmProcess

{
public: 

  emNISTElasticScatter(const G4String& processName
			        = "emNISTElasticScatter",
			      G4ProcessType type
			        = fElectromagnetic);

  virtual ~emNISTElasticScatter();

  virtual G4bool IsApplicable(const G4ParticleDefinition&);
  
  virtual void PrintInfo();

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*);

private:

};

#endif
