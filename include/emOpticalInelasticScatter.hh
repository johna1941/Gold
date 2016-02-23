#ifndef emOpticalInelasticScatter_h
#define emOpticalInelasticScatter_h 1

#include "G4VEmProcess.hh"
#include "G4Electron.hh"

// Available models
#include "emOpticalInelasticScatterModel.hh"

class emOpticalInelasticScatter : public G4VEmProcess

{
public: 

  emOpticalInelasticScatter(const G4String& processName
			        = "emOpticalInelasticScatter",
			      G4ProcessType type
			        = fElectromagnetic);

  virtual ~emOpticalInelasticScatter();

  virtual G4bool IsApplicable(const G4ParticleDefinition&);
  
  virtual void PrintInfo();

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*);

private:
};

#endif
