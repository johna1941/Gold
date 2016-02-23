#include "emNISTElasticScatter.hh"
#include "emNISTElasticScatterModel.hh"
#include "G4DummyModel.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4TrackStatus.hh"

#include "G4VEmModel.hh"
#include "G4ParticleChangeForLoss.hh"

#include "G4DataSet.hh"
#include "G4DataVector.hh"

#include "G4CrossSectionDataSet.hh"
#include "G4LinInterpolation.hh"

#include "G4ProcessType.hh"

#include "G4Threading.hh"
#include "G4ThreadLocalSingleton.hh"

#include <iostream>
#include <fstream>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;

G4ThreadLocal G4bool NISTisInitialised = false;

emNISTElasticScatter::emNISTElasticScatter(const G4String& processName,
                                           G4ProcessType type)
  :G4VEmProcess (processName, type)
{
  SetProcessSubType(51);
}
 
emNISTElasticScatter::~emNISTElasticScatter()
{}

G4bool emNISTElasticScatter::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Electron::Electron());
}

void emNISTElasticScatter::InitialiseProcess(const G4ParticleDefinition* p)
{
  if(!NISTisInitialised) {
    NISTisInitialised = true;
    SetBuildTableFlag(false);
    G4String name = p->GetParticleName();

    if(!EmModel(1)) SetEmModel(new emNISTElasticScatterModel, 1);
    AddEmModel(2, EmModel(1));   
  } 
}

void emNISTElasticScatter::PrintInfo()
{
}         
