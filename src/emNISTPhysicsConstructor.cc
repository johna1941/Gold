// JA 26-Feb-22016

#include "emNISTPhysicsConstructor.hh"

#include "emNISTElasticScatter.hh"
#include "emOpticalInelasticScatter.hh"

void emNISTPhysicsConstructor::ConstructProcess()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  aParticleIterator->reset();
  while( (*aParticleIterator)() ){
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4String particleName = particle->GetParticleName();
    if (particleName == "e-") {
      ph->RegisterProcess(new emNISTElasticScatter, particle);
      ph->RegisterProcess(new emOpticalInelasticScatter, particle);
    }
  }
}
