// JA 26-Feb-22016

#ifndef EMNISTPHYSICSCONSTRUCTOR_HH
#define EMNISTPHYSICSCONSTRUCTOR_HH

#include "G4VPhysicsConstructor.hh"

class emNISTPhysicsConstructor: public G4VPhysicsConstructor
{
public:
  emNISTPhysicsConstructor()
  : G4VPhysicsConstructor("emNISTPhysics",7777) {}  // Physics type 7777????
  ~emNISTPhysicsConstructor() {};
  virtual void ConstructParticle() {}  // Nothing to be done
  virtual void ConstructProcess();
};

#endif
