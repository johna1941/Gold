// *******************************************************************
// * License and Disclaimer                                          *
// *                                                                 *
// * This software is copyright of Geant4 Associates International   *
// * Ltd (hereafter 'G4AI'). It is provided under the terms and      *
// * conditions described in the file 'LICENSE' included in the      *
// * software system.                                                *
// * Neither the authors of this software system nor G4AI make any   *
// * representation or warranty, express or implied, regarding this  *
// * software system or assume any liability for its use.            *
// * Please see the file 'LICENSE' for full disclaimer and the       *
// * limitation of liability.                                        *
// *******************************************************************
// $Id$
// John Allison  28th August 2015

#ifndef Gold1SENSITIVEDETECTOR_HH
#define Gold1SENSITIVEDETECTOR_HH

#include "G4VSensitiveDetector.hh"

class Gold1SensitiveDetector : public G4VSensitiveDetector
{
public:
  Gold1SensitiveDetector(const G4String& name);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
};

#endif
