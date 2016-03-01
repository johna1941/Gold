#include "emNISTElasticScatterModel.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4NistManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

#include "G4DataSet.hh"
#include "G4DataVector.hh"

#include "G4CrossSectionDataSet.hh"
#include "G4LinInterpolation.hh"

#include "G4ParticleChangeForGamma.hh"

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

G4ThreadLocal G4ParticleChangeForGamma* NISTparticleChange = nullptr;

G4CrossSectionDataSet * angularDistributionDataset;
G4CrossSectionDataSet * transportDataset;

G4bool isInitialisedNIST = false;

G4NistManager* manNist = G4NistManager::Instance();

emNISTElasticScatterModel::emNISTElasticScatterModel(const G4ParticleDefinition*,
                                                     const G4String& nam)
  :G4VEmModel(nam)
{
  if (NISTparticleChange == nullptr) NISTparticleChange = new G4ParticleChangeForGamma;
}

emNISTElasticScatterModel::~emNISTElasticScatterModel()
{
  delete NISTparticleChange;
  G4cout << "NISTElasticScatter" << G4endl;
}

void emNISTElasticScatterModel::Initialise(const G4ParticleDefinition* particle,
                                           const G4DataVector& /*cuts*/)
{
  if (isInitialisedNIST) { return; }

  const G4ParticleDefinition* particleLoc = particle;
  
  //////////////////////////////////////////////////
  // Load cross section tables
  //////////////////////////////////////////////////
  
  G4LinInterpolation * algo = new G4LinInterpolation;

  double bohrRadius = 0.0529177211*nm;
  angularDistributionDataset = new G4CrossSectionDataSet(algo->Clone(),
							 eV,
							 pow(bohrRadius,2));
  angularDistributionDataset->LoadData("NewCompliedTable");

  transportDataset = new G4CrossSectionDataSet(algo->Clone(), eV, m);
  transportDataset->LoadData("new64_AU_EMFP");

  //////////////////////////////////////////////////
  // Clean up
  ////////////////////////////////////////////////// 

  isInitialisedNIST = true;

  delete algo;
  
  G4cout << "emNISTElasticScatterModel::Initialise" << G4endl;
}

G4double emNISTElasticScatterModel::
CrossSectionPerVolume(const G4Material* material,
		      const G4ParticleDefinition* particleDefinition,
		      G4double ekin,
		      G4double emin,
		      G4double emax)
{
  //G4Material* gold = manNist->FindOrBuildMaterial("G4_Au");
  G4double ekinLoc = ekin;
  G4double particleEnergy = ekinLoc;
  G4double lowerLimit = DBL_MIN;
  
  //if (material != gold) {return lowerLimit;}
  //else if (particleEnergy >= 1550*eV) {return lowerLimit;}
  if (particleEnergy >= 1550*eV) {return lowerLimit;}
  else if (particleEnergy <= 50*eV) {return lowerLimit;}
  
  G4double meanFreePath = transportDataset->FindValue(particleEnergy,0);

  return (2000/(meanFreePath));
}

void emNISTElasticScatterModel::
SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
		  const G4MaterialCutsCouple* couple,
		  const G4DynamicParticle* particle, G4double, G4double) {
  //G4cout << "NIST elastic scattering called" << G4endl;
  const G4DynamicParticle* particleLoc = particle;
  
  //G4ThreadLocal G4ParticleChangeForGamma NISTparticleChange;// = new G4ParticleChangeForGamma;
  NISTparticleChange = GetParticleChangeForGamma();
  
  G4double particleEnergy = particleLoc->GetKineticEnergy();
  
  // Choose scattering angles
  G4double thetaDegrees = RandomTheta(particleEnergy);
  G4double phiDegrees = (G4UniformRand()*360)-180;

  G4double theta = thetaDegrees*degree;//(thetaDegrees/180)*pi;
  //  G4double phi = -pi+(phiDegrees/180)*pi;
  G4double phi = phiDegrees*degree;//(phiDegrees/180)*pi;

  //G4double theta = pi/2;//RandomTheta(particleEnergy);
  //G4double phi = pi/2;//G4UniformRand()*360;
  
  G4ThreeVector currentMomentumDirection = particleLoc->GetMomentumDirection();
  G4ThreeVector proposedMomentumDirection = currentMomentumDirection;

  /*  G4ThreeVector pmdZ = particleLoc->GetMomentumDirection();
  G4ThreeVector pmdX = pmdZ.orthogonal();
  G4ThreeVector pmdY = pmdZ.cross(pmdX);

  G4double xDir = std::sqrt(1. - pow(cos(theta),2));
  G4double yDir = xDir;
  xDir *= std::cos(phi);
  yDir *= std::sin(phi);

  G4ThreeVector proposedMomentumDirection(xDir*pmdX + yDir*pmdY + std::cos(theta)*pmdZ);*/

  G4double currentTheta = currentMomentumDirection.theta();
  G4double currentPhi = currentMomentumDirection.phi();
  G4double newTheta = currentTheta + theta;
  G4double newPhi = currentPhi + phi;

  if (newTheta>pi) newTheta -= pi;
  if (newPhi>(2*pi)) newPhi -= 2*pi;
  proposedMomentumDirection.setTheta(newTheta);
  proposedMomentumDirection.setPhi(phi);//newPhi);

  /*G4cout << "angle in: Theta " << currentMomentumDirection.theta()/degree
	 << " Phi " << currentMomentumDirection.phi()/degree << G4endl;
  G4cout << "angle out: Theta " << proposedMomentumDirection.theta()/degree
	 << " Phi " << proposedMomentumDirection.phi()/degree << G4endl;
  G4cout << "angle change: Theta "
	 << (proposedMomentumDirection.theta()/degree)
	 -(currentMomentumDirection.theta()/degree) << G4endl;*/

  NISTparticleChange->ProposeMomentumDirection(proposedMomentumDirection.unit());
  
  //if (!NISTparticleChange->CheckIt(*NISTparticleChange->GetCurrentTrack())) {
  //G4cout << "oh dear NIST " << G4endl;
  //  abort();
  //}
}

G4double emNISTElasticScatterModel::RandomTheta(G4double E)
{
  G4double particleEnergy = E;
  
  //////////////////////////////////////////////////
  // Obtain pdf at particle energy
  //////////////////////////////////////////////////

  // Find total across angles
  G4double totalCrossSection = angularDistributionDataset->FindValue(particleEnergy/eV,0);

  // Find each value as a fraction of the total
  G4DataVector pdf;
  //pdf.clear(); 
  for (G4int i = 0; i!=180; ++i) {
      G4double value = angularDistributionDataset->GetComponent(i)
	                                         ->FindValue(particleEnergy/eV,0);
      pdf.push_back(value/totalCrossSection);
    }

  //////////////////////////////////////////////////
  // Obtain random angle
  //////////////////////////////////////////////////

  G4double pdfMax = *max_element(pdf.begin(),pdf.end()); // find pdf maximum
  G4double randAngle = floor(G4UniformRand()*180); // randomly select angle
  G4double randY = G4UniformRand()*pdfMax; // Generate random number between 0 and pdfMax
  while(pdf[randAngle]<randY) {
    randAngle = floor(G4UniformRand()*180); 
    randY = G4UniformRand()*pdfMax;
    //cout << "Sample reject" << endl;
  }
  //cout << "Sample accept" << endl;
  G4double theta = randAngle;
  
  return theta;
}
