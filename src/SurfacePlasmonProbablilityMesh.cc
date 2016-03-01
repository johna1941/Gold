#include "SurfacePlasmonProbabilityMesh.hh"

#include <iostream>
#include <fstream>
//#include <iomanip>

//#include <stdlib.h> 

#include <cmath>
#include <array>
#include <complex>

#define m_e 1
#define pi 3.14159265
#define c 2.99792E8

#define au_length 5.29177E-11
#define au_velocity 2.18769E6
#define au_energy 4.35974E-18
#define si_energy 1.60217E-19

using namespace std;

// function declarations
complex<G4double> calcRc (G4double omega, G4double v, G4double s,
			complex<G4double> epsilon, G4double k_s);
complex<G4double> calcDielectricFunction(G4double k, G4double omega);
G4double calcPs (G4double k_s, G4double v, G4double alpha, G4double k,
	       complex<G4double> epsilon, complex<G4double> Rc,
	       array<G4double,3> v_vec, array<G4double,3> k_vec,
	       G4double omega);
  
// local variables
array <G4double,3> v_vec; // incoming particle velocity
G4double s; // sample thickness (distance from the e- to the surface??)
G4double alpha ; // incoming particle angle to surface normal

array <G4double,3> k_0; // incoming partile momentum
G4double E; // incoming particle energy
G4double v; // scalar incoming velocity

G4double omega; // energy loss
array <G4double,3> k_vec; // vector momentum loss
G4double k; // scalar momentum loss
G4double k_s; // parallel component of k along surface
array <G4double,3> k_1; // momentum vector of electron after scatter
G4double k_1_mag; // magnitude of post-scatter electron vector
 G4double k_0_mag; // magnitude of pre-scatter electron vector

G4double theta; // direction of scattered electron
G4double phi; // direction of scattered electron
G4double theta_k; // direction of momentum loss vector
G4double phi_k; // direction of momentum loss vector

complex<G4double> Rc; // sample width dependant term
complex<G4double> epsilon; // dielectric constant

G4double epsilon_1; // real part of dielectric function
G4double epsilon_2; // imaginary part of dielectric function

G4double Ps; // probability of a surface excitation occuring

const G4double epsilon_b = 1; // background dielectric constant

// Oscillator strength, damping coefficient and critial point energy
//array <G4double,12> A_i = {79, 9, 36, 17, 60, 100, 120, 155, 145, 280, 360, 183};
//array <G4double,12> gamma_i = {0.1, 1, 1.9, 2.3, 4, 9, 10, 6, 7.2, 20, 28, 26};
//array <G4double,12> omega_i = {0.1, 3.1, 4.1, 5.3, 8.17, 12, 14, 21.3, 29.5, 38.5, 63, 100};
const array <G4double,12> A_i = {0.10669,0.01215,0.04862,0.02296,0.08103,0.13505,0.16206,0.20933,0.19582,0.37814,0.48618,0.24714};
const array <G4double,12> gamma_i = {0.0037,0.0367,0.0698,0.0845,0.1470,0.3307,0.3675,0.2205,0.2646,0.7350,1.0290,0.9555};
const array <G4double,12> omega_i = {0.00367,0.11392,0.15067,0.19477,0.30024,0.44099,0.51449,0.78276,1.08410,1.41484,2.31520,3.67492};

const array <G4double,3> surface_vec = {0,1,0}; // unit vector of surface

G4double omegaOut, thetaOut, phiOut;

SurfacePlasmonProbabilityMesh::SurfacePlasmonProbabilityMesh (G4double E_eV,
							      G4double s_m,
							      G4double alpha_in) {
  // Convert to atomic units
  E = (E_eV*si_energy)/au_energy;
  s = s_m/au_length;

  /*for(int i = 0 ; i < sizeof(omega_i) ; ++i) {
    // first to SI
    A_i[i] *= pow(si_energy,2);
    gamma_i[i] *= si_energy;
    omega_i[i] *= si_energy;
    
    A_i[i] /= pow(au_energy,2);
    gamma_i[i] /= au_energy;
    omega_i[i] /= au_energy;
  }*/

  // Calculations from initial conditions
  v = (c*sqrt(E)*sqrt(2*pow(c,2)*m_e+E))/(pow(c,2)*m_e+E);
  v_vec[0] = v*sin(alpha);
  v_vec[1] = v*cos(alpha);
  v_vec[2] = 0;

  G4double lorentzFactor = 1./sqrt(1-(pow(v,2)/pow(c/au_velocity,2)));
  
  k_0[0] = m_e*lorentzFactor*v_vec[0];
  k_0[1] = m_e*lorentzFactor*v_vec[1];
  k_0[2] = m_e*lorentzFactor*v_vec[2];

  k_0_mag = sqrt(pow(k_0[0],2)+pow(k_0[1],2)+pow(k_0[2],2));

  do {
    // Define random variables
    G4double R1 = G4UniformRand();
    G4double R2 = G4UniformRand();
    G4double R3 = G4UniformRand();
    
    //// Step 1: Generate theta_k, phi_k, k
    k = v*pow(R1,(1./3.));
    theta_k = acos(1-(2*R2));
    phi_k = 2*pi*R3;
    
    //// Step 2: Calculate k_vec, omega, k_s
    k_vec[0] = k*sin(theta_k)*cos(phi_k);
    k_vec[1] = k*sin(theta_k)*sin(phi_k);
    k_vec[2] = k*cos(theta_k);
    
    G4double dotVvecKvec = (v_vec[0]*k_vec[0]) +
      (v_vec[1]*k_vec[1]) +
      (v_vec[2]*k_vec[2]);
    
    omega = dotVvecKvec-((pow(k,2))/2);
    if (omega<0 || omega>E) continue;

    // k_s = surface unit vector dot k
    k_s = (surface_vec[0]*k_vec[0]) +
      (surface_vec[1]*k_vec[1]) +
      (surface_vec[2]*k_vec[2]);
    
    //// Step 3: Translate to global frame
    k_1[0] = k_0[0]-k_vec[0];
    k_1[1] = k_0[1]-k_vec[1];
    k_1[2] = k_0[2]-k_vec[2];
    
    //// Step 4: Calculate abs(k_1), theta, phi
    k_1_mag = sqrt(pow(k_1[0],2)+pow(k_1[1],2)+pow(k_1[2],2));
     G4double dotKK_0 = (k_0[0]*k_vec[0]) +
      (k_0[1]*k_vec[1]) +
      (k_0[2]*k_vec[2]);
    G4double k_mag = sqrt(pow(k_vec[0],2)+pow(k_vec[1],2)+pow(k_vec[2],2));
    if(dotKK_0/(k_mag*k_0_mag)>1||dotKK_0/(k_mag*k_0_mag)<-1) continue; 
    theta = acos(dotKK_0/(k_mag*k_0_mag));//acos(k_1[2]/k_1_mag);
    phi = atan2(k_1[1],k_1[0]);
    
    //// Step 5: Compute Ps
    epsilon = calcDielectricFunction(k, omega);
    //Rc = calcRc(omega, v, s, epsilon, k_s);
    Rc = 0;
    
    complex<G4double> test = (epsilon*(1.+epsilon));
    if (test.real()==0||test.imag()==0) continue;
    
    //Ps = calcPs(k_s,v,alpha,k,epsilon,Rc,v_vec,k_vec,omega); 
    Ps = calcPs(k_s,v,alpha,k,epsilon,Rc,v_vec,k_vec,omega);

    if (Ps<0 || omega<0) continue;
    //G4cout << omega << " " << Ps*1000 << G4endl;
    //G4cout << "Ps: " << Ps << G4endl;
    //if (G4UniformRand()<Ps) break;
  } while (G4UniformRand()>Ps);

  // might be getting interferance from other threads

  omegaOut = omega;
  thetaOut = theta;
  phiOut = phi;
  return;
}

SurfacePlasmonProbabilityMesh::~SurfacePlasmonProbabilityMesh ()
{
  //G4cout << "SurfacePlasmonProbabiltyMesh" << G4endl;
}

complex<G4double> SurfacePlasmonProbabilityMesh::calcRc(G4double omega_loc,
						      G4double v_loc,
						      G4double s_loc,
						      complex<G4double> epsilon_loc,
						      G4double k_s_loc) {
  return (pow(sin((s_loc*omega_loc)/(2*v_loc)),2)/
	  (epsilon_loc+tanh((s_loc*k_s_loc)/2))) +
    (pow(cos((s_loc*omega_loc)/(2*v_loc)),2)/
     (epsilon_loc+(1/tanh((s_loc*k_s_loc)/2))));
}

complex<G4double> SurfacePlasmonProbabilityMesh::calcDielectricFunction (G4double k_loc, G4double omega_loc) {
  G4double epsilon_1_sum = 0;
  G4double epsilon_2_sum = 0;
  for (int i = 0 ; i < omega_i.size() ; ++i) {
    if(((pow(pow(omega_loc,2)-pow((omega_i[i]+
				   (pow(k_loc,2)/2)),2),2)+
	 (pow(omega_loc,2)*pow(gamma_i[i],2))))==0) continue;
    
    epsilon_1_sum += ((A_i[i]*(pow(omega_loc,2)-
			       pow((omega_i[i]+(pow(k_loc,2)/2)),2)))
		      /(pow(pow(omega_loc,2)-
			    pow((omega_i[i]+(pow(k_loc,2)/2)),2),2)
			+(pow(omega_loc,2)*pow(gamma_i[i],2))));
      
    epsilon_2_sum += (A_i[i]*gamma_i[i]*omega_loc)
      /(pow(pow(omega_loc,2)-pow((omega_i[i]+(pow(k_loc,2)/2)),2),2)
	+(pow(omega_loc,2)*pow(gamma_i[i],2)));


    ////////////////////// EPSILON 2 IS LIKELY WRONG
    // BRACKETING GETS A POWER 2 IN THE WRONG PLACE
  }
  
  epsilon_1 = epsilon_b - epsilon_1_sum;
  epsilon_2 = epsilon_2_sum;

  return complex<G4double>(epsilon_1,epsilon_2);
}

G4double SurfacePlasmonProbabilityMesh::calcPs (G4double k_s_loc, G4double v_loc,
					      G4double alpha_loc, G4double k_loc,
					      complex<G4double> epsilon_loc,
					      complex<G4double> Rc_loc,
					      array<G4double,3> v_vec_loc,
					      array<G4double,3> k_vec_loc,
					      G4double omega_loc) {
   G4double dotVvecKvec_loc = (v_vec_loc[0]*k_vec_loc[0]) +
                        (v_vec_loc[1]*k_vec_loc[1]) +
                        (v_vec_loc[2]*k_vec_loc[2]);
   
   G4double Ps_term1 = (2*abs(k_s_loc))/
     (pow(pi,2)*v_vec_loc[1]*pow(k_loc,4));
   
   complex<double> Ps_term2_comp = ((pow((epsilon_loc-1.),2))/
				    (epsilon_loc*(1.+epsilon_loc)));
   G4double Ps_term2 = Ps_term2_comp.imag();

   // all the negative ps have a negative term 2 
   /*if (Ps_term2<0
       && omega>=0 && omega<=E // check omega is physical
       && omega<((pow(k,2))/2)
       ) cout  << std::imag(pow((epsilon_loc-1.),2)) << ","
	       << omega << ","
	       << pow((k_loc/2),2)
	       << endl;*/
   
   return (Ps_term1 * Ps_term2);
}

G4double SurfacePlasmonProbabilityMesh::getOmega() {
  return omegaOut;
}

G4double SurfacePlasmonProbabilityMesh::getTheta() {
  return thetaOut;
}

G4double SurfacePlasmonProbabilityMesh::getPhi() {
  return phiOut;
}
