#include "OpticalInelasticScatterModel.hh"

#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <iomanip>

#include <stdlib.h> 
#include <fstream>

#include <cmath>
#include <array>
#include <complex>

#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define m_e 1
#define pi 3.14159265

#define au_length 5.29177E-11
#define au_velocity 2.18769E6
#define au_energy 4.35974E-18
#define si_energy 1.60217E-19

using namespace std;

// function declarations
double calcFirstBorn(double E, double omega);
double calcSecondBorn(double E, double omega);
double calcExchange(double E, double omega);
double calcELF (double k, double omega);
complex<double> calcDielectricFunction(double k, double omega);

// Oscillator strength, damping coefficient and critial point energy
//const array <double,12> A_i = {79, 9, 36, 17, 60, 100, 120, 155, 145, 280, 360, 183};
//const array <double,12> gamma_i = {0.1, 1, 1.9, 2.3, 4, 9, 10, 6, 7.2, 20, 28, 26};
//const array <double,12> omega_i = {0, 3.1, 4.1, 5.3, 8.17, 12, 14, 21.3, 29.5, 38.5, 63, 100};

const array <double,12> A_i = {0.1066898493, 0.0121545398, 0.0486181592, 0.0229585752, 0.0810302653, 0.1350504422, 0.1620605306, 0.2093281854, 0.1958231412, 0.3781412382, 0.4861815919, 0.2471423092};
const array <double,12> gamma_i = {0.003674921, 0.0367492098, 0.0698234986, 0.0845231826, 0.1469968393, 0.3307428883, 0.3674920982, 0.2204952589, 0.2645943107, 0.7349841963, 1.0289778748, 0.9554794552};
const array <double,12> omega_i = {0, 0.1139225504, 0.1506717602, 0.194770812, 0.3002410442, 0.4409905178, 0.5144889374, 0.7827581691, 1.0841016896, 1.4148445779, 2.3152002184, 3.6749209815};

// global variables
//double E = 0; // kinetic energy of subject electron
//double v = 0; // scalar velocity of subject electron

//double k = 0; // momentum change of subject electron
//double kMax = 0; // maximum momentum change of subject electron
//double kMin = 0; // minimum momentum change of subject electron
//double omega = 0; // energy loss of subject electron through inelastic scattering

//complex<double> epsilon (0,0); // dielectric function
const double epsilon_b = 1; // background dielectric constant
//double epsilon_1 = 0; // real part of dielectric function
//double epsilon_2 = 0; // imaginary part of dielectric function

const int samplingPoints = 1000;
//array <double,samplingPoints> ELF; // the energy loss function
//array <double,samplingPoints> OELF; // the optical energy loss function

//double DIMFP = 0; // the differential inverse mean free path
//double DIMFPFirstBorn = 0; // first Born approximation contribution to DIMFP
//double DIMFPSecondBorn = 0; // second Born approximatin correction to DIMFP
//double DIMFPExchange = 0; // exchange-corrected contribution to DIMFP

//double IMFP = 0; // the inverse mean free path
//double MFP = 0; // the mean free path

// second born approximation parameters
//double chi = 0; 
//double I = 0;
//double L = 0;

OpticalInelasticScatterModel::OpticalInelasticScatterModel() {
  // convert dielectric function parameters to atomic units
  /*for(int i = 0 ; i < omega_i.size() ; ++i) {
    A_i[i] = (A_i[i]*pow(si_energy,2))/pow(au_energy,2);
    gamma_i[i] = (gamma_i[i]*si_energy)/au_energy;
    omega_i[i] = (omega_i[i]*si_energy)/au_energy;
    }*/
}

OpticalInelasticScatterModel::~OpticalInelasticScatterModel() {
  G4cout << "OpticalInelasticScatterModelCalc" << G4endl;
}

double OpticalInelasticScatterModel::computeMFP(double energy) {
  double DIMFP = 0; // the differential inverse mean free path
  double DIMFPFirstBorn = 0; // first Born approximation contribution to DIMFP
  double DIMFPSecondBorn = 0; // second Born approximatin correction to DIMFP
  double DIMFPExchange = 0; // exchange-corrected contribution to DIMFP

  double IMFP = 0; // the inverse mean free path
  double MFP = 0; // the mean free path
  
  //////////////////////////////////////////////////
  // Perform calculation initialisation
  //////////////////////////////////////////////////
  
  double E = (energy*si_energy)/au_energy; // in eV
  double Emin = (35*si_energy)/au_energy;

  //////////////////////////////////////////////////
  // Begin calculations - interating over omega
  //////////////////////////////////////////////////

  double dOmega = E/samplingPoints;
  for (double omega = dOmega ; omega < E ; omega += dOmega) {
  //if (E > Emin) {
  //for (double omega = Emin ; omega < E ; omega += dOmega) {
      //// Step 1: Calculate DIMFP for first born approximation
      DIMFPFirstBorn = calcFirstBorn(E,omega);
      
      //// Step 2: Calculate DIMFP correction for second born approximation
      DIMFPSecondBorn = calcSecondBorn(E,omega);
      
      //// Step 3: Calculate DIMFPexchange correction term
      DIMFPExchange = calcExchange(E,omega);
      
      //// Step 4: Combine DIMFP contributions
      DIMFP = DIMFPFirstBorn + DIMFPSecondBorn;// + DIMFPExchange
      //- sqrt((DIMFPFirstBorn + DIMFPSecondBorn)*DIMFPExchange);
      
      /*G4cout << "DIdebug "
	     << E << " "
	     << Emin << " "
	     << omega << "   "
	<< DIMFPFirstBorn << "   "
	<< DIMFPSecondBorn << "   "
	<< DIMFPExchange << "   "
	<< DIMFP << G4endl;*/

      IMFP += DIMFP*dOmega;
  }

  //////////////////////////////////////////////////
  // Convert IMFP to MFP (in SI) and output
  //////////////////////////////////////////////////
   
  MFP = (1/IMFP)*au_length;

  return MFP;
  //} else return DBL_MAX;
}

double OpticalInelasticScatterModel::sampleELF(double E) {
  //////////////////////////////////////////////////
  // Perform calculation initialisation
  //////////////////////////////////////////////////
  
  double E_loc = E*si_energy/au_energy; 

  //////////////////////////////////////////////////
  // Calculate ELF
  //////////////////////////////////////////////////
  
  double v_loc = sqrt(2*E_loc);

  // integration over momentum transfer
  double dOmega = E_loc/samplingPoints;
  double ELFX[samplingPoints-1];
  double ELFY[samplingPoints-1];
  
  double ELFmax = 0;
  int i = 0;

  for (double omega_loc = dOmega ; omega_loc < E_loc ; omega_loc += dOmega) {
    double kMin_loc = sqrt(2*E_loc)-sqrt(2*(E_loc-omega_loc));
    double kMax_loc = sqrt(2*E_loc)+sqrt(2*(E_loc-omega_loc));
    double dK = (kMax_loc-kMin_loc)/samplingPoints;

    double ELFpoint = 0;
    for (double k = kMin_loc ; k <= kMax_loc ; k += dK) {
      ELFpoint += calcELF(k,omega_loc)*dK;
    }
    if (ELFpoint > ELFmax) ELFmax = ELFpoint;
    ELFX[i] = omega_loc;
    ELFY[i] = ELFpoint;
    ++i;
  }

  //////////////////////////////////////////////////
  // Sample ELF
  //////////////////////////////////////////////////

  gsl_interp_accel *accelerator = gsl_interp_accel_alloc();
  gsl_spline *energyLossFunction = gsl_spline_alloc(gsl_interp_cspline, samplingPoints-1);
  
  gsl_spline_init (energyLossFunction, ELFX, ELFY, samplingPoints-1);

  double yValue = 0;
  double omegaOut = 0;
  do {
    double R1 = G4UniformRand();
    omegaOut = R1*(E_loc-dOmega-dOmega)+dOmega;
    yValue = gsl_spline_eval(energyLossFunction, omegaOut, accelerator);
  } while ((G4UniformRand()*ELFmax)>yValue);

  gsl_spline_free (energyLossFunction);
  gsl_interp_accel_free (accelerator);
  
  //  delete energyLossFunction;
  //  delete accelerator;
  
  return omegaOut*au_energy/si_energy;
}

double OpticalInelasticScatterModel::calcELFpoint(double E, double omega) {
  //////////////////////////////////////////////////
  // Perform calculation initialisation
  //////////////////////////////////////////////////
  
  double E_loc = E*si_energy/au_energy;
  double omega_loc = omega*si_energy/au_energy;

  //////////////////////////////////////////////////
  // Calculate ELF
  //////////////////////////////////////////////////
  
  double v_loc = sqrt(2*E_loc);

  // integration over momentum transfer

  double kMin_loc = sqrt(2*E_loc)-sqrt(2*(E_loc-omega_loc));
  double kMax_loc = sqrt(2*E_loc)+sqrt(2*(E_loc-omega_loc));
  double dK = (kMax_loc-kMin_loc)/samplingPoints;

  double ELFpoint = 0;
  if (kMin_loc<kMax_loc) {
    for (double k = kMin_loc ; k <= kMax_loc ; k += dK) {
      ELFpoint += calcELF(k,omega_loc)*dK;
    }
  } else {
    //ELFpoint = calcELF(kMin_loc,omega_loc);
  }
  //G4cout << ELFpoint << G4endl;      
  return ELFpoint;
}

double OpticalInelasticScatterModel::calcFirstBorn(double E, double omega) {
  double E_loc = E;
  double omega_loc = omega;
  
  // calculate constant (A)
  double v_loc = sqrt(2*E_loc);
  double A = 2/(pi*pow(v_loc,2));

  // integration over momentum transfer to obtain second term (B)
  double kMin_loc = sqrt(2*E_loc)-sqrt(2*(E_loc-omega_loc));
  double kMax_loc = sqrt(2*E_loc)+sqrt(2*(E_loc-omega_loc));

  double B = 0;
  double dK = (kMax_loc-kMin_loc)/samplingPoints;
  for (double k = kMin_loc ; k <= kMax_loc ; k += dK) {
    double ELF_loc = calcELF(k, omega_loc);
    B += (1/k)*ELF_loc*dK;
  }
  return A*B;  
}

double OpticalInelasticScatterModel::calcSecondBorn(double E, double omega) {
  double E_loc = E;
  double omega_loc = omega;
  
  // calculate constant (A)
  double v_loc = sqrt(2*E_loc);
  double A = -2/(pi*pow(v_loc,2));

  // retrieve OELF
  double OELF_loc = calcELF(0, omega_loc);
  
  // calculate second born approximation parameters (L)
  double chi_loc = 1/sqrt(2*omega_loc);
  double I_loc = -(3/2)*pi*log((chi_loc*omega_loc)/v_loc)-2.4;
  double L_loc = (omega_loc/pow(v_loc,3))*I_loc;

  /*cout << "2nd born "
       << chi_loc << " "
       << I_loc << " "
       << L_loc << "   " 
       << A*OELF_loc*L_loc << " " 
       << ((chi_loc * omega_loc) / v_loc)
       << G4endl;*/
  
  return A*OELF_loc*L_loc;
}

double OpticalInelasticScatterModel::calcExchange(double E, double omega) {
  double E_loc = E;
  double omega_loc = omega;

  double exchangeOmega = E_loc - omega_loc;

  /*  G4cout << "calcExchange " << E_loc
	 << " " << omega_loc
	 << " " << exchangeOmega
	 << G4endl;*/
  double DIMFPFirstBorn_loc = calcFirstBorn(E_loc,exchangeOmega);
  double DIMFPSecondBorn_loc = calcSecondBorn(E_loc,exchangeOmega);
  //cout << DIMFPFirstBorn_loc << " " << DIMFPSecondBorn_loc << endl;
  return DIMFPFirstBorn_loc+DIMFPSecondBorn_loc;
}

double OpticalInelasticScatterModel::calcELF (double k, double omega) {
  double k_loc = k;
  double omega_loc = omega;

  complex<double> epsilon = calcDielectricFunction(k_loc,omega_loc);
   complex<double> inverseEpsilon = -1./epsilon;

  return inverseEpsilon.imag();
}

complex<double> OpticalInelasticScatterModel::calcDielectricFunction (double k, double omega) {
  double k_loc = k;
  double omega_loc = omega;
  
  double epsilon_1_sum = 0;
  double epsilon_2_sum = 0;
  for (int i = 0 ; i < omega_i.size() ; ++i) {
    epsilon_1_sum += ((A_i[i]*(pow(omega_loc,2)-pow((omega_i[i]+(pow(k_loc,2)/2)),2)))/((pow(pow(omega_loc,2)-pow((omega_i[i]+(pow(k_loc,2)/2)),2),2)+(pow(omega_loc,2)*pow(gamma_i[i],2)))));

    epsilon_2_sum += ((A_i[i]*gamma_i[i]*omega_loc)/((pow(pow(omega_loc,2)-pow((omega_i[i]+(pow(k_loc,2)/2)),2),2)+(pow(omega_loc,2)*pow(gamma_i[i],2)))));
  }

  double epsilon_1 = epsilon_b - epsilon_1_sum;
  double epsilon_2 = epsilon_2_sum;

  return complex<double>(epsilon_1,epsilon_2);
}
