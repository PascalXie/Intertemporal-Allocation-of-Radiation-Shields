//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  //fParticleGun->SetParticleEnergy(185.7*keV);

  // reading
  // Atmosphere or Ground
  ReadIsotope("../data/129I.txt",  Istp_129I_energies,  Istp_129I_weights);
  ReadIsotope("../data/134Cs.txt", Istp_134Cs_energies, Istp_134Cs_weights);
  ReadIsotope("../data/137Cs.txt", Istp_137Cs_energies, Istp_137Cs_weights);
  // Well
  ReadIsotope("../data/131I.txt",  Istp_energies_W[0],  Istp_weights__W[0]);
  ReadIsotope("../data/134Cs.txt", Istp_energies_W[1],  Istp_weights__W[1]);
  ReadIsotope("../data/137Cs.txt", Istp_energies_W[2],  Istp_weights__W[2]);
  ReadIsotope("../data/140Ba.txt", Istp_energies_W[3],  Istp_weights__W[3]);
  ReadIsotope("../data/140La.txt", Istp_energies_W[4],  Istp_weights__W[4]);
  ReadIsotope("../data/89Sr.txt",  Istp_energies_W[5],  Istp_weights__W[5]);
  ReadIsotope("../data/90Sr.txt",  Istp_energies_W[6],  Istp_weights__W[6]);

  double Area_working = PI_*25.*25.; // m2
  double Area_well = PI_*11.*11.; // m2
  double Area_Ground = Area_working-Area_well;
  G4cout<<"Area_working "<<Area_working<<endl;
  G4cout<<"Area_well "<<Area_well<<endl;
  G4cout<<"Area_Ground "<<Area_Ground<<endl;

  // set activities
  // Ground
  Act_grd[0] = Area_Ground*71.7e1;   // MBq, 129I
  Act_grd[1] = Area_Ground*29.0e7;    // MBq, 134Cs
  Act_grd[2] = Area_Ground*28.67e7;   // MBq, 137Cs
  Act_grd_Total = (Act_grd[0] + Act_grd[1] + Act_grd[2]);

  G4cout<<"Act_grd[0] "<<Act_grd[0]<<endl;
  G4cout<<"Act_grd[1] "<<Act_grd[1]<<endl;
  G4cout<<"Act_grd[2] "<<Act_grd[2]<<endl;


  // Atmosphere
  double scale_Act = 1.;
  Act_atm[0] = scale_Act * Act_grd[0]; // MBq, 129I
  Act_atm[1] = scale_Act * Act_grd[1]; // MBq, 134Cs
  Act_atm[2] = scale_Act * Act_grd[2]; // MBq, 137Cs
  Act_atm_Total = Act_atm[0] + Act_atm[1] + Act_atm[2];

  // Well
  Act_wel[0] = 1.36e9;  // MBq, 131I 
  Act_wel[1] = 3.16e8;  // MBq, 134Cs
  Act_wel[2] = 4.21e8;  // MBq, 137Cs
  Act_wel[3] = 8.93e6;  // MBq, 140Ba
  Act_wel[4] = 6.52e11; // MBq, 140La
  Act_wel[5] = 2.36e5;  // MBq, 89Sr
  Act_wel[6] = 5.53e4;  // MBq, 90Sr

  //Act_wel[0] = 0;  // MBq, 131I 
  //Act_wel[1] = 0;  // MBq, 134Cs
  //Act_wel[2] = 0;  // MBq, 137Cs
  //Act_wel[3] = 0;  // MBq, 140Ba
  //Act_wel[4] = 0; // MBq, 140La
  //Act_wel[5] = 0;  // MBq, 89Sr
  //Act_wel[6] = 0;  // MBq, 90Sr
  for(int i=0;i<7;i++)
  {
    Act_wel_Total += Act_wel[i];
  }


  double Act_Total = Act_atm_Total + Act_grd_Total + Act_wel_Total;
  G4cout<<"Act_atm_Total "<<Act_atm_Total<<endl;
  G4cout<<"Act_grd_Total "<<Act_grd_Total<<endl;
  G4cout<<"Act_wel_Total "<<Act_wel_Total<<endl;
  G4cout<<"Act_Total "<<Act_Total<<endl;

  double N_sim = 50000000;
  double time_real = 4.*N_sim/(Act_Total*1e6); // Act_total : Bq
  G4cout<<"time_real "<<time_real<<endl;

  double mass = 0.5*0.16*0.04 * 2328.3; // kg
  G4cout<<"mass "<<mass<<" kg"<<endl;

  // Geometry
  radius_W   = 11.*m;
  radius_A_G = 25.*m;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //
  // step 1 : choose an environment 
  int EnvrID = GetEnvironmentFrom3();
  //G4cout<<"EnvrID "<<EnvrID<<G4endl;

  // step 2 : choose an isotope
  int IstpID = 0;
  if(EnvrID==0 || EnvrID==1)
    IstpID = GetIsotope_G_or_A();
  else
    IstpID = GetIsotope_W();

  // step 3 : choose a location
  G4double x0 = 0.; 
  G4double y0 = 0.; 
  G4double z0 = 0.; 

  if(EnvrID==0) // ground, uniform
  {
      GetLocation_Ground(x0,y0,z0);
      //cout<<"ground, uniform"<<endl;
      //cout<<"Location "<<x0/m<<" m, "<<y0/m<<" m, "<<z0/m<<" m, "<<endl;
  }
  else if(EnvrID==1)// atmosphere, Gaussian
  {
      GetLocation_Atmosphere(x0,y0,z0);
      //cout<<"atmosphere, Gaussian"<<endl;
      //cout<<"Location "<<x0/m<<" m, "<<y0/m<<" m, "<<z0/m<<" m, "<<endl;
  }
  else if(EnvrID==2)// Well, uniform
  {
      GetLocation_Well(x0,y0,z0);
  }
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  // step 4 : choose a direction
  double dir[3];
  GetDirection(dir[0],dir[1],dir[2]);
  // Debug
  //G4cout<<"Direction : "<<dir[0]<<", "<<dir[1]<<",  "<<dir[2]<<G4endl;
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dir[0],dir[1],dir[2]));

  // step 5 : choose an energy
  double Energy = 0.;
  if(EnvrID==0) // ground, uniform
  {
      Energy = GetEnergy_A_G(IstpID);
  }
  else if(EnvrID==1) // atmosphere, Gaussian
  {
      Energy = GetEnergy_A_G(IstpID);
  }
  else if(EnvrID==2)// Well, uniform
  {
      Energy = GetEnergy_W(IstpID);
  }
  fParticleGun->SetParticleEnergy(Energy);

  // step 6 : generate an event
  fParticleGun->GeneratePrimaryVertex(anEvent);


  //
  // well
  //

  //
  // ground and atmosphere
  //
//  // step 1 : choose an isotope
//  int IstpID = GetIsotope_G_or_A();
//  //cout<<"IstpID "<<IstpID<<endl;
//
//  // step 2 : choose an environment 
//  int EnvrID = GetEnvironment();
//  //cout<<"EnvrID "<<EnvrID<<endl;
//
//  // step 3 : choose a location
//  G4double x0 = 0.; 
//  G4double y0 = 0.; 
//  G4double z0 = 0.; 
//
//  if(EnvrID==0) // ground, uniform
//  {
//      GetLocation_Ground(x0,y0,z0);
//      //cout<<"ground, uniform"<<endl;
//      //cout<<"Location "<<x0/m<<" m, "<<y0/m<<" m, "<<z0/m<<" m, "<<endl;
//  }
//  else // atmosphere, Gaussian
//  {
//      GetLocation_Atmosphere(x0,y0,z0);
//      //cout<<"atmosphere, Gaussian"<<endl;
//      //cout<<"Location "<<x0/m<<" m, "<<y0/m<<" m, "<<z0/m<<" m, "<<endl;
//  }
//  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
//
//  // step 4 : choose a direction
//  double dir[3];
//  GetDirection(dir[0],dir[1],dir[2]);
//  // Debug
//  //G4cout<<"Direction : "<<dir[0]<<", "<<dir[1]<<",  "<<dir[2]<<G4endl;
//  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dir[0],dir[1],dir[2]));
//
//  // step 5 : choose an energy
//  fParticleGun->SetParticleEnergy(GetEnergy(IstpID));
//
//  // step 6 : generate an event
//  fParticleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int B1PrimaryGeneratorAction:: GetEnvironmentFrom3() 
{
    double total = Act_atm_Total + Act_grd_Total + Act_wel_Total;
    double uniAct[3]; // the uniformed activity
    uniAct[0] = Act_atm_Total/total; // atmosphere
    uniAct[1] = Act_grd_Total/total; // ground
    uniAct[2] = Act_wel_Total/total; // well

    //G4cout<<"uniAct[0] "<<uniAct[0]<<G4endl;
    //G4cout<<"uniAct[1] "<<uniAct[1]<<G4endl;
    //G4cout<<"uniAct[2] "<<uniAct[2]<<G4endl;

    int EnID = 0;
    double r1 = G4UniformRand();
    if(r1<uniAct[0])                        EnID=0; // atmosphere
    else if(r1>=uniAct[0] && r1<uniAct[0]+uniAct[1])  EnID=1; // ground
    else                                    EnID=2; // well

    return EnID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void B1PrimaryGeneratorAction:: GetLocation_Well(double &x, double &y, double &z) 
{
    int isLocationGood = false;
    while(!isLocationGood)
    {
        double r1 = G4UniformRand();
        x = radius_W*2.*(G4UniformRand()-0.5);
        y = radius_W*2.*(G4UniformRand()-0.5);
        z = -0.5*G4UniformRand() *m;

        if(x<-2.*m || y<-2.*m) continue;

        double d2 = x*x + y*y;
        double d  = sqrt(d2);
        if(d<radius_W) break;
    }

    //cout<<"x, y, z: "<<x/m<<", "<<y/m<<", "<<z/m<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void B1PrimaryGeneratorAction:: GetLocation_Atmosphere(double &x, double &y, double &z) 
{

  int isLocationGood = false;
  while(!isLocationGood)
  {
    double SdtDev = 8.*m;
    x = G4RandGauss::shoot(0.,SdtDev);
    y = G4RandGauss::shoot(0.,SdtDev);

    double r1 = G4UniformRand();
    z = 0.5*G4UniformRand() *m;

    if(x>-2.*m && y>-2.*m) break;

    //cout<<"Atm: x, y, z: "<<x/m<<", "<<y/m<<", "<<z/m<<endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void B1PrimaryGeneratorAction:: GetLocation_Ground(double &x, double &y, double &z) 
{
    int isLocationGood = false;
    while(!isLocationGood)
    {
        double r1 = G4UniformRand();
        x = radius_A_G*2.*(G4UniformRand()-0.5);
        y = radius_A_G*2.*(G4UniformRand()-0.5);
        z = -0.5*G4UniformRand() *m;

        if(x<-2.*m || y<-2.*m) continue;

        double d2 = x*x + y*y;
        double d  = sqrt(d2);
        if(d<radius_A_G) break;
    }

    //cout<<"x, y, z: "<<x/m<<", "<<y/m<<", "<<z/m<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void B1PrimaryGeneratorAction:: GetDirection(double &x, double &y, double &z) 
{
    //double phi        = (G4UniformRand()/2.-0.5)*M_PI; // fly to the detector directly
    double phi      = G4UniformRand()*M_PI; 
    double theta    = G4UniformRand()*2.*M_PI;

    x = sin(phi) * cos(theta);
    y = sin(phi) * sin(theta);
	z = cos(phi);

	// Debug
	//G4cout<<"GetDirection : "<<x<<", "<<y<<",  "<<z<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int B1PrimaryGeneratorAction::GetIsotope_W()
{
    int ID = 4;

    // step 1
    double act[7];
    double total = 0;
    for(int i=0;i<7;i++)
    {
        act[i] = Act_wel[i]/Act_wel_Total;
        total += act[i];
    }

    // step 2
    double ran_Istp = G4UniformRand();
    double accumulated_act = 0;
    for(int i=0;i<7;i++)
    {
        double lowerBound = accumulated_act;
        double upperBound = accumulated_act + act[i];
        //cout<<i<<" "<<lowerBound<<" "<<upperBound<<endl;

        if(ran_Istp>=lowerBound && ran_Istp<=upperBound)
        {
            ID = i;
            break;
        }

        accumulated_act += act[i];
    }

    //cout<<"ID "<<ID<<endl;

    return ID;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int B1PrimaryGeneratorAction::GetIsotope_G_or_A()
{
    int ID = 4;

    // step 1
    double act[3];
    act[0] = Act_grd[0]/Act_grd_Total;
    act[1] = Act_grd[1]/Act_grd_Total;
    act[2] = Act_grd[2]/Act_grd_Total;
    //cout<<"act[0] "<<act[0]<<endl;
    //cout<<"act[1] "<<act[1]<<endl;
    //cout<<"act[2] "<<act[2]<<endl;

    // step 2
    double ran_Istp = G4UniformRand();
    double accumulated_act = 0;
    for(int i=0;i<3;i++)
    {
        double lowerBound = accumulated_act;
        double upperBound = accumulated_act + act[i];
        //cout<<i<<" "<<lowerBound<<" "<<upperBound<<endl;

        if(ran_Istp>=lowerBound && ran_Istp<=upperBound)
        {
            ID = i;
            break;
        }

        accumulated_act += act[i];
    }

    return ID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int B1PrimaryGeneratorAction::GetEnvironment()
{
    int ID = 0;

    // step 1
    double total = Act_grd_Total + Act_atm_Total;
    double act[2];
    act[0] = Act_grd_Total/total; // ground
    act[1] = Act_atm_Total/total; // atmosphere
    //cout<<"act[0]"<<act[0]<<endl;
    //cout<<"act[1]"<<act[1]<<endl;

    // step 2
    double ran = G4UniformRand();
    if(ran<act[0]) ID=0;
    else ID=1;

    return ID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double B1PrimaryGeneratorAction::GetEnergy_W(int IstpID)
{
    // step 1 : get isotope
    vector<double> energies;
    vector<double> weights;

    energies = Istp_energies_W[IstpID];
    weights  = Istp_weights__W[IstpID];

    // step 2 : get energy
    double r2 = G4UniformRand();
    //double energy = energies[0];
    double energy = 0;
    double accumulatedW = 0;

    for(int i=0;i<weights.size();i++)
    {
        double lowerBound = accumulatedW;
        double upperBound = accumulatedW + weights[i];
        //cout<<"r2 "<<r2<<", accumulatedW "<<accumulatedW<<", lowerBound "<<lowerBound<<endl;

        if(r2>=lowerBound && r2<=upperBound)
        {
            energy = energies[i];
            break;
        }

        accumulatedW += weights[i];
    }

    //cout<<"energy "<<energy<<" keV"<<endl;

    return energy*keV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double B1PrimaryGeneratorAction::GetEnergy_A_G(int IstpID)
{
    // step 1 : get isotope
    vector<double> energies;
    vector<double> weights;

    if(IstpID==0)
    {
        energies = Istp_129I_energies;
        weights  = Istp_129I_weights;
    }
    else if(IstpID==1)
    {
        energies = Istp_134Cs_energies;
        weights  = Istp_134Cs_weights;
    }
    else
    {
        energies = Istp_137Cs_energies;
        weights  = Istp_137Cs_weights;
    }


    // step 2 : get energy
    double r2 = G4UniformRand();
    //double energy = energies[0];
    double energy = 0;
    double accumulatedW = 0;

    for(int i=0;i<weights.size();i++)
    {
        double lowerBound = accumulatedW;
        double upperBound = accumulatedW + weights[i];
        //cout<<"r2 "<<r2<<", accumulatedW "<<accumulatedW<<", lowerBound "<<lowerBound<<endl;

        if(r2>=lowerBound && r2<=upperBound)
        {
            energy = energies[i];
            break;
        }

        accumulatedW += weights[i];
    }

    //cout<<"energy "<<energy<<" keV"<<endl;

    return energy*keV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void B1PrimaryGeneratorAction::ReadIsotope(string filename, vector<double> &energies, vector<double> &weights)
{
  energies.clear();
  weights. clear();

  // reading
  ifstream read(filename);
  if(read.fail())
  {
      cout<<"Data file dose not exist"<<endl;
      return ;
  }

  double e, w, w2;
  while(!read.eof())
  {
      read>>e>>w>>w2;

      if(read.eof()) break;

      //cout<<e<<" "<<w<<" "<<w2<<endl;

      energies.push_back(e);
      weights. push_back(w);
  }

  read.close();

  // step 2 : 
  double totalW = 0;
  for(int i=0;i<weights.size();i++)
  {
        totalW += weights[i];
  }

  for(int i=0;i<weights.size();i++)
  {
        weights[i] = weights[i]/totalW;
  }

  // step 3 :
  for(int i=0;i<weights.size();i++)
  {
      cout<<filename<<", "<<i<<" "<<energies[i]<<" "<<weights[i]<<endl;
  }


  return;
}
