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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <strstream>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{

  // UGV locations
  for(int i=0;i<6;i++)
  {
    double x = 12.*m + double(i)*2.*m;
    UGV_xs.push_back(x);
    UGV_ys.push_back(0.);
    G4cout<<"UGV ID"<<i<<", x "<<UGV_xs[i]/m<<" m"<<G4endl;
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  env_mat = nist->FindOrBuildMaterial("G4_AIR");
  concrete_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
  device_mat = nist->FindOrBuildMaterial("G4_Si");
  shield_mat = nist->FindOrBuildMaterial("G4_Pb");
  water_mat = nist->FindOrBuildMaterial("G4_WATER");
  steel_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 52.*m;
  G4double world_sizeZ  = 5.*m;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //     
  // World Fake
  //
  G4double worldFake_sizeXY = world_sizeXY - 2.*mm;
  G4double worldFake_sizeZ  = world_sizeZ  - 2.*mm;
                     
  G4Box* solidWorldFake =    
    new G4Box("WorldFake",                       //its name
       0.5*worldFake_sizeXY, 0.5*worldFake_sizeXY, 0.5*worldFake_sizeZ);     //its size
                
  //G4LogicalVolume* logicWorldFake =                         
  logicWorldFake =                         
    new G4LogicalVolume(solidWorldFake,          //its solid
                        world_mat,           //its material
                        "WorldFake");            //its name

  G4VPhysicalVolume* physWorldFake = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorldFake,            //its logical volume
                      "WorldFake",               //its name
                      logicWorld,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //     
  // Ground
  //
  G4double Ground_RMin = 11.*m;
  G4double Ground_RMax = 25.*m;
  G4double Ground_Dz   = 2.*m;
  G4double Ground_SPhi = 0.;
  G4double Ground_DPhi = 360.*deg;

  double Ground_Locx = 0.*m;
  double Ground_Locy = 0.;
  double Ground_Locz = -1.*Ground_Dz/2.;


  G4Tubs* solidGround = 
   new G4Tubs("Ground",
                Ground_RMin,
                Ground_RMax,
                Ground_Dz/2., 
                Ground_SPhi,
                Ground_DPhi
                );

  G4LogicalVolume* logicGround =                         
    new G4LogicalVolume(solidGround,          //its solid
                        concrete_mat,           //its material
                        "Ground");            //its name

  G4VPhysicalVolume* physGround = 
    new G4PVPlacement(0,                     //no rotation
		              G4ThreeVector(Ground_Locx,Ground_Locy,Ground_Locz),       //at (0,0,0)
                      logicGround,            //its logical volume
                      "Ground",               //its name
                      logicWorldFake,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  //     
  // water well
  //
  G4double Well_RMin = 0.*m;
  G4double Well_RMax = 11.*m;
  G4double Well_Dz   = 2.*m;
  G4double Well_SPhi = 0.;
  G4double Well_DPhi = 360.*deg;

  double Well_Locx = 0.*m;
  double Well_Locy = 0.;
  double Well_Locz = -1.*Well_Dz/2.;


  G4Tubs* solidWell = 
   new G4Tubs("Well",
                Well_RMin,
                Well_RMax,
                Well_Dz/2., 
                Well_SPhi,
                Well_DPhi
                );

  G4LogicalVolume* logicWell =                         
    new G4LogicalVolume(solidWell,          //its solid
                        water_mat,           //its material
                        "Well");            //its name

  G4VPhysicalVolume* physWell = 
    new G4PVPlacement(0,                     //no rotation
		              G4ThreeVector(Well_Locx,Well_Locy,Well_Locz),       //at (0,0,0)
                      logicWell,            //its logical volume
                      "Well",               //its name
                      logicWorldFake,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //     
  // steel pressure vessel
  //
  G4double Vessel_RMin = 0.*m;
  G4double Vessel_RMax = 11.*m;
  G4double Vessel_Dz   = 2.*m;
  G4double Vessel_SPhi = 0.;
  G4double Vessel_DPhi = 360.*deg;

  double Vessel_Locx = 0.*m;
  double Vessel_Locy = 0.;
  double Vessel_Locz = Vessel_Dz/2.;


  G4Tubs* solidVessel = 
   new G4Tubs("Vessel",
                Vessel_RMin,
                Vessel_RMax,
                Vessel_Dz/2., 
                Vessel_SPhi,
                Vessel_DPhi
                );

  G4LogicalVolume* logicVessel =                         
    new G4LogicalVolume(solidVessel,          //its solid
                        steel_mat,           //its material
                        "Vessel");            //its name

  G4VPhysicalVolume* physVessel = 
    new G4PVPlacement(0,                     //no rotation
		              G4ThreeVector(Vessel_Locx,Vessel_Locy,Vessel_Locz),       //at (0,0,0)
                      logicVessel,            //its logical volume
                      "Vessel",               //its name
                      logicWorldFake,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //     
  // Wall 
  //
  G4double Wall_RMin = 24.*m;
  G4double Wall_RMax = 25.*m;
  G4double Wall_Dz   = 2.*m;
  G4double Wall_SPhi = 0.;
  G4double Wall_DPhi = 360.*deg;

  double Wall_Locx = 0.*m;
  double Wall_Locy = 0.;
  double Wall_Locz = Wall_Dz/2.;


  G4Tubs* solidWall = 
   new G4Tubs("Wall",
                Wall_RMin,
                Wall_RMax,
                Wall_Dz/2., 
                Wall_SPhi,
                Wall_DPhi
                );

  G4LogicalVolume* logicWall =                         
    new G4LogicalVolume(solidWall,          //its solid
                        concrete_mat,           //its material
                        "Wall");            //its name

  G4VPhysicalVolume* physWall = 
    new G4PVPlacement(0,                     //no rotation
		              G4ThreeVector(Wall_Locx,Wall_Locy,Wall_Locz),       //at (0,0,0)
                      logicWall,            //its logical volume
                      "Wall",               //its name
                      logicWorldFake,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  // UGV and shields
  int UGVID = 0;
  for(int i=0;i<6;i++)
  {
    int UGVID = i;
    int isUGVGood = LocateUGV(UGVID);
  }



  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int B1DetectorConstruction::LocateUGV(int UGVID)
{
  string UGVID_s;
  strstream ss;
  ss << UGVID;
  ss >> UGVID_s;
  UGVID_s = "UGVID_" + UGVID_s + "_ShieldNum";
  G4cout<<"B1DetectorConstruction UGVname "<<UGVID_s<<G4endl;


  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

	//     
	// UGV devices
	//
	double UGVSizex = 50.*cm;
	double UGVSizey = 16.*cm;
	double UGVSizez = 4.*cm;
	double UGVLocx = UGV_xs[UGVID];
	double UGVLocy = 0.;
	double UGVLocz = 0.2*m;

	G4Box* solidUGV =
	  new G4Box("UGV",              //its name
	            UGVSizex/2, UGVSizey/2, UGVSizez/2);     //its size

	G4LogicalVolume* logicUGV =                         
	  new G4LogicalVolume(solidUGV,          //its solid
	                      device_mat,           //its material
	                      "UGV");            //its name

	//G4VPhysicalVolume* physDet = 
	//	  new G4PVPlacement(0,                     //no rotation
	//	                    G4ThreeVector(UGVLocx,UGVLocy,UGVLocz),       //at (0,0,0)
	//	                    logicUGV,            //its logical volume
	//	                    "UGV",               //its name
	//	                    logicWorldFake,                     //its mother logical volume
	//	                    false,                 //no boolean operation
	//	                    0,                     //copy number
	//	                    checkOverlaps);        //overlaps checking

	//     
	// UGV shields
	//
    vector<double> SSizexs, SSizeys, SSizezs;
    vector<double> SLocxs, SLocys, SLoczs;
    int NumberShields = 6;

	//     
	// UGV shield : S1
	//
	double S1Sizex = 500.*mm;
	double S1Sizey = 50. *mm;
	double S1Sizez = 3.  *mm;
	double S1Locx  = UGVLocx;
	double S1Locy  = -1.*(S1Sizey + 1.6*cm);
	double S1Locz  = UGVLocz - 5.2*cm/2.;

    SSizexs.push_back(S1Sizex);
    SSizeys.push_back(S1Sizey);
    SSizezs.push_back(S1Sizez);
    SLocxs .push_back(S1Locx);
    SLocys .push_back(S1Locy);
    SLoczs .push_back(S1Locz);

	G4Box* solidS1 =
	  new G4Box("S1",              //its name
	            S1Sizex/2, S1Sizey/2, S1Sizez/2);     //its size

	G4LogicalVolume* logicS1 =                         
	  new G4LogicalVolume(solidS1,          //its solid
	                      shield_mat,           //its material
	                      "S1");            //its name

    //if(NumberShields>=1)
    //{
	//    G4VPhysicalVolume* physS1 = 
	//    	  new G4PVPlacement(0,                     //no rotation
	//    	                    G4ThreeVector(S1Locx,S1Locy,S1Locz),       //at (0,0,0)
	//    	                    logicS1,            //its logical volume
	//    	                    "S1",               //its name
	//    	                    logicWorldFake,                     //its mother logical volume
	//    	                    false,                 //no boolean operation
	//    	                    0,                     //copy number
	//    	                    checkOverlaps);        //overlaps checking
    //}

	//     
	// UGV shield : S2
	//
	double S2Sizex = 500.*mm;
	double S2Sizey = 50. *mm;
	double S2Sizez = 3.  *mm;
	double S2Locx  = UGVLocx;
	double S2Locy  = 0.;
	double S2Locz  = UGVLocz - 5.2*cm/2.;

    SSizexs.push_back(S2Sizex);
    SSizeys.push_back(S2Sizey);
    SSizezs.push_back(S2Sizez);
    SLocxs .push_back(S2Locx);
    SLocys .push_back(S2Locy);
    SLoczs .push_back(S2Locz);

	G4Box* solidS2 =
	  new G4Box("S2",              //its name
	            S2Sizex/2, S2Sizey/2, S2Sizez/2);     //its size

	G4LogicalVolume* logicS2 =                         
	  new G4LogicalVolume(solidS2,          //its solid
	                      shield_mat,           //its material
	                      "S2");            //its name

    //if(NumberShields>=2)
    //{
	//    G4VPhysicalVolume* physS2 = 
	//    	  new G4PVPlacement(0,                     //no rotation
	//    	                    G4ThreeVector(S2Locx,S2Locy,S2Locz),       //at (0,0,0)
	//    	                    logicS2,            //its logical volume
	//    	                    "S2",               //its name
	//    	                    logicWorldFake,                     //its mother logical volume
	//    	                    false,                 //no boolean operation
	//    	                    0,                     //copy number
	//    	                    checkOverlaps);        //overlaps checking
    //}

	//     
	// UGV shield : S3
	//
	double S3Sizex = 500.*mm;
	double S3Sizey = 50. *mm;
	double S3Sizez = 3.  *mm;
	double S3Locx  = UGVLocx;
	double S3Locy  = S3Sizey + 1.6*cm;
	double S3Locz  = UGVLocz - 5.2*cm/2.;

    SSizexs.push_back(S3Sizex);
    SSizeys.push_back(S3Sizey);
    SSizezs.push_back(S3Sizez);
    SLocxs .push_back(S3Locx);
    SLocys .push_back(S3Locy);
    SLoczs .push_back(S3Locz);

	G4Box* solidS3 =
	  new G4Box("S3",              //its name
	            S3Sizex/2, S3Sizey/2, S3Sizez/2);     //its size

	G4LogicalVolume* logicS3 =                         
	  new G4LogicalVolume(solidS3,          //its solid
	                      shield_mat,           //its material
	                      "S3");            //its name

    //if(NumberShields>=3)
    //{
	//    G4VPhysicalVolume* physS3 = 
	//    	  new G4PVPlacement(0,                     //no rotation
	//    	                    G4ThreeVector(S3Locx,S3Locy,S3Locz),       //at (0,0,0)
	//    	                    logicS3,            //its logical volume
	//    	                    "S3",               //its name
	//    	                    logicWorldFake,                     //its mother logical volume
	//    	                    false,                 //no boolean operation
	//    	                    0,                     //copy number
	//    	                    checkOverlaps);        //overlaps checking
    //}

	//     
	// UGV shield : S4
	//
	double S4Sizex = 500.*mm;
	double S4Sizey = 3.  *mm;
	double S4Sizez = 50. *mm;
	double S4Locx  = UGVLocx;
	double S4Locy  = -1.*(1.5*50.*mm + 1.6*cm + 0.2*cm);
	double S4Locz  = UGVLocz;

    SSizexs.push_back(S4Sizex);
    SSizeys.push_back(S4Sizey);
    SSizezs.push_back(S4Sizez);
    SLocxs .push_back(S4Locx);
    SLocys .push_back(S4Locy);
    SLoczs .push_back(S4Locz);

	G4Box* solidS4 =
	  new G4Box("S4",              //its name
	            S4Sizex/2, S4Sizey/2, S4Sizez/2);     //its size

	G4LogicalVolume* logicS4 =                         
	  new G4LogicalVolume(solidS4,          //its solid
	                      shield_mat,           //its material
	                      "S4");            //its name

    //if(NumberShields>=4)
    //{
	//    G4VPhysicalVolume* physS4 = 
	//    	  new G4PVPlacement(0,                     //no rotation
	//    	                    G4ThreeVector(S4Locx,S4Locy,S4Locz),       //at (0,0,0)
	//    	                    logicS4,            //its logical volume
	//    	                    "S4",               //its name
	//    	                    logicWorldFake,                     //its mother logical volume
	//    	                    false,                 //no boolean operation
	//    	                    0,                     //copy number
	//    	                    checkOverlaps);        //overlaps checking
    //}


	//     
	// UGV shield : S5
	//
	double S5Sizex = 500.*mm;
	double S5Sizey = 3.  *mm;
	double S5Sizez = 50. *mm;
	double S5Locx  = UGVLocx;
	double S5Locy  = 1.5*50.*mm + 1.6*cm + 0.2*cm;
	double S5Locz  = UGVLocz;

    SSizexs.push_back(S5Sizex);
    SSizeys.push_back(S5Sizey);
    SSizezs.push_back(S5Sizez);
    SLocxs .push_back(S5Locx);
    SLocys .push_back(S5Locy);
    SLoczs .push_back(S5Locz);

	G4Box* solidS5 =
	  new G4Box("S5",              //its name
	            S5Sizex/2, S5Sizey/2, S5Sizez/2);     //its size

	G4LogicalVolume* logicS5 =                         
	  new G4LogicalVolume(solidS5,          //its solid
	                      shield_mat,           //its material
	                      "S5");            //its name

    //if(NumberShields>=5)
    //{
	//    G4VPhysicalVolume* physS5 = 
	//    	  new G4PVPlacement(0,                     //no rotation
	//    	                    G4ThreeVector(S5Locx,S5Locy,S5Locz),       //at (0,0,0)
	//    	                    logicS5,            //its logical volume
	//    	                    "S5",               //its name
	//    	                    logicWorldFake,                     //its mother logical volume
	//    	                    false,                 //no boolean operation
	//    	                    0,                     //copy number
	//    	                    checkOverlaps);        //overlaps checking
    //}

	//     
	// UGV shield : S6
	//
	double S6Sizex = 500.*mm;
	double S6Sizey = 50. *mm;
	double S6Sizez = 3.  *mm;
	double S6Locx  = UGVLocx;
	double S6Locy  = -1.*(50.*mm + 1.6*cm);
	double S6Locz  = UGVLocz + 5.2*cm/2.;

    SSizexs.push_back(S6Sizex);
    SSizeys.push_back(S6Sizey);
    SSizezs.push_back(S6Sizez);
    SLocxs .push_back(S6Locx);
    SLocys .push_back(S6Locy);
    SLoczs .push_back(S6Locz);

	G4Box* solidS6 =
	  new G4Box("S6",              //its name
	            S6Sizex/2, S6Sizey/2, S6Sizez/2);     //its size

	G4LogicalVolume* logicS6 =                         
	  new G4LogicalVolume(solidS6,          //its solid
	                      shield_mat,           //its material
	                      "S6");            //its name

    //if(NumberShields>=6)
    //{
	//    G4VPhysicalVolume* physS6 = 
	//    	  new G4PVPlacement(0,                     //no rotation
	//    	                    G4ThreeVector(S6Locx,S6Locy,S6Locz),       //at (0,0,0)
	//    	                    logicS6,            //its logical volume
	//    	                    "S6",               //its name
	//    	                    logicWorldFake,                     //its mother logical volume
	//    	                    false,                 //no boolean operation
	//    	                    0,                     //copy number
	//    	                    checkOverlaps);        //overlaps checking
    //}

	//     
	// UGV shield : S7
	//
	double S7Sizex = 500.*mm;
	double S7Sizey = 50. *mm;
	double S7Sizez = 3.  *mm;
	double S7Locx  = UGVLocx;
	double S7Locy  = 0.;
	double S7Locz  = UGVLocz + 5.2*cm/2.;

    SSizexs.push_back(S7Sizex);
    SSizeys.push_back(S7Sizey);
    SSizezs.push_back(S7Sizez);
    SLocxs .push_back(S7Locx);
    SLocys .push_back(S7Locy);
    SLoczs .push_back(S7Locz);

	G4Box* solidS7 =
	  new G4Box("S7",              //its name
	            S7Sizex/2, S7Sizey/2, S7Sizez/2);     //its size

	G4LogicalVolume* logicS7 =                         
	  new G4LogicalVolume(solidS7,          //its solid
	                      shield_mat,           //its material
	                      "S7");            //its name

    //if(NumberShields>=7)
    //{
	//    G4VPhysicalVolume* physS7 = 
	//    	  new G4PVPlacement(0,                     //no rotation
	//    	                    G4ThreeVector(S7Locx,S7Locy,S7Locz),       //at (0,0,0)
	//    	                    logicS7,            //its logical volume
	//    	                    "S7",               //its name
	//    	                    logicWorldFake,                     //its mother logical volume
	//    	                    false,                 //no boolean operation
	//    	                    0,                     //copy number
	//    	                    checkOverlaps);        //overlaps checking
    //}



	//     
	// UGV shield : S8
	//
	double S8Sizex = 500.*mm;
	double S8Sizey = 50. *mm;
	double S8Sizez = 3.  *mm;
	double S8Locx  = UGVLocx;
	double S8Locy  = 50.*mm + 1.6*cm;
	double S8Locz  = UGVLocz + 5.2*cm/2.;

    SSizexs.push_back(S8Sizex);
    SSizeys.push_back(S8Sizey);
    SSizezs.push_back(S8Sizez);
    SLocxs .push_back(S8Locx);
    SLocys .push_back(S8Locy);
    SLoczs .push_back(S8Locz);

	G4Box* solidS8 =
	  new G4Box("S8",              //its name
	            S8Sizex/2, S8Sizey/2, S8Sizez/2);     //its size

	G4LogicalVolume* logicS8 =                         
	  new G4LogicalVolume(solidS8,          //its solid
	                      shield_mat,           //its material
	                      "S8");            //its name

    //if(NumberShields>=8)
    //{
	//    G4VPhysicalVolume* physS8 = 
	//    	  new G4PVPlacement(0,                     //no rotation
	//    	                    G4ThreeVector(S8Locx,S8Locy,S8Locz),       //at (0,0,0)
	//    	                    logicS8,            //its logical volume
	//    	                    "S8",               //its name
	//    	                    logicWorldFake,                     //its mother logical volume
	//    	                    false,                 //no boolean operation
	//    	                    0,                     //copy number
	//    	                    checkOverlaps);        //overlaps checking
    //}

    //
    double shieldNum = 0.;
    double theta = 0*deg + 10.*deg*shieldNum;
    double x = UGVLocx*cos(theta);
    double y = UGVLocx*sin(theta);


    double yshift = 3.*50.*cm;

	//     
	// UGV, 0 shield
	//
    yshift = -5.*50.*cm;

    shieldNum = 0.;
    theta = 0*deg + 10.*deg*shieldNum;
    x = UGVLocx*cos(theta);
    y = UGVLocx*sin(theta);

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,y,UGVLocz),       //at (0,0,0)
	                  logicUGV,            //its logical volume
	                  UGVID_s+"_0",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking


	//     
	// UGV, 1 shield
	//
    yshift = -4.*50.*cm;

    shieldNum = 1.;
    theta = 0*deg + 10.*deg*shieldNum;
    x = UGVLocx*cos(theta);
    y = UGVLocx*sin(theta);

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,y,UGVLocz),       //at (0,0,0)
	                  logicUGV,            //its logical volume
	                  UGVID_s+"_1",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S1Locy+y,S1Locz),       //at (0,0,0)
	                  logicS1,            //its logical volume
	                  "S1",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking



	//     
	// UGV, 2 shields
	//
    yshift = -3.*50.*cm;

    shieldNum = 2.;
    theta = 0*deg + 10.*deg*shieldNum;
    x = UGVLocx*cos(theta);
    y = UGVLocx*sin(theta);

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,y,UGVLocz),       //at (0,0,0)
	                  logicUGV,            //its logical volume
	                  UGVID_s+"_2",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S1Locy+y,S1Locz),       //at (0,0,0)
	                  logicS1,            //its logical volume
	                  "S1",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S2Locy+y,S2Locz),       //at (0,0,0)
	                  logicS2,            //its logical volume
	                  "S2",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	//     
	// UGV, 3 shields
	//
    yshift = -2.*50.*cm;

    shieldNum = 3.;
    theta = 0*deg + 10.*deg*shieldNum;
    x = UGVLocx*cos(theta);
    y = UGVLocx*sin(theta);

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x, y,UGVLocz),       //at (0,0,0)
	                  logicUGV,            //its logical volume
	                  UGVID_s+"_3",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S1Locy+y,S1Locz),       //at (0,0,0)
	                  logicS1,            //its logical volume
	                  "S1",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S2Locy+y,S2Locz),       //at (0,0,0)
	                  logicS2,            //its logical volume
	                  "S2",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S3Locy+y,S3Locz),       //at (0,0,0)
	                  logicS3,            //its logical volume
	                  "S3",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	//     
	// UGV, 4 shields
	//
    yshift = -1.*50.*cm;

    shieldNum = 4.;
    theta = 0*deg + 10.*deg*shieldNum;
    x = UGVLocx*cos(theta);
    y = UGVLocx*sin(theta);

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x, y,UGVLocz),       //at (0,0,0)
	                  logicUGV,            //its logical volume
	                  UGVID_s+"_4",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S1Locy+y,S1Locz),       //at (0,0,0)
	                  logicS1,            //its logical volume
	                  "S1",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S2Locy+y,S2Locz),       //at (0,0,0)
	                  logicS2,            //its logical volume
	                  "S2",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S3Locy+y,S3Locz),       //at (0,0,0)
	                  logicS3,            //its logical volume
	                  "S3",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S4Locy+y,S4Locz),       //at (0,0,0)
	                  logicS4,            //its logical volume
	                  "S4",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	//     
	// UGV, 5 shields
	//
    yshift = 1.*50.*cm;

    shieldNum = 5.;
    theta = 0*deg + 10.*deg*shieldNum;
    x = UGVLocx*cos(theta);
    y = UGVLocx*sin(theta);

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,y,UGVLocz),       //at (0,0,0)
	                  logicUGV,            //its logical volume
	                  UGVID_s+"_5",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S1Locy+y,S1Locz),       //at (0,0,0)
	                  logicS1,            //its logical volume
	                  "S1",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S2Locy+y,S2Locz),       //at (0,0,0)
	                  logicS2,            //its logical volume
	                  "S2",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S3Locy+y,S3Locz),       //at (0,0,0)
	                  logicS3,            //its logical volume
	                  "S3",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S4Locy+y,S4Locz),       //at (0,0,0)
	                  logicS4,            //its logical volume
	                  "S4",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S5Locy+y,S5Locz),       //at (0,0,0)
	                  logicS5,            //its logical volume
	                  "S5",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking


	//     
	// UGV, 6 shields
	//
    yshift = 2.*50.*cm;

    shieldNum = 6.;
    theta = 0*deg + 10.*deg*shieldNum;
    x = UGVLocx*cos(theta);
    y = UGVLocx*sin(theta);

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,y,UGVLocz),       //at (0,0,0)
	                  logicUGV,            //its logical volume
	                  UGVID_s+"_6",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S1Locy+y,S1Locz),       //at (0,0,0)
	                  logicS1,            //its logical volume
	                  "S1",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S2Locy+y,S2Locz),       //at (0,0,0)
	                  logicS2,            //its logical volume
	                  "S2",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S3Locy+y,S3Locz),       //at (0,0,0)
	                  logicS3,            //its logical volume
	                  "S3",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S4Locy+y,S4Locz),       //at (0,0,0)
	                  logicS4,            //its logical volume
	                  "S4",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S5Locy+y,S5Locz),       //at (0,0,0)
	                  logicS5,            //its logical volume
	                  "S5",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S6Locy+y,S6Locz),       //at (0,0,0)
	                  logicS6,            //its logical volume
	                  "S6",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking


	//     
	// UGV, 7 shields
	//
    yshift = 3.*50.*cm;

    shieldNum = 7.;
    theta = 0*deg + 10.*deg*shieldNum;
    x = UGVLocx*cos(theta);
    y = UGVLocx*sin(theta);

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,y,UGVLocz),       //at (0,0,0)
	                  logicUGV,            //its logical volume
	                  UGVID_s+"_7",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S1Locy+y,S1Locz),       //at (0,0,0)
	                  logicS1,            //its logical volume
	                  "S1",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S2Locy+y,S2Locz),       //at (0,0,0)
	                  logicS2,            //its logical volume
	                  "S2",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S3Locy+y,S3Locz),       //at (0,0,0)
	                  logicS3,            //its logical volume
	                  "S3",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S4Locy+y,S4Locz),       //at (0,0,0)
	                  logicS4,            //its logical volume
	                  "S4",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S5Locy+y,S5Locz),       //at (0,0,0)
	                  logicS5,            //its logical volume
	                  "S5",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S6Locy+y,S6Locz),       //at (0,0,0)
	                  logicS6,            //its logical volume
	                  "S6",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S7Locy+y,S7Locz),       //at (0,0,0)
	                  logicS7,            //its logical volume
	                  "S7",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	//     
	// UGV, 8 shield
	//
    yshift = 4.*50.*cm;

    shieldNum = 8.;
    theta = 0*deg + 10.*deg*shieldNum;
    x = UGVLocx*cos(theta);
    y = UGVLocx*sin(theta);

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,y,UGVLocz),       //at (0,0,0)
	                  logicUGV,            //its logical volume
	                  UGVID_s+"_8",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S1Locy+y,S1Locz),       //at (0,0,0)
	                  logicS1,            //its logical volume
	                  "S1",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S2Locy+y,S2Locz),       //at (0,0,0)
	                  logicS2,            //its logical volume
	                  "S2",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S3Locy+y,S3Locz),       //at (0,0,0)
	                  logicS3,            //its logical volume
	                  "S3",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S4Locy+y,S4Locz),       //at (0,0,0)
	                  logicS4,            //its logical volume
	                  "S4",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S5Locy+y,S5Locz),       //at (0,0,0)
	                  logicS5,            //its logical volume
	                  "S5",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S6Locy+y,S6Locz),       //at (0,0,0)
	                  logicS6,            //its logical volume
	                  "S6",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S7Locy+y,S7Locz),       //at (0,0,0)
	                  logicS7,            //its logical volume
	                  "S7",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

	new G4PVPlacement(0,                     //no rotation
	                  G4ThreeVector(x,S8Locy+y,S8Locz),       //at (0,0,0)
	                  logicS8,            //its logical volume
	                  "S8",               //its name
	                  logicWorldFake,                     //its mother logical volume
	                  false,                 //no boolean operation
	                  0,                     //copy number
	                  checkOverlaps);        //overlaps checking

  return 1;
}
