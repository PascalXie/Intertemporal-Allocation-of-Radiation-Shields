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
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <strstream>
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)
{ 
    for(int i=0;i<54;i++)
    {
        UGVEdeps.push_back(0.);
    }

    for(int i=0;i<6;i++)
    for(int j=0;j<9;j++)
    {
        string UGVID;
        strstream ss;
        ss << i;
        ss >> UGVID;

        string ShieldNum;
        strstream ss2;
        ss2 << j;
        ss2 >> ShieldNum;

        string UGVName_current = "UGVID_" + UGVID + "_ShieldNum_" + ShieldNum;
        UGVNames.push_back(UGVName_current);
        G4cout<<"UGVName "<<UGVName_current<<G4endl;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
//  // inform the runManager to save random number seed
//  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
//
//  // reset accumulables to their initial values
//  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
//  accumulableManager->Reset();

	//std::ofstream write("output.txt");
	//write.close();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
	std::ofstream write("output.txt");
    for(int i=0;i<UGVEdeps.size();i++)
    {
        write<<UGVNames[i]<<"  "<<UGVEdeps[i]<<endl;
    }

	write.close();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdepUGV(string UGVname, G4double edep)
{
    for(int i=0;i<UGVNames.size();i++)
    {
        if(UGVname==UGVNames[i]) 
        {
            UGVEdeps[i] = UGVEdeps[i] + edep/keV;
            //G4cout<<"i "<<i<<", name "<<UGVname<<", edep "<<edep/keV<<" keV"<<endl;
            break;
        }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

