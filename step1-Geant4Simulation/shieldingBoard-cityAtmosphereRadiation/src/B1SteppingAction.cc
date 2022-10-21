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
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1RunAction.hh"
#include "B1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1RunAction * runAction, B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fRunAction(runAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
	G4Track * track = step->GetTrack();

	G4StepPoint * point1 = step->GetPreStepPoint();
	G4StepPoint * point2 = step->GetPostStepPoint();

	//Use a handle to get the volume
	G4TouchableHandle  touch1 = point1->GetTouchableHandle();
	G4TouchableHandle  touch2 = point2->GetTouchableHandle();

	// get position
	G4ThreeVector pos1 = point1->GetPosition();
	G4ThreeVector pos2 = point2->GetPosition();

	double pos1x = pos1.x()/mm;
	double pos1y = pos1.y()/mm;
	double pos1z = pos1.z()/mm;

	// get energy deposit
	double eDep = step->GetTotalEnergyDeposit();
	double kineticEnergy = track->GetKineticEnergy();

	// get EventID
	double eventID = fEventAction->GetEventID();

	// get TrackID
	double trackID = track->GetTrackID();

	// get parentID
	double parentID = track -> GetParentID();

	// get step ID
	double stepID = track -> GetCurrentStepNumber();

  	// get particle name
  	const G4DynamicParticle *particle = track->GetDynamicParticle();
  	G4String particleName = particle->GetDefinition()->GetParticleName();

	// get Volume Name
	G4VPhysicalVolume * volume1 = touch1->GetVolume();
	G4String VolumeName1 = volume1->GetName();

	if(VolumeName1!="World")
	{
		// get Volume Name
		G4VPhysicalVolume * volume2 = touch2->GetVolume();
		G4String VolumeName2 = volume2->GetName();

		double pos2x = pos2.x()/mm;
		double pos2y = pos2.y()/mm;
		double pos2z = pos2.z()/mm;


		if(VolumeName1.at(0)=='U' && VolumeName2.at(0)=='U')
		{
            //G4cout<<VolumeName1<<G4endl;
            fRunAction->AddEdepUGV(VolumeName1,eDep);
			//std::ofstream write("output.txt",std::ios::app);

            //std::cout<<VolumeName1<<" "<<eventID<<"	"<<stepID<<"	"<<particleName<<std::endl;

			/*
			if(kineticEnergy!=185.7)
			std::cout<<eventID<<"	"<<stepID<<"	"<<particleName<<"	"<<pos2x<<"	"<<pos2y<<"	"<<pos2z<<" "<<kineticEnergy<<std::endl;
			*/

			// much more data
			//write<<eventID<<"	"<<trackID<<"	"<<parentID<<"	"<<stepID<<" "<<particleName<<"	"<<VolumeName1<<"	"<<eDep<<"	"<<pos1x<<"	"<<pos1y<<"	"<<pos1z<<std::endl;

			// simple 
			//write<<eventID<<" "<<VolumeName1<<" "<<particleName<<"	"<<eDep<<std::endl;
			//write.close();
		}
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

