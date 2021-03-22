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
/// \file optical/LXe/src/OpticFiberSD.cc
/// \brief Implementation of the OpticFiberSD class
//
//
#include "OpticFiberSD.hh"
#include "OpticFiberHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4OpticalPhoton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpticFiberSD::OpticFiberSD(G4String name)
  : G4VSensitiveDetector(name)
{
  foptCollection = nullptr;
  collectionName.insert("OpticFiberHitsCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpticFiberSD::~OpticFiberSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpticFiberSD::Initialize(G4HCofThisEvent* hitsCE){
  foptCollection = new OpticFiberHitsCollection
                      (SensitiveDetectorName,collectionName[0]);
  //A way to keep all the hits of this event in one place if needed
  static G4int hitsCID = -1;
  if(hitsCID<0){
    hitsCID = GetCollectionID(0);
  }
  hitsCE->AddHitsCollection( hitsCID, foptCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool OpticFiberSD::ProcessHits(G4Step* step,G4TouchableHistory* ){

  G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4StepPoint* postStepPoint = step->GetPostStepPoint();
  G4Track* track = step->GetTrack();
	G4ParticleDefinition* particleType = track->GetDefinition();
	if (particleType != G4OpticalPhoton::OpticalPhotonDefinition()){
	return false;
	}
	if ( preStepPoint->GetStepStatus() == fGeomBoundary){ //&& preStepPoint->GetTouchableHandle()->GetVolume()->GetName()=="crystal";){
	OpticFiberHit* hit = new OpticFiberHit(); 

    	hit->IncPhotonCount(); 
    	
	if(track->GetCreatorProcess()->GetProcessName() == "Scintillation"){
		hit->IncCerenkovCount();
		}
	      if(track->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
		hit->IncScintillationCount();
		
		    foptCollection->insert( hit );
		}
		}
                   
           return true;;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


