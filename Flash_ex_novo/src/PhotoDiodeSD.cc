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
/// \file optical/LXe/src/PhotoDiodeSD.cc
/// \brief Implementation of the PhotoDiodeSD class
//
//
#include "PhotoDiodeSD.hh"
#include "PhotoDiodeHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhotoDiodeSD::PhotoDiodeSD(G4String name)
  : G4VSensitiveDetector(name)
{
  fphotoCollection = nullptr;
  collectionName.insert("PhtoDiodeCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhotoDiodeSD::~PhotoDiodeSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhotoDiodeSD::Initialize(G4HCofThisEvent* hitsCE){
  fphotoCollection = new PhotoDiodeHitsCollection
                      (SensitiveDetectorName,collectionName[0]);
  //A way to keep all the hits of this event in one place if needed
  static G4int hitsCID = -1;
  if(hitsCID<0){
    hitsCID = GetCollectionID(0);
  }
  hitsCE->AddHitsCollection( hitsCID, fphotoCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PhotoDiodeSD::ProcessHits(G4Step* ,G4TouchableHistory* ){
  return false;
}
G4bool PhotoDiodeSD::ProcessHits_constStep(G4Step* aStep,G4TouchableHistory* ){

	G4Track* track = aStep->GetTrack();
	G4ParticleDefinition* particleType = track->GetDefinition();
	if (particleType != G4OpticalPhoton::OpticalPhotonDefinition()){
	return false;
	}
	
	hit_photo = new PhotoDiodeHit(); 

    	hit_photo->IncPhotonCount_photo(); 
    	
	if(track->GetCreatorProcess()->GetProcessName() == "Scintillation"){
		hit_photo->IncCerenkovCount_photo();
		}
	      if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
		hit_photo->IncScintillationCount_photo();
		}
                   
           return true;
	}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


