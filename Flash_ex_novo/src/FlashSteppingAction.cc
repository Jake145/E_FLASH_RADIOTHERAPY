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
/// \file FlashSteppingAction.cc
/// \brief Implementation of the FlashSteppingAction class

#include "FlashSteppingAction.hh"
#include "FlashEventAction.hh"
#include "FlashDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
//#include "FlashAnalysis.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashSteppingAction::FlashSteppingAction(FlashEventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashSteppingAction::~FlashSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FlashSteppingAction::UserSteppingAction(const G4Step*)
{

/*
 G4String volumeName = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName();
 //if(volumeName == "CrystalLV"){
 G4int eventid = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
 G4double kineticEnergy = aStep->GetTrack()->GetKineticEnergy();
 G4int trackID = aStep->GetTrack()->GetTrackID();
 
  G4AnalysisManager* analysisManager =   G4AnalysisManager::Instance();
 analysisManager->FillNtupleIColumn(2,0, trackID);
   analysisManager->FillNtupleDColumn(2,1, kineticEnergy); 
      analysisManager->FillNtupleIColumn(2,2, eventid); 
            analysisManager->FillNtupleSColumn(2,3, volumeName); 
    analysisManager->AddNtupleRow();
 //};*/
  }

 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

