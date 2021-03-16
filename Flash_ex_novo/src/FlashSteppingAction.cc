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
#include "G4Threading.hh"
#include "G4OpticalPhoton.hh"
#include <string>
#include <sstream>
//#include "FlashAnalysis.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashSteppingAction::FlashSteppingAction(FlashEventAction* eventAction)
: G4UserSteppingAction(),
  ThreadNumber(G4Threading::G4GetThreadId()),fEventAction(eventAction){
  //G4int ThreadNumber=G4Threading::G4GetThreadId();
std::ostringstream oss;
oss << "Kinetic_E_crystal_" << ThreadNumber << ".csv" ;
std::string filename = oss.str();
  //G4String filename=("Kinetic_E_crystal_%d.csv",ThreadNumber);
  KinEnFile.open(filename,std::ios_base::app);
  
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashSteppingAction::~FlashSteppingAction()
{KinEnFile.close();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FlashSteppingAction::UserSteppingAction(const G4Step* aStep)
{

 G4StepPoint* preStep = aStep->GetPreStepPoint();
  if ( preStep->GetStepStatus() == fGeomBoundary ){
 G4String volumeName = preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
 if(volumeName == "CrystalLV"&&aStep->GetTrack()->GetTrackID()==1){
 G4int eventid = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
 G4double kineticEnergy = aStep->GetTrack()->GetKineticEnergy();
 G4int trackID = aStep->GetTrack()->GetTrackID();
 
 G4cout      << "Event ID--->"<<  " " <<  eventid<< " "<< "kineticEnergy--->"<<  " " <<kineticEnergy<< " "<< "Logical Volume --->"<< "  " << volumeName<< " "<< G4endl;
 
 /*if(aStep->GetTrack()->GetTrackID()==1){
 		if(aStep->GetPreStepPoint()->GetProcessDefinedStep()!=0){
 		if(aStep->GetTrack()->GetNextVolume()->GetName()=="crystal" && aStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName()=="Transportation"){
 		G4int eventid = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
 G4double kineticEnergy = aStep->GetTrack()->GetKineticEnergy();
 G4int trackID = aStep->GetTrack()->GetTrackID();
  G4String volumeName = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName();
 		G4cout      << "Event ID--->"<<  " " <<  eventid<< " "<< "kineticEnergy--->"<<  " " <<kineticEnergy<< " "<< "Logical Volume --->"<< "  " << volumeName<< " "<< G4endl;*/
 
 if(KinEnFile.is_open()){
 
 KinEnFile<< eventid<< "\t" << kineticEnergy << "\t" <<trackID<<G4endl;
 }
 	}
 
 
  }
G4Track* track = aStep->GetTrack();
G4ParticleDefinition* particleType = track->GetDefinition();
/*
if (track->GetTrackID()!=1 && particleType == G4OpticalPhoton::OpticalPhotonDefinition()){if(preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName()!="CrystalLV"&&preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName()!="OF_core_LV"){ track->SetTrackStatus(fStopAndKill); }
  }*/

 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

