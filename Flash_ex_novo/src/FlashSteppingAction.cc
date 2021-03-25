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
#include "FlashRunAction.hh"
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
std::string filename_1 = oss.str();

  KinEnFile.open(filename_1,std::ios_base::app);
  
std::ostringstream oss_of;
oss_of << "Optic_fiber" << ThreadNumber << ".csv" ;
std::string filename_2 = oss_of.str();

  OpticFiber.open(filename_2,std::ios_base::app);
  
std::ostringstream oss_pd;
oss_pd << "Photodetector_" << ThreadNumber << ".csv" ;
std::string filename_3 = oss_pd.str();

  Photodiode.open(filename_3,std::ios_base::app);
  
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashSteppingAction::~FlashSteppingAction()
{
KinEnFile.close();
OpticFiber.close();
Photodiode.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FlashSteppingAction::UserSteppingAction(const G4Step* aStep)
{
//analizziamo i primari che incidono sul cristallo
 G4int eventid = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
 G4StepPoint* postStep = aStep->GetPostStepPoint();
  G4int trackID = aStep->GetTrack()->GetTrackID();
  if ( postStep->GetStepStatus() == fGeomBoundary ){

   G4String volumeName = postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
 if(volumeName == "CrystalLV"&&aStep->GetTrack()->GetTrackID() == 1){

 G4double kineticEnergy = aStep->GetTrack()->GetKineticEnergy();

 
 //G4cout      << "Event ID--->"<<  " " <<  eventid<< " "<< "track ID--->"<<  " " <<  trackID<< " "<< "kineticEnergy--->"<<  " " <<kineticEnergy<< " "<< "Logical Volume --->"<< "  " << volumeName<< " "<< G4endl;
 
 
 if(KinEnFile.is_open()){
 
 KinEnFile<< eventid<< "\t" << kineticEnergy << "\t" <<trackID<<G4endl;
 }
 	}
 if(volumeName == "OF_core_LV" && aStep->GetTrack()->GetTrackID()!= 1){
 if(aStep->GetTrack()->GetCreatorProcess()->GetProcessName() == "Scintillation"){
  //G4cout<< eventid<< "\t" << "Scintillation in core" << "\t" <<trackID<<G4endl;
if(OpticFiber.is_open()){
 
 OpticFiber<< eventid<< "\t" << "Scintillation in core" << "\t" <<trackID<<G4endl;
 }
 }
else if(aStep->GetTrack()->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
 //G4cout<< eventid<< "\t" << "Cerenkov in core" << "\t" <<trackID<<G4endl;
if(OpticFiber.is_open()){
 
 OpticFiber<< eventid<< "\t" << "Cerenkov in core" << "\t" <<trackID<<G4endl;
 }
 }
 }
 if(volumeName == "photocath_surf"  && aStep->GetTrack()->GetTrackID()!= 1){
 if(aStep->GetTrack()->GetCreatorProcess()->GetProcessName() == "Scintillation"){
  G4cout<< eventid<< "\t" << "Scintillation in Photodiode" << "\t" <<trackID<<G4endl;
if(Photodiode.is_open()){
 
 Photodiode<< eventid<< "\t" << "Scintillation in Photodiode" << "\t" <<trackID<<G4endl;
 }
 }
else if(aStep->GetTrack()->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
 G4cout<< eventid<< "\t" << "cerenkov in Photodiode" << "\t" <<trackID<<G4endl;
if(Photodiode.is_open()){
 
 Photodiode<< eventid<< "\t" << "cerenkov in Photodiode" << "\t" <<trackID<<G4endl;
 }
 }
 else{ 
  G4cout<< eventid<< "\t" << "Secondary in Photodiode" << "\t" <<trackID<<G4endl;
 if(Photodiode.is_open()){
 
 Photodiode<< eventid<< "\t" << "Secondary in Photodiode" << "\t" <<trackID<<G4endl;
 } }
 }
  }
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

