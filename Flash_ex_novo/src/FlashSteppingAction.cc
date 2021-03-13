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

G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
const G4VProcess* CurrentProcess=preStepPoint->GetProcessDefinedStep();

 //G4String volumeName = preStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetName();

 
 if (CurrentProcess != 0) {
const G4String & StepProcessName = CurrentProcess->GetProcessName();
G4String volumePos = aStep->GetTrack()->GetNextVolume()->GetName();

 if(StepProcessName== "Transportation" && volumePos == "crystal"){
 G4int eventid = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
 G4double kineticEnergy = aStep->GetTrack()->GetKineticEnergy();
 G4int trackID = aStep->GetTrack()->GetTrackID();
 
 G4cout      << "Event ID--->"<<  " " <<  eventid<< " "<< "kineticEnergy--->"<<  " " <<kineticEnergy<< " "<< " Volume --->"<< "  " << volumePos<< " "<< G4endl;
 if(KinEnFile.is_open()){
 
 KinEnFile<< eventid<< "\t" << kineticEnergy << "\t" <<trackID<<"\t" <<volumePos<<G4endl;
 };
  };
  };
  }

 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

