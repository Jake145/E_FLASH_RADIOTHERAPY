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
/// \file FlashStackingAction.cc
/// \brief Implementation of the FlashStackingAction class

#include "FlashStackingAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include <sstream>
#include <string>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashStackingAction::FlashStackingAction()
    : G4UserStackingAction(), fScintillationCounter(0), fCerenkovCounter(0),fBremstralung(0),fFluo(0) {
  std::ostringstream oss_1;
  oss_1 << "Optic_" << G4Threading::G4GetThreadId() << ".csv";
  std::string filename_1 = oss_1.str();
  // G4String filename=("Kinetic_E_crystal_%d.csv",ThreadNumber);
  OpticFile.open(filename_1, std::ios_base::app);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashStackingAction::~FlashStackingAction() { OpticFile.close(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
FlashStackingAction::ClassifyNewTrack(const G4Track *aTrack) {

/*if(aTrack->GetDefinition()==G4Gamma::GammaDefinition()){
if((aTrack->GetVolume()->GetLogicalVolume()->GetName() == "OF_core_LV"||aTrack->GetVolume()->GetLogicalVolume()->GetName() == "OF_clad_LV"||aTrack->GetVolume()->GetLogicalVolume()->GetName() == "OF_cladding_LV") || aTrack->GetVolume()->GetLogicalVolume()->GetName() == "CrystalLV" ){

if(aTrack->GetCreatorProcess()->GetProcessName()== "eBrem")
{
fBremstralung++;
}
else if (aTrack->GetCreatorProcess()->GetProcessName()== "eIoni")
{
fFluo++;
}
}
}*/
  if (aTrack->GetDefinition() ==
      G4OpticalPhoton::OpticalPhotonDefinition()) { 
      if ((aTrack->GetVolume()->GetLogicalVolume()->GetName() == "OF_core_LV"||aTrack->GetVolume()->GetLogicalVolume()->GetName() == "OF_clad_LV"||aTrack->GetVolume()->GetLogicalVolume()->GetName() == "OF_cladding_LV") || aTrack->GetVolume()->GetLogicalVolume()->GetName() == "CrystalLV" ){
        
      

        if (aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation") {
          /*if(aTrack->GetParentID()==0) {fScintillationCounter++;}
          else if(aTrack->GetParentID()>0) {fBremstralung++;}*/
          fScintillationCounter++;
           //G4cout<<"found a scintillation in:"<<""<<aTrack->GetVolume()->GetName()<<G4endl;
        }
        if (aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov") {
          /*if(aTrack->GetParentID()==0) {fCerenkovCounter++;}
                    else if(aTrack->GetParentID()>0) {fFluo++;}*/
                    fCerenkovCounter++;
          //G4cout << "found a cherenkov in : "<<" "<<aTrack->GetVolume()->GetName()<<G4endl;
        }
      }
      else 
      	return fKill; 
    
  }

  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FlashStackingAction::NewStage()

{
  if (fScintillationCounter != 0 || fCerenkovCounter != 0) {
   
    if (OpticFile.is_open()) {

      OpticFile
          << "event ID: "
          << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
          << "\t"
          << "Primary scintillation events:" << fScintillationCounter << "\t"
          << "Primary Cerenkov events: " << fCerenkovCounter << "\t"<< G4endl;
          //<< "Secondary Scintillation events: "<<fBremstralung<<"\t"
          //<<"Secondary Cerenkov Events: "<<fFluo<< G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FlashStackingAction::PrepareNewEvent() {
  fScintillationCounter = 0;
  fCerenkovCounter = 0;
}
