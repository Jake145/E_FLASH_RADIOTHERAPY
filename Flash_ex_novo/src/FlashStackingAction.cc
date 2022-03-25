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
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include <sstream>
#include <string>
#include "G4UnitsTable.hh"

FlashStackingAction::FlashStackingAction()
    : G4UserStackingAction(), fScintillationCounter(0), fCerenkovCounter(0),
      fBremstralung(0), fFluo(0), Annihil_x(0), Annihil_y(0), Annihil_z(0) {
  std::ostringstream oss_1;
  oss_1 << "Electrons_time_distr" << G4Threading::G4GetThreadId() << ".csv";
  std::string filename_1 = oss_1.str();

  OpticFile.open(filename_1, std::ios_base::app);
}

FlashStackingAction::~FlashStackingAction() { OpticFile.close(); }

G4ClassificationOfNewTrack
FlashStackingAction::ClassifyNewTrack(const G4Track *aTrack) {

G4double globaltime;

  /*if(aTrack->GetDefinition()==G4Gamma::GammaDefinition()){
  if(aTrack->GetVolume()->GetLogicalVolume()->GetName() == "phantomLog"||aTrack->GetVolume()->GetLogicalVolume()->GetName() == "SupportLog"||aTrack->GetVolume()->GetLogicalVolume()->GetName() == "phantomLog_2"){

  if(aTrack->GetCreatorProcess()->GetProcessName()== "annihil")
  {
  Annihil_x = aTrack->GetPosition().x();
  Annihil_y = aTrack->GetPosition().y();
  Annihil_z = aTrack->GetPosition().z();
  globaltime= aTrack->GetGlobalTime();
  
 if (OpticFile.is_open()) {

      OpticFile // event ID ; X ; Y ; Z #all units are mm (if you want use G4BestUnits(x, "Length"))
         
          << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
          << "\t"
          << Annihil_x << "\t"
          << Annihil_y<< "\t"
          << Annihil_z <<"\t"<<globaltime<< G4endl;
    }
  }
  }}*/
  
if (aTrack->GetParentID() == 0){
globaltime=aTrack->GetGlobalTime();

if (OpticFile.is_open()) {

      OpticFile 
         
          << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
          << "\t"
          <<globaltime<< G4endl;
    }

} //this is a test
  

  return fUrgent;
}

void FlashStackingAction::NewStage()

{
  
}

void FlashStackingAction::PrepareNewEvent() {
  
}
