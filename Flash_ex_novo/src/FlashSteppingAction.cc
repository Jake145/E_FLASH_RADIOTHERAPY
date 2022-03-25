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
#include "FlashDetectorConstruction.hh"
#include "FlashEventAction.hh"
#include "FlashRunAction.hh"
#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include <sstream>
#include <string>

FlashSteppingAction::FlashSteppingAction(FlashEventAction *eventAction)
    : G4UserSteppingAction(), ThreadNumber(G4Threading::G4GetThreadId()),
      fEventAction(eventAction) {

  std::ostringstream oss;
  oss << "Annihil_s_surf_det_1_" << ThreadNumber << ".csv";
  
  std::ostringstream oss_2;
  oss_2 << "Annihil_s_surf_det_2_" << ThreadNumber << ".csv";
  
  std::string filename_1 = oss.str();
  
  std::string filename_2 = oss_2.str();

  KinEnFile.open(filename_1, std::ios_base::app);
  
  KinEnFile_2.open(filename_2, std::ios_base::app);
}

FlashSteppingAction::~FlashSteppingAction() { KinEnFile.close(); KinEnFile_2.close(); }

void FlashSteppingAction::UserSteppingAction(const G4Step *aStep) {

  G4int eventid =
      G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  G4StepPoint *postStep = aStep->GetPostStepPoint();
  G4StepPoint *preStep = aStep->GetPreStepPoint();
  G4int trackID = aStep->GetTrack()->GetTrackID();

  //===========================INCIDENT
  //PARTICLES===============================================

  //==============================PRIMARIES ENERGY
  //SPECTRUM=======================================
  if (postStep->GetStepStatus() == fGeomBoundary) {

    G4String volumeName =
        postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
    G4String prevolumeName =
        preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
    if (aStep->GetTrack()->GetDefinition()==G4Gamma::GammaDefinition() && volumeName == "s_surf_log" && aStep->GetTrack()->GetCreatorProcess()->GetProcessName()== "annihil") {

      G4double kineticEnergy = aStep->GetTrack()->GetKineticEnergy();
      G4double kin_en_2= aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy(); //test
      G4double mom_axial = aStep->GetTrack()->GetMomentumDirection().x();
      G4double mom_y = aStep->GetTrack()->GetMomentumDirection().y();
      G4double mom_z = aStep->GetTrack()->GetMomentumDirection().z();
      G4double pos_x = aStep->GetTrack()->GetPosition().x();
      G4double pos_y = aStep->GetTrack()->GetPosition().y();
      G4double pos_z = aStep->GetTrack()->GetPosition().z();
      G4double time = aStep->GetTrack()->GetGlobalTime();
      if (KinEnFile.is_open()) {

        KinEnFile << eventid << "\t" << kineticEnergy << "\t" <<kin_en_2 << "\t" <<mom_axial<< "\t" << mom_y<< "\t" << mom_z<< "\t" <<pos_x<< "\t" << pos_y<< "\t" << pos_z<<"\t" << time<<G4endl;
       
      }
    }

 if (aStep->GetTrack()->GetDefinition()==G4Gamma::GammaDefinition() && volumeName == "s_surf_log_2" && aStep->GetTrack()->GetCreatorProcess()->GetProcessName()== "annihil") {

      G4double kineticEnergy = aStep->GetTrack()->GetKineticEnergy();
      G4double kin_en_2= aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy(); //test
      G4double mom_axial = aStep->GetTrack()->GetMomentumDirection().x();
      G4double mom_y = aStep->GetTrack()->GetMomentumDirection().y();
      G4double mom_z = aStep->GetTrack()->GetMomentumDirection().z();
      G4double pos_x = aStep->GetTrack()->GetPosition().x();
      G4double pos_y = aStep->GetTrack()->GetPosition().y();
      G4double pos_z = aStep->GetTrack()->GetPosition().z();
      G4double time = aStep->GetTrack()->GetGlobalTime();
      if (KinEnFile_2.is_open()) {

        KinEnFile_2 << eventid << "\t" << kineticEnergy << "\t" <<kin_en_2 << "\t" <<mom_axial<< "\t" << mom_y<< "\t" << mom_z<< "\t" <<pos_x<< "\t" << pos_y<< "\t" << pos_z<<"\t" << time<<G4endl;
       
      }
    }
    
  }
}
