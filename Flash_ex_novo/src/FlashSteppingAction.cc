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
  oss << "Kinetic_E_crystal_" << ThreadNumber << ".csv";
  std::string filename_1 = oss.str();

  KinEnFile.open(filename_1, std::ios_base::app);
}

FlashSteppingAction::~FlashSteppingAction() { KinEnFile.close(); }

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
    if (volumeName == "CrystalLV" && aStep->GetTrack()->GetTrackID() == 1) {
      G4String procName = postStep->GetProcessDefinedStep()->GetProcessName();

      G4double kineticEnergy = aStep->GetTrack()->GetKineticEnergy();

      if (KinEnFile.is_open()) {

        KinEnFile << "Incoming Energy"
                  << "\t" << eventid << "\t" << kineticEnergy << "\t" << trackID
                  << "\t" << procName << G4endl;
      }
    }

    if ((volumeName == "CrystalLV" && aStep->GetTrack()->GetDefinition() ==
                                          G4Electron::ElectronDefinition())) {

      if (KinEnFile.is_open()) {

        KinEnFile << "Incident Electron"
                  << "\t" << eventid << "\t" << trackID << "\t" << G4endl;
      }
    }
  }
}
