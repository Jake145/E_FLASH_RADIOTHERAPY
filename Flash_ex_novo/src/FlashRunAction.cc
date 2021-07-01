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
/// \file FlashRunAction.cc
/// \brief Implementation of the FlashRunAction class

#include "FlashRunAction.hh"
#include "FlashPrimaryGeneratorAction.hh"
//#include "FlashDetectorConstruction.hh"
// #include "FlashRun.hh"

#include "FlashAnalysis.hh"
#include "G4Accumulable.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashRunAction::FlashRunAction() : G4UserRunAction(), fSumEdep(0.),fEvents(0) {
  // Register accumulable to the accumulable manager
  G4AccumulableManager *accumulableManager = G4AccumulableManager::Instance();

  accumulableManager->RegisterAccumulable(fSumEdep);
    accumulableManager->RegisterAccumulable(fEvents);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashRunAction::~FlashRunAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void FlashRunAction::BeginOfRunAction(const G4Run* run)
void FlashRunAction::BeginOfRunAction(const G4Run *run) {
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  // oooooooooooooooooooOOOOOOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOOOOOOOOooooooooooooOOOOOOOOo
  /* //This is for creating the Ntuple. It is not working as expected
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  // Open an output file
  //analysisManager->OpenFile("FLASH_G4");

  // Creation of ntuple
  analysisManager->CreateNtuple("MyNtuple", "Kinetic Energy in Detector");
  // X = D in CreateNtupleXColumn stands for G4double (I,F,D,S)
  analysisManager->CreateNtupleIColumn("Particle_ID");
  analysisManager->CreateNtupleDColumn("Kinetic_Energy");
  analysisManager->CreateNtupleIColumn("Event");
  analysisManager->CreateNtupleSColumn("LogicalVolume");
  analysisManager->FinishNtuple(); */
  // OOOOOOOOOOOOOOOOOOOOOoooooooooooooOOOoooOOooooooooooooooooooOOOOOOOOOOooooooooOOoooOOOoooOOoooo
  // reset accumulables to their initial values
  G4AccumulableManager *accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void FlashRunAction::EndOfRunAction(const G4Run* run)
void FlashRunAction::EndOfRunAction(const G4Run *run) {
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0)
    return;

  // Merge accumulables
  G4AccumulableManager *accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Print results
  //
  if (IsMaster()) {
    G4cout << G4endl
           << "--------------------End of Global Run-----------------------"
           << G4endl << "  The run was " << nofEvents << " events ";
  } else {
    G4cout << G4endl
           << "--------------------End of Local Run------------------------"
           << G4endl << "  The run was " << nofEvents << " events ";
  }
  G4cout

      << " Total Energy in crystal : "
      << G4BestUnit(fSumEdep.GetValue(), "Energy") << "\n"<<"Total events: " <<fEvents.GetValue()<<G4endl
      << "------------------------------------------------------------"
      << G4endl << G4endl;
}
