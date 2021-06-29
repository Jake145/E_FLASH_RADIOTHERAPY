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
/// \file FlashEventAction.cc
/// \brief Implementation of the FlashEventAction class

#include "FlashEventAction.hh"
#include "FlashRunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashEventAction::FlashEventAction(FlashRunAction *runAction)
    : G4UserEventAction(), fRunAction(runAction), fCollID_cryst(-1) {}
// inizializzo il costruttore dando gli argomenti appropriati alle funzioni

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashEventAction::~FlashEventAction() // distruttore
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FlashEventAction::BeginOfEventAction(const G4Event *) {
  G4cout << "Starting event: "
         << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void FlashEventAction::EndOfEventAction(const G4Event* evt)
void FlashEventAction::EndOfEventAction(const G4Event *evt) {
 /* G4HCofThisEvent *HCE = evt->GetHCofThisEvent();
  if (!HCE)
    return;
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  fCollID_cryst = SDMan->GetCollectionID("crystalSD/edep");

  G4THitsMap<G4double> *evtMap =
      (G4THitsMap<G4double> *)(HCE->GetHC(fCollID_cryst));

  std::map<G4int, G4double *>::iterator itr;
  G4double edep = 0.;
  G4int Number_of_events=0;
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    /// G4int copyNb  = (itr->first);
    edep = *(itr->second);
  };
  if (edep > 0.)
    fRunAction->SumEdep(edep);

  G4cout << "Ending event: "
         << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
         << G4endl;*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
