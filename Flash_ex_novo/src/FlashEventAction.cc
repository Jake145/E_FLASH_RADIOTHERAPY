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

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <string>
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashEventAction::FlashEventAction(FlashRunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),fCollID_cryst(-1),collectionID_cryst(-1),collectionID_of(-1),
  collectionID_pd(-1){
  fAbsorptionCount = 0;
  fBoundaryAbsorptionCount = 0;
  std::ostringstream oss;
oss << "Kinetic_E_crystal_" << G4Threading::G4GetThreadId() << ".csv" ;
std::string filename = oss.str();

  KinEnFile.open(filename,std::ios_base::app);
  }
   //inizializzo il costruttore dando gli argomenti appropriati alle funzioni
 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashEventAction::~FlashEventAction() //distruttore
{KinEnFile.close();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FlashEventAction::BeginOfEventAction(const G4Event*)
{    
  

  }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void FlashEventAction::EndOfEventAction(const G4Event* evt)
void FlashEventAction::EndOfEventAction(const G4Event* evt)
{     
G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  G4SDManager* SDMan = G4SDManager::GetSDMpointer();  
    fCollID_cryst   = SDMan->GetCollectionID("crystalSD/edep");
     
   G4THitsMap<G4double>* evtMap = 
                     (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst));
                     
  std::map<G4int,G4double*>::iterator itr;
        G4double edep = 0.;           
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++){
    ///G4int copyNb  = (itr->first);
     edep = *(itr->second);};
    if (edep > 0.) fRunAction->SumEdep(edep);
    //Look at optic events in crystal
    
    G4cout << ">>> Event " << evt->GetEventID() << G4endl;

 
G4SDManager * SDman_cryst = G4SDManager::GetSDMpointer();

  if(collectionID_cryst<0){

    G4String colNam;

    collectionID_cryst = SDman_cryst->GetCollectionID("CrystSD/FlashHitsCollection");
    
    G4HCofThisEvent * HCE_cryst = evt->GetHCofThisEvent();

  FlashHitsCollection* HC_cryst = 0;

  if(HCE_cryst)

  {

    HC_cryst = (FlashHitsCollection*)(HCE_cryst->GetHC(collectionID_cryst));

  }

 

  if ( HC_cryst ) {

    int n_hit = HC_cryst->entries();

    for ( int i = 0 ; i < n_hit; i++){

      
      G4double      energy   = (*HC_cryst)[i]->GetEnergy();
if(KinEnFile.is_open()){
 
 KinEnFile<< evt->GetEventID()<< "\t" << energy <<G4endl;
 }

    }

  } 
}


G4SDManager * SDman_of = G4SDManager::GetSDMpointer();

  if(collectionID_of<0){

    G4String colNam_of;

    collectionID_of = SDman_of->GetCollectionID("OpticFiberSD/OpticFiberHitsCollection");
    
    G4HCofThisEvent * HCE_of = evt->GetHCofThisEvent();

  OpticFiberHitsCollection* HC_of = 0;

  if(HCE_of)

  {

    HC_of = (OpticFiberHitsCollection*)(HCE_of->GetHC(collectionID_of));

  }

 

  if ( HC_of ) {

    int n_hit_of = HC_of->entries();

    for ( int i = 0 ; i < n_hit_of; i++){

      
      G4int      photons_of   = (*HC_of)[i]->GetPhotonCount();
      G4int 	 cerenkov_of = (*HC_of)[i]->GetCerenkovCount();
      G4int	 scintillation_of = (*HC_of)[i]->GetScintillationCount();
      
      if (photons_of > 0) fRunAction->CountPhotons_of(photons_of);
      if (cerenkov_of > 0) fRunAction->CountCerenkov_of(cerenkov_of);
      if (scintillation_of> 0) fRunAction->CountScintillation_of(scintillation_of);
    }

  } }



G4SDManager * SDman_pd = G4SDManager::GetSDMpointer();

  if(collectionID_pd<0){

    G4String colNam_pd;

    collectionID_of = SDman_pd->GetCollectionID("PhotoDiodeSD/PhotoDiodeHitsCollection");
    
    G4HCofThisEvent * HCE_pd = evt->GetHCofThisEvent();

  PhotoDiodeHitsCollection* HC_pd = 0;

  if(HCE_pd)

  {

    HC_pd = (PhotoDiodeHitsCollection*)(HCE_pd->GetHC(collectionID_pd));

  }

 

  if ( HC_pd ) {

    int n_hit_pd = HC_pd->entries();

    for ( int i = 0 ; i < n_hit_pd; i++){

      
      G4int      photons_pd   = (*HC_pd)[i]->GetPhotonCount_photo();
      G4int 	 cerenkov_pd = (*HC_pd)[i]->GetCerenkovCount_photo();
      G4int	 scintillation_pd = (*HC_pd)[i]->GetScintillationCount_photo();
      
      if (photons_pd > 0) fRunAction->CountPhotons_pd(photons_pd);
      if (cerenkov_pd > 0) fRunAction->CountCerenkov_pd(cerenkov_pd);
      if (scintillation_pd > 0) fRunAction->CountScintillation_pd(scintillation_pd);
    }

  } }
}
  
  


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
