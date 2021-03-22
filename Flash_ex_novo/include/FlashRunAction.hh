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
/// \file FlashRunAction.hh
/// \brief Definition of the FlashRunAction class

#ifndef FlashRunAction_h
#define FlashRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4RunManager.hh"
#include "globals.hh"
#include "G4Accumulable.hh"
class G4Run; // Dichiaro la classe della run

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class FlashRunAction : public G4UserRunAction //dichiaro il run action
{
  public:
    FlashRunAction(); //costruttore
    virtual ~FlashRunAction();//distruttore

    
    virtual void BeginOfRunAction(const G4Run*);//inizia la run
    virtual void   EndOfRunAction(const G4Run*);//finisce la run

    //void CountEvent()           { fGoodEvents += 1; };
    void SumEdep(G4double Edep) { fSumEdep += Edep; };  
    void CountPhotons_of(G4int count_of) {fcount_of +=count_of;};
    void CountCerenkov_of(G4int cerenkov_of) {fcerenkov_of+=cerenkov_of;};
    void CountScintillation_of(G4int scintillation_of) {fscintillation_of+=scintillation_of;};
    void CountPhotons_pd(G4int count_pd) {fcount_pd +=count_pd;};
    void CountCerenkov_pd(G4int cerenkov_pd) {fcerenkov_pd+=cerenkov_pd;};
    void CountScintillation_pd(G4int scintillation_pd) {fscintillation_pd+=scintillation_pd;};
    
 
private:
    //G4Accumulable<G4int>    fGoodEvents;
    G4Accumulable<G4double> fSumEdep;  
    G4Accumulable<G4int> fcount_of;
        G4Accumulable<G4int> fcerenkov_of;
            G4Accumulable<G4int> fscintillation_of;
            
    G4Accumulable<G4int> fcount_pd;
        G4Accumulable<G4int> fcerenkov_pd;
             G4Accumulable<G4int> fscintillation_pd;
   
};

#endif

