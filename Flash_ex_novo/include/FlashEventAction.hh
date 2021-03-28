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
/// \file FlashEventAction.hh
/// \brief Definition of the FlashEventAction class

#ifndef FlashEventAction_h
#define FlashEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class FlashRunAction; // definisco la classe Run Action

/// Event action class
///

class FlashEventAction
    : public G4UserEventAction // definisco la Event Action che eredita da
                               // //G4UserEventAction
{
public:
  FlashEventAction(FlashRunAction *runAction); // definisco il costruttore che
                                               // ha in argomento un
                                               // //puntatore runAction
  virtual ~FlashEventAction(); // definisco il distruttore

  virtual void BeginOfEventAction(
      const G4Event *event); // definisco la funzione che inizializza
                             // //l'evento
  virtual void
  EndOfEventAction(const G4Event *event); // similmente quello che lo termina

private:
  FlashRunAction *fRunAction;
  G4int fCollID_cryst;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
