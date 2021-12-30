//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Appaboration.  It is provided  under  the terms  and *
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
//

#ifndef Collimator_H
#define Collimator_H 1

#include "FlashDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

class G4VPhysicalVolume;

class VHEE_collimator {
public:
  VHEE_collimator(G4VPhysicalVolume *);
  ~VHEE_collimator();

  void ConstructCollimator();


  

private:
G4Material *Bromuro_SI;
  G4VPhysicalVolume *motherPhys;
  G4double X_center;
   G4double y_center;
    G4double z_center;
     G4double radious_field;
      G4int number_of_leafs;
  void SetMaterial();
  void ConstructColl();

  G4VisAttributes *blue;
  

  G4LogicalVolume * LeafLV;
};
#endif
