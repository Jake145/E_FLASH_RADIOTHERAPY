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

#ifndef Applicator_H
#define Applicator_H 1

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

class Applicator {
public:
  Applicator(G4VPhysicalVolume *);
  ~Applicator();

  void ConstructCollimator(G4VPhysicalVolume *);

  void FlashBeamLineVacuumSource();

  void FlashBeamLineTitaniumWindows();
  void FlashBeamLineFirstApplicator();

  void FlashBeamLineJunctions();
  void FlashBeamLineFinalApplicator();

  G4double initial_pos;
  G4double finalApplicatorXPositionFlash;
  G4double hightFinalApplicatorFlash;

private:
  G4VPhysicalVolume *motherPhys;
  void SetDefaultDimensions();
  void ConstructApplicator();

  G4VisAttributes *blue;
  G4VisAttributes *gray;
  G4VisAttributes *white;
  G4VisAttributes *red;
  G4VisAttributes *yellow;
  G4VisAttributes *green;
  G4VisAttributes *darkGreen;
  G4VisAttributes *darkOrange3;
  G4VisAttributes *skyBlue;
  G4VisAttributes *magenta;
  G4double Final_Attachment_X;

  G4double innerRadiusFirstApplicatorFlash, OuterRadiusFirstApplicatorFlash;
  G4Tubs *solidFirstApplicatorFlash;
  G4VPhysicalVolume *physiFirstApplicatorFlash;
  G4Material *firstApplicatorMaterialFlash;

  G4double innerRadiusFinalApplicatorFlash, OuterRadiusFinalApplicatorFlash;
  G4Tubs *solidFinalApplicatorFlash;
  G4VPhysicalVolume *physiFinalApplicatorFlash;
  G4Material *finalApplicatorMaterialFlash;

  G4Tubs *solidGiunz1FinalAppFlash;
  G4VPhysicalVolume *physiGiunz1FinalAppFlash;
  G4Material *Giunz1FinalAppMaterialFlash;

  G4Tubs *solidGiunz2FinalAppFlash;
  G4VPhysicalVolume *physiGiunz2FinalAppFlash;
  G4Material *Giunz2FinalAppMaterialFlash;

  // Junction 3 FINAL COLLIMATOR Flash
  G4Tubs *solidGiunz3FinalAppFlash;
  G4VPhysicalVolume *physiGiunz3FinalAppFlash;
  G4Material *Giunz3FinalAppMaterialFlash;

  // Junction 3 FINAL COLLIMATOR intFlash
  G4Cons *solidGiunz3FinalAppIntFlash;
  G4Material *Giunz3FinalAppMaterialIntFlash;
  G4VPhysicalVolume *physiGiunz3FinalAppIntFlash;

  // Junction 4 FINAL COLLIMATOR Flash
  G4Tubs *solidGiunz4FinalAppFlash;
  G4VPhysicalVolume *physiGiunz4FinalAppFlash;
  G4Material *Giunz4FinalAppMaterialFlash;

  // Junction 5 FINAL COLLIMATOR Flash
  G4Tubs *solidGiunz5FinalAppFlash;
  G4VPhysicalVolume *physiGiunz5FinalAppFlash;
  G4Material *Giunz5FinalAppMaterialFlash;

  // protector1
  G4Tubs *protector1;
  G4VPhysicalVolume *physiprotector1;
  G4Material *protector1material;

  // protector2
  G4Tubs *protector2;
  G4VPhysicalVolume *physiprotector2;
  G4Material *protector2material;

  // protector3
  G4Tubs *protector3;
  G4VPhysicalVolume *physiprotector3;
  G4Material *protector3material;

  // protector4
  G4Tubs *protector4;
  G4VPhysicalVolume *physiprotector4;
  G4Material *protector4material;

  // cover 1
  G4Tubs *cover1;
  G4VPhysicalVolume *physicover1Flash;
  G4Material *cover1material;

  // cover 2
  G4Tubs *cover2;
  G4VPhysicalVolume *physicover2Flash;
  G4Material *cover2material;

  // cover 3
  G4Tubs *cover3;
  G4VPhysicalVolume *physicover3Flash;
  G4Material *cover3material;

  //  Titanium Window
  G4Tubs *solidFTFlash;
  G4VPhysicalVolume *physiFTFlash;
  G4Material *FTFlashMaterialFlash;

  //  Vacuum Source
  G4Tubs *solidVSFlash;
  G4VPhysicalVolume *physiVSFlash;
  G4Material *VSFlashMaterialFlash;
};
#endif
