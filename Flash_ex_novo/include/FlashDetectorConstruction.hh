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
/// \file FlashDetectorConstruction.hh
/// \Classe per costruire il detector

#ifndef FlashDetectorConstruction_h
#define FlashDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4Material.hh"
#include "tls.hh"

#include "G4UserLimits.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class Applicator;
class VHEE_collimator;
class G4VSensitiveDetector;
class G4NistManager;
class G4Tubs;
class G4Box;
class G4Element;
class G4VisAttributes;

class FlashDetectorConstruction : public G4VUserDetectorConstruction {
public:
  G4VPhysicalVolume *physicalTreatmentRoom;
  G4LogicalVolume *logicTreatmentRoom;

  G4VPhysicalVolume *ConstructPhantom_Support(G4double CollPos, G4double Cx,
                                              G4double Cy, G4double Cz,
                                              G4double d, G4double Oz,
                                              G4bool plastic_bool);
  G4VPhysicalVolume *ConstructPhantom(G4double CollPos);
  G4VPhysicalVolume *BuildDetector(G4double dX, G4double dY, G4double dZ,
                                   G4double fPTFEThickness,
                                   G4double opticfiber_core_dx,
                                   G4bool plastic_bool);
  FlashDetectorConstruction();
  virtual ~FlashDetectorConstruction();

  virtual G4VPhysicalVolume *Construct();
  virtual void ConstructSDandField();

private:
  G4LogicalVolume *AirBox;
  G4Material *airNist;
  G4Material *TEFLON;
  Applicator *Collimator;
  VHEE_collimator *COLL;
  G4VPhysicalVolume *detector_physical;
  G4LogicalVolume *opticfiber_core_log;
  G4LogicalVolume *logicCryst;
  G4Box *phantom;
  void DefineMaterials();
  G4VisAttributes *skyBlue;
  G4VisAttributes *red;
  G4VisAttributes *blue;
  G4VisAttributes *green;

  G4LogicalVolume *phantomLogicalVolume;
  G4LogicalVolume *DetectorSupport;
  G4VPhysicalVolume *phant_phys;
  G4VPhysicalVolume *phantom_physical;
  G4UserLimits *fStepLimit;
  G4bool fCheckOverlaps;
  G4bool select_EJ212;
  G4bool Detector_builder;
  G4bool VHEE;
  G4double supp_coordinateX;
  G4double support_x;
  G4VPhysicalVolume *phys_cryst;
  G4double dX_;
  G4double dY_;
  G4double dZ_;
  G4double opticfiber_core_dx_;
  G4double fPTFEThickness_;
  G4Element *fN;
  G4Element *fO;

  G4Element *fC;
  G4Element *fH;
  G4NistManager *nist;
  G4Material *fPethylene2;
  G4Material *cryst_mat;
  G4Material *PMMA;
  G4Material *PMMA_optic;
  G4Material *PE;
  G4Material *ej212;
  G4Material *EJ212;
};

#endif
