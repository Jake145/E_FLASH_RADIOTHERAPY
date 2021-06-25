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
//#include "G4Element.hh"
#include "G4Material.hh"
#include "tls.hh"
//#include "G4tubs.hh"
//#include "G4Box.hh"
#include "G4UserLimits.hh"

class G4VPhysicalVolume; // definisco la classe del volume fisico
class G4LogicalVolume;   // definisco la classe del volume logico
class Applicator80BeamLine;
class G4VSensitiveDetector;
class G4NistManager;
class G4Tubs;
class G4Box;
class G4Element;
class G4VisAttributes;
// class FSensitiveDetector;
// class OpticFiberSD;
// class PhotoDiodeSD;
class FlashDetectorConstruction
    : public G4VUserDetectorConstruction // classe della costruzione del
                                         // detector e del fantoccio
{
public:
  G4VPhysicalVolume *physicalTreatmentRoom;
G4LogicalVolume *logicTreatmentRoom;
  // FlashDetectorConstruction();
  G4VPhysicalVolume * ConstructPhantom(G4double Cx,G4double Cy,G4double Cz,G4double d,G4double Oz);
  FlashDetectorConstruction(); // costruttore ci passo il puntatore al mondo
  virtual ~FlashDetectorConstruction(); // distruttore

  virtual G4VPhysicalVolume *
  Construct(); // dichiaro il metodo per costruire il costruttore

  virtual void
  ConstructSDandField(); // metodo per costruire il sensitive detector,
                         // necessario per contare gli eventi

private:
  // FSensitiveDetector* sd_of_cryst;
  // OpticFiberSD* sd_of_of;
  // PhotoDiodeSD* sd_of_pd;
  G4LogicalVolume* AirBox;
   G4Material *airNist;
  G4Material* TEFLON;
  Applicator80BeamLine *Collimator;
  G4LogicalVolume *opticfiber_core_log;
  G4LogicalVolume *logicCryst;
  G4Box *phantom;
  void DefineMaterials(); // metodo per definire i materiali
  G4VisAttributes *skyBlue;
  G4VisAttributes *red;
    G4VisAttributes *blue;
  G4VisAttributes *green;
  // G4VPhysicalVolume* motherPhys;

  G4LogicalVolume *phantomLogicalVolume;
    G4LogicalVolume *DetectorSupport;
  G4VPhysicalVolume *phant_phys;
  G4VPhysicalVolume *phantom_physical;
  G4UserLimits *fStepLimit;
  G4bool fCheckOverlaps; // booleano per vedere se vi sono degli overlap nella
                         // geometria (preso da B3a)

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
  G4Tubs *fPhotocath;

  G4LogicalVolume *fPhotocath_log;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
