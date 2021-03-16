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
#include "tls.hh"
#include "G4Box.hh"
#include "G4UserLimits.hh"

class G4VPhysicalVolume; //definisco la classe del volume fisico
class G4LogicalVolume; //definisco la classe del volume logico
class Applicator80BeamLine;
class G4VSensitiveDetector;
class FlashSensitiveDetector;
class FlashDetectorConstruction : public G4VUserDetectorConstruction //classe della costruzione del detector e del fantoccio
{
  public:
      G4VPhysicalVolume* physicalTreatmentRoom;
  //FlashDetectorConstruction();
    void ConstructPhantom();
    FlashDetectorConstruction(); //costruttore ci passo il puntatore al mondo
    virtual ~FlashDetectorConstruction(); //distruttore

    virtual G4VPhysicalVolume* Construct(); //dichiaro il metodo per costruire il costruttore
    
  
    virtual void ConstructSDandField();// metodo per costruire il sensitive detector, necessario per contare gli eventi
    
 
    private:
     Applicator80BeamLine* Collimator;
	G4LogicalVolume* opticfiber_core_log;
	G4LogicalVolume* logicCryst;
	G4Box *phantom;
    void DefineMaterials(); //metodo per definire i materiali
	G4VisAttributes* skyBlue;
		G4VisAttributes* red;
				G4VisAttributes* green;
	//G4VPhysicalVolume* motherPhys;

	G4LogicalVolume* phantomLogicalVolume;
	G4VPhysicalVolume* phant_phys;
	G4UserLimits* fStepLimit; 
    G4bool  fCheckOverlaps; //booleano per vedere se vi sono degli overlap nella geometria (preso da B3a)
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

