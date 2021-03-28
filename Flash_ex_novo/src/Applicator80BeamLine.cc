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
// This is the first version of Flash, a Geant4-based application
//
//
//////////////////////////////////////////////////////////////////////////////////////////////

#include "Applicator80BeamLine.hh"
#include "FlashDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistElementBuilder.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4RunManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
//#include "Applicator80BeamLineMessenger.hh"

Applicator80BeamLine::Applicator80BeamLine(
    G4VPhysicalVolume *physicalTreatmentRoom)
    : // physicalTreatmentRoom(0),
      motherPhys(physicalTreatmentRoom),

      solidFirstApplicatorFlash(0), physiFirstApplicatorFlash(0),
      solidFinalApplicatorFlash(0), physiFinalApplicatorFlash(0),

      solidGiunz1FinalAppFlash(0), physiGiunz1FinalAppFlash(0),

      solidGiunz2FinalAppFlash(0), physiGiunz2FinalAppFlash(0),

      solidGiunz3FinalAppFlash(0), physiGiunz3FinalAppFlash(0),

      solidGiunz3FinalAppIntFlash(0), physiGiunz3FinalAppIntFlash(0),

      solidGiunz4FinalAppFlash(0), physiGiunz4FinalAppFlash(0),

      solidGiunz5FinalAppFlash(0), physiGiunz5FinalAppFlash(0),

      protector1(0), physiprotector1(0),

      protector2(0), physiprotector2(0),

      protector3(0), physiprotector3(0),

      protector4(0), physiprotector4(0),

      cover1(0), physicover1Flash(0),

      cover2(0), physicover2Flash(0),

      cover3(0), physicover3Flash(0),

      solidFTFlash(0), physiFTFlash(0) {
  ConstructCollimator(motherPhys);
  // Messenger to change parameters of the applicator80BeamLine geometry
  // applicatorMessenger = new Applicator80BeamLineMessenger(this);
}

Applicator80BeamLine::~Applicator80BeamLine() {}

void Applicator80BeamLine::ConstructCollimator(G4VPhysicalVolume *) {
  // Sets default geometry and materials
  SetDefaultDimensions();

  // Construct the whole Applicator Beam Line
  ConstructApplicator80BeamLine();

  // FlashDetectorConstruction builds ONLY the phantom and the detector with its
  // associated Geometry
  // flashDetectorConstruction = new
  // FlashDetectorConstruction(physicalTreatmentRoom);

  // return physicalTreatmentRoom;
}

// DEFAULT MATERIAL ARE ALSO PROVIDED
// and COLOURS ARE ALSO DEFINED
// ----------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////
void Applicator80BeamLine::SetDefaultDimensions() {

  // Set of coulors that can be used
  white = new G4VisAttributes(G4Colour());
  white->SetVisibility(true);
  // white -> SetForceSolid(true);

  blue = new G4VisAttributes(G4Colour(0., 0., 1.));
  blue->SetVisibility(true);
  // blue -> SetForceSolid(true);

  gray = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  gray->SetVisibility(true);
  // gray-> SetForceSolid(true);

  red = new G4VisAttributes(G4Colour(1., 0., 0.));
  red->SetVisibility(true);
  // red-> SetForceSolid(true);

  yellow = new G4VisAttributes(G4Colour(1., 1., 0.));
  yellow->SetVisibility(true);
  // yellow-> SetForceSolid(true);

  green = new G4VisAttributes(G4Colour(25 / 255., 255 / 255., 25 / 255.));
  green->SetVisibility(true);
  // green -> SetForceSolid(true);

  darkGreen = new G4VisAttributes(G4Colour(0 / 255., 100 / 255., 0 / 255.));
  darkGreen->SetVisibility(true);
  // darkGreen -> SetForceSolid(true);

  darkOrange3 =
      new G4VisAttributes(G4Colour(205 / 255., 102 / 255., 000 / 255.));
  darkOrange3->SetVisibility(true);
  // darkOrange3 -> SetForceSolid(true);

  skyBlue = new G4VisAttributes(G4Colour(135 / 255., 206 / 255., 235 / 255.));
  skyBlue->SetVisibility(true);
  // skyBlue -> SetForceSolid(true);
  magenta = new G4VisAttributes(G4Colour(255 / 255., 0 / 255., 255 / 255.));
  magenta->SetVisibility(true);

  // Geometry FIRST APPLICATOR DEFAULTS

  G4double defaultOuterRadiusFirstApplicatorFlash = 55. * mm;
  OuterRadiusFirstApplicatorFlash = defaultOuterRadiusFirstApplicatorFlash;

  G4double defaultinnerRadiusFirstApplicatorFlash = 50. * mm;
  innerRadiusFirstApplicatorFlash = defaultinnerRadiusFirstApplicatorFlash;

  // Geometry FINAL APPLICATOR DEFAULTS

  G4double defaultOuterRadiusFinalApplicatorFlash = 55. * mm;
  OuterRadiusFinalApplicatorFlash = defaultOuterRadiusFinalApplicatorFlash;

  G4double defaultinnerRadiusFinalApplicatorFlash = 50. * mm;
  innerRadiusFinalApplicatorFlash = defaultinnerRadiusFinalApplicatorFlash;

  // DEFAULT DEFINITION OF THE MATERIALS
  // All elements and compound definition follows the NIST database

  // ELEMENTS
  G4bool isotopes = false;
  G4Material *aluminumNist =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Al", isotopes);
  // G4Material* tantalumNist =
  // G4NistManager::Instance()->FindOrBuildMaterial("G4_Ta", isotopes);
  // G4Material* copperNistAsMaterial =
  // G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu", isotopes);
  G4Element *zincNist = G4NistManager::Instance()->FindOrBuildElement("Zn");
  G4Element *copperNist = G4NistManager::Instance()->FindOrBuildElement("Cu");

  // COMPOUND
  // G4Material* airNist =
  // G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  // G4Material* kaptonNist =
  // G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON", isotopes);
  G4Material *galacticNist =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic", isotopes);
  G4Material *PMMANist =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);
  // G4Material* mylarNist =
  // G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR", isotopes);
  G4Material *titanioNist =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Ti", isotopes);

  G4double d;            // Density
  G4int nComponents;     // Number of components
  G4double fractionmass; // Fraction in mass of an element in a material

  d = 8.40 * g / cm3; // brass
  nComponents = 2;
  G4Material *brass = new G4Material("Brass", d, nComponents);
  brass->AddElement(zincNist, fractionmass = 30 * perCent);
  brass->AddElement(copperNist, fractionmass = 70 * perCent);

  // MATERIAL ASSIGNMENT

  // Material of the FINAL APPLICATOR Flash
  firstApplicatorMaterialFlash = PMMANist;

  // Material of the FINAL APPLICATOR Flash
  finalApplicatorMaterialFlash = PMMANist;

  // monitor 1 FINAL APPLICATOR Flash
  Giunz1FinalAppMaterialFlash = aluminumNist;

  // monitor 2 FINAL APPLICATOR Flash
  Giunz2FinalAppMaterialFlash = aluminumNist;

  // part of monitor
  cover1material = aluminumNist;

  // part of monitor
  cover2material = aluminumNist;

  // part of applicator
  cover3material = PMMANist;

  // all protectors encase the monitors and are made of pmma
  protector1material = PMMANist;
  protector2material = PMMANist;
  protector3material = PMMANist;
  protector4material = PMMANist;

  // monitor 3 FINAL APPLICATOR Int Flash
  Giunz3FinalAppMaterialFlash = aluminumNist;

  // monitor 4 FINAL APPLICATOR Flash
  Giunz4FinalAppMaterialFlash = aluminumNist;

  // Junction 5 FINAL APPLICATOR Flash
  Giunz5FinalAppMaterialFlash = PMMANist;

  // Superior Final Part Monitor Chambers Material
  FTFlashMaterialFlash = titanioNist;

  // Vacuum Source
  VSFlashMaterialFlash = galacticNist;
}

/////////////////////////////////////////////////////////////////////////////
void Applicator80BeamLine::ConstructApplicator80BeamLine() {

  // Components of the Applicator Beam Line

  FlashBeamLineVacuumSource();
  FlashBeamLineTitaniumWindows();
  FlashBeamLineFirstApplicator();

  FlashBeamLineJunctions();
  FlashBeamLineFinalApplicator();
}

void Applicator80BeamLine::FlashBeamLineVacuumSource() {
  // ---------------------------------------------------------------//
  //                     Vacuum Source                             //
  // ---------------------------------------------------------------//

  G4double phi1 = 90. * deg;

  G4RotationMatrix rm1;
  rm1.rotateY(phi1);

  const G4double outRadiusVSFlash = 44.75 * mm;
  const G4double innRadiusVSFlash = 0. * mm;
  const G4double hightVSFlash = 0.5 * mm;
  const G4double startAngleVSFlash = 0. * deg;
  const G4double spanningAngleVSFlash = 360. * deg;
  const G4double XPositionVSFlash = -1000.0 * mm;

  solidVSFlash =
      new G4Tubs("VSFlash", innRadiusVSFlash, outRadiusVSFlash, hightVSFlash,
                 startAngleVSFlash, spanningAngleVSFlash);

  G4LogicalVolume *logVSFlash = new G4LogicalVolume(
      solidVSFlash, VSFlashMaterialFlash, "VSFlash", 0, 0, 0);

  physiVSFlash = new G4PVPlacement(
      G4Transform3D(rm1, G4ThreeVector((XPositionVSFlash), 0., 0.)), "VSFlash",
      logVSFlash, motherPhys, false, 0);

  logVSFlash->SetVisAttributes(green);
}
/////////
void Applicator80BeamLine::FlashBeamLineTitaniumWindows() {
  // ---------------------------------------------------------------//
  //                     Titanium Window                        //
  // ---------------------------------------------------------------//
  // with just this piece ssd=1.6cm
  G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);

  const G4double outRadiusFTFlash = 44.75 * mm;
  const G4double innRadiusFTFlash = 8.5 * mm;
  const G4double hightFTFlash = 0.03 * mm;
  const G4double startAngleFTFlash = 0. * deg;
  const G4double spanningAngleFTFlash = 360. * deg;
  const G4double XPositionFTFlash = -999.47 * mm;

  solidFTFlash =
      new G4Tubs("FTFlash", innRadiusFTFlash, outRadiusFTFlash, hightFTFlash,
                 startAngleFTFlash, spanningAngleFTFlash);

  G4LogicalVolume *logFTFlash = new G4LogicalVolume(
      solidFTFlash, FTFlashMaterialFlash, "FTFlash", 0, 0, 0);

  physiFTFlash = new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPositionFTFlash), 0., 0.)), "FTFlash",
      logFTFlash, motherPhys, false, 0);

  logFTFlash->SetVisAttributes(yellow);
}
void Applicator80BeamLine::FlashBeamLineFirstApplicator() {

  // -----------------------//
  // FIRST APPLICATOR Flash  //
  //------------------------//

  // const G4double outRadiusFirstApplicatorFlash = 45. *mm;
  // const G4double innRadiusFirstApplicatorFlash = 40. *mm;
  const G4double hightFirstApplicatorFlash = 100. * mm;
  const G4double startAngleFirstApplicatorFlash = 0. * deg;
  const G4double spanningAngleFirstApplicatorFlash = 360. * deg;
  const G4double firstApplicatorXPositionFlash = -799.44 * mm;

  G4double phi6 = 90. * deg;

  G4RotationMatrix rm6;
  rm6.rotateY(phi6);

  solidFirstApplicatorFlash = new G4Tubs(
      "FirstApplicatorFlash", innerRadiusFirstApplicatorFlash,
      OuterRadiusFirstApplicatorFlash, hightFirstApplicatorFlash,
      startAngleFirstApplicatorFlash, spanningAngleFirstApplicatorFlash);

  G4LogicalVolume *logFirstApplicatorFlash = new G4LogicalVolume(
      solidFirstApplicatorFlash, firstApplicatorMaterialFlash,
      "FirstApplicatorFlash", 0, 0, 0);

  physiFirstApplicatorFlash = new G4PVPlacement(
      G4Transform3D(rm6,
                    G4ThreeVector((firstApplicatorXPositionFlash), 0., 0.)),
      "FirstApplicatorFlash", logFirstApplicatorFlash, motherPhys, false, 0);

  //  logFirstApplicatorFlash ->
  //  SetVisAttributes(G4VisAttributes::GetInvisible());
  logFirstApplicatorFlash->SetVisAttributes(magenta);
}

void Applicator80BeamLine::FlashBeamLineJunctions() {

  G4double phi5 = 90. * deg;

  G4RotationMatrix rm5;
  rm5.rotateY(phi5);
  // --------------------------------- //
  // Junction 5 FIRST APPLICATOR Flash //
  // --------------------------------- //

  const G4double outRadiusGiunz5FinalAppFlash = 60 * mm;
  const G4double innRadiusGiunz5FinalAppFlash = 55 * mm;
  const G4double hightGiunz5FinalAppFlash = 10. * mm;
  const G4double startAngleGiunz5FinalAppFlash = 0. * deg;
  const G4double spanningAngleGiunz5FinalAppFlash = 360. * deg;
  const G4double Giunz5FinalAppXPositionFlash = -689.44 * mm;

  solidGiunz5FinalAppFlash = new G4Tubs(
      "Giunz5FinalAppFlash", innRadiusGiunz5FinalAppFlash,
      outRadiusGiunz5FinalAppFlash, hightGiunz5FinalAppFlash,
      startAngleGiunz5FinalAppFlash, spanningAngleGiunz5FinalAppFlash);

  G4LogicalVolume *logGiunz5FinalAppFlash =
      new G4LogicalVolume(solidGiunz5FinalAppFlash, Giunz5FinalAppMaterialFlash,
                          "Giunz5FinalAppFlash", 0, 0, 0);

  physiGiunz5FinalAppFlash = new G4PVPlacement(
      G4Transform3D(rm5, G4ThreeVector((Giunz5FinalAppXPositionFlash), 0., 0.)),
      "Giunz5FinalAppFlash", logGiunz5FinalAppFlash, motherPhys, false, 0);

  logGiunz5FinalAppFlash->SetVisAttributes(yellow);

  // --------------------------------- //
  // Junction 4 FINAL APPLICATOR Flash //
  // --------------------------------- //

  const G4double innRadiusGiunz4FinalAppFlash = 22.25 * mm;
  const G4double outRadiusGiunz4FinalAppFlash = 27.25 * mm;
  const G4double hightGiunz4FinalAppFlash = 9.35 * mm;
  const G4double startAngleGiunz4FinalAppFlash = 0. * deg;
  const G4double spanningAngleGiunz4FinalAppFlash = 360. * deg;
  const G4double Giunz4FinalAppXPositionFlash = -908.79 * mm;

  solidGiunz4FinalAppFlash = new G4Tubs(
      "Giunz4FinalAppFlash", innRadiusGiunz4FinalAppFlash,
      outRadiusGiunz4FinalAppFlash, hightGiunz4FinalAppFlash,
      startAngleGiunz4FinalAppFlash, spanningAngleGiunz4FinalAppFlash);

  G4LogicalVolume *logGiunz4FinalAppFlash =
      new G4LogicalVolume(solidGiunz4FinalAppFlash, Giunz4FinalAppMaterialFlash,
                          "Giunz4FinalAppFlash", 0, 0, 0);

  physiGiunz4FinalAppFlash = new G4PVPlacement(
      G4Transform3D(rm5, G4ThreeVector((Giunz4FinalAppXPositionFlash), 0., 0.)),
      "Giunz4FinalAppFlash", logGiunz4FinalAppFlash, motherPhys, false, 0);

  logGiunz4FinalAppFlash->SetVisAttributes(blue);
  // protector4
  const G4double outRadiusprotector4Flash = 32 * mm;
  const G4double innRadiusprotector4Flash = 27.25 * mm;
  const G4double hightprotector4Flash = 9.35 * mm;
  const G4double startAngleprotector4Flash = 0. * deg;
  const G4double spanningAngleprotector4Flash = 360. * deg;
  const G4double protector4XPositionFlash = -908.79 * mm;

  protector4 =
      new G4Tubs("protector4", innRadiusprotector4Flash,
                 outRadiusprotector4Flash, hightprotector4Flash,
                 startAngleprotector4Flash, spanningAngleprotector4Flash);

  G4LogicalVolume *logprotector4Flash = new G4LogicalVolume(
      protector4, protector4material, "protector4", 0, 0, 0);

  physiprotector4 = new G4PVPlacement(
      G4Transform3D(rm5, G4ThreeVector((protector4XPositionFlash), 0., 0.)),
      "protector4", logprotector4Flash, motherPhys, false, 0);

  // --------------------------------- //
  // monitor 3 FINAL APPLICATOR Flash //
  // --------------------------------- //

  const G4double outRadiusGiunz3FinalAppFlash = 18.25 * mm;
  const G4double innRadiusGiunz3FinalAppFlash = 13.75 * mm;
  const G4double hightGiunz3FinalAppFlash = 4.65 * mm;
  const G4double startAngleGiunz3FinalAppFlash = 0. * deg;
  const G4double spanningAngleGiunz3FinalAppFlash = 360. * deg;
  const G4double Giunz3FinalAppXPositionFlash = -922.79 * mm;

  solidGiunz3FinalAppFlash = new G4Tubs(
      "Giunz3FinalAppFlash", innRadiusGiunz3FinalAppFlash,
      outRadiusGiunz3FinalAppFlash, hightGiunz3FinalAppFlash,
      startAngleGiunz3FinalAppFlash, spanningAngleGiunz3FinalAppFlash);

  G4LogicalVolume *logicsolidGiunz3FinalAppFlash =
      new G4LogicalVolume(solidGiunz3FinalAppFlash, Giunz3FinalAppMaterialFlash,
                          "Giunz3FinalAppFlash", 0, 0, 0);

  physiGiunz3FinalAppFlash = new G4PVPlacement(
      G4Transform3D(rm5, G4ThreeVector((Giunz3FinalAppXPositionFlash), 0., 0.)),
      "Giunz3FinalAppFlash", logicsolidGiunz3FinalAppFlash, motherPhys, false,
      0);

  logicsolidGiunz3FinalAppFlash->SetVisAttributes(yellow);

  // protector3
  const G4double outRadiusprotector3Flash = 25 * mm;
  const G4double innRadiusprotector3Flash = 18.25 * mm;
  const G4double hightprotector3Flash = 4.65 * mm;
  const G4double startAngleprotector3Flash = 0. * deg;
  const G4double spanningAngleprotector3Flash = 360. * deg;
  const G4double protector3XPositionFlash = -922.79 * mm;

  protector3 =
      new G4Tubs("protector3", innRadiusprotector3Flash,
                 outRadiusprotector3Flash, hightprotector3Flash,
                 startAngleprotector3Flash, spanningAngleprotector3Flash);

  G4LogicalVolume *logprotector3Flash = new G4LogicalVolume(
      protector3, protector3material, "protector3", 0, 0, 0);

  physiprotector3 = new G4PVPlacement(
      G4Transform3D(rm5, G4ThreeVector((protector3XPositionFlash), 0., 0.)),
      "protector3", logprotector3Flash, motherPhys, false, 0);

  logprotector3Flash->SetVisAttributes(gray);
  // cover disk between monitor piece 3 and 4

  const G4double innRadiuscover2Flash = 18.25 * mm;
  const G4double outRadiuscover2Flash = 22.25 * mm;
  const G4double hightcover2Flash = 0.01 * mm;
  const G4double startAnglecover2Flash = 0. * deg;
  const G4double spanningAnglecover2Flash = 360. * deg;
  const G4double cover2XPositionFlash = -918.14 * mm;

  cover2 = new G4Tubs("cover2", innRadiuscover2Flash, outRadiuscover2Flash,
                      hightcover2Flash, startAnglecover2Flash,
                      spanningAnglecover2Flash);

  G4LogicalVolume *logcover2 =
      new G4LogicalVolume(cover2, cover2material, "cover2", 0, 0, 0);

  physicover2Flash = new G4PVPlacement(
      G4Transform3D(rm5, G4ThreeVector((cover2XPositionFlash), 0., 0.)),
      "cover2", logcover2, motherPhys, false, 0);

  logcover2->SetVisAttributes(gray);
  // --------------------------------- //
  // monitor 2 FINAL APPLICATOR Flash //
  // --------------------------------- //

  const G4double innRadiusGiunz2FinalAppFlash = 15.0 * mm;
  const G4double outRadiusGiunz2FinalAppFlash = 20.0 * mm;
  const G4double hightGiunz2FinalAppFlash = 30 * mm;
  const G4double startAngleGiunz2FinalAppFlash = 0. * deg;
  const G4double spanningAngleGiunz2FinalAppFlash = 360. * deg;
  const G4double Giunz2FinalAppXPositionFlash = -957.44 * mm;

  solidGiunz2FinalAppFlash = new G4Tubs(
      "Giunz2FinalAppFlash", innRadiusGiunz2FinalAppFlash,
      outRadiusGiunz2FinalAppFlash, hightGiunz2FinalAppFlash,
      startAngleGiunz2FinalAppFlash, spanningAngleGiunz2FinalAppFlash);

  G4LogicalVolume *logGiunz2FinalAppFlash =
      new G4LogicalVolume(solidGiunz2FinalAppFlash, Giunz2FinalAppMaterialFlash,
                          "Giunz2FinalAppFlash", 0, 0, 0);

  physiGiunz2FinalAppFlash = new G4PVPlacement(
      G4Transform3D(rm5, G4ThreeVector((Giunz2FinalAppXPositionFlash), 0., 0.)),
      "Giunz2FinalAppFlash", logGiunz2FinalAppFlash, motherPhys, false, 0);

  logGiunz2FinalAppFlash->SetVisAttributes(white);

  // protector2
  const G4double outRadiusprotector2Flash = 25 * mm;
  const G4double innRadiusprotector2Flash = 20.0 * mm;
  const G4double hightprotector2Flash = 30 * mm;
  const G4double startAngleprotector2Flash = 0. * deg;
  const G4double spanningAngleprotector2Flash = 360. * deg;
  const G4double protector2XPositionFlash = -957.44 * mm;

  protector2 =
      new G4Tubs("protector2", innRadiusprotector2Flash,
                 outRadiusprotector2Flash, hightprotector2Flash,
                 startAngleprotector2Flash, spanningAngleprotector2Flash);

  G4LogicalVolume *logprotector2Flash = new G4LogicalVolume(
      protector2, protector2material, "protector2", 0, 0, 0);

  physiprotector2 = new G4PVPlacement(
      G4Transform3D(rm5, G4ThreeVector((protector2XPositionFlash), 0., 0.)),
      "protector2", logprotector2Flash, motherPhys, false, 0);

  logprotector2Flash->SetVisAttributes(gray);

  // --------------------------------- //
  // monitor 1 FINAL APPLICATOR Flash //
  // --------------------------------- //

  const G4double innRadiusGiunz1FinalAppFlash = 10.0 * mm;
  const G4double outRadiusGiunz1FinalAppFlash = 15.0 * mm;
  const G4double hightGiunz1FinalAppFlash = 6.0 * mm;
  const G4double startAngleGiunz1FinalAppFlash = 0. * deg;
  const G4double spanningAngleGiunz1FinalAppFlash = 360. * deg;
  const G4double Giunz1FinalAppXPositionFlash = -993.44 * mm;

  solidGiunz1FinalAppFlash = new G4Tubs(
      "Giunz1FinalAppFlash", innRadiusGiunz1FinalAppFlash,
      outRadiusGiunz1FinalAppFlash, hightGiunz1FinalAppFlash,
      startAngleGiunz1FinalAppFlash, spanningAngleGiunz1FinalAppFlash);

  G4LogicalVolume *logGiunz1FinalAppFlash =
      new G4LogicalVolume(solidGiunz1FinalAppFlash, Giunz1FinalAppMaterialFlash,
                          "Giunz1FinalAppFlash", 0, 0, 0);

  physiGiunz1FinalAppFlash = new G4PVPlacement(
      G4Transform3D(rm5, G4ThreeVector((Giunz1FinalAppXPositionFlash), 0., 0.)),
      "Giunz1FinalAppFlash", logGiunz1FinalAppFlash, motherPhys, false, 0);

  logGiunz1FinalAppFlash->SetVisAttributes(gray);

  // protector1
  const G4double outRadiusprotector1Flash = 20 * mm;
  const G4double innRadiusprotector1Flash = 15 * mm;
  const G4double hightprotector1Flash = 6.0 * mm;
  const G4double startAngleprotector1Flash = 0. * deg;
  const G4double spanningAngleprotector1Flash = 360. * deg;
  const G4double protector1XPositionFlash = -993.44 * mm;

  protector1 =
      new G4Tubs("protector1", innRadiusprotector1Flash,
                 outRadiusprotector1Flash, hightprotector1Flash,
                 startAngleprotector1Flash, spanningAngleprotector1Flash);

  G4LogicalVolume *logprotector1Flash = new G4LogicalVolume(
      protector1, protector1material, "protector1", 0, 0, 0);

  physiprotector1 = new G4PVPlacement(
      G4Transform3D(rm5, G4ThreeVector((protector1XPositionFlash), 0., 0.)),
      "protector1", logprotector1Flash, motherPhys, false, 0);

  logprotector1Flash->SetVisAttributes(gray);

  // cover disk between monitor piece 1 and 2

  const G4double outRadiuscover1Flash = 20. * mm;
  const G4double innRadiuscover1Flash = 15 * mm;
  const G4double hightcover1Flash = 0.01 * mm;
  const G4double startAnglecover1Flash = 0. * deg;
  const G4double spanningAnglecover1Flash = 360. * deg;
  const G4double cover1XPositionFlash = -987.44 * mm;

  cover1 = new G4Tubs("cover1", innRadiuscover1Flash, outRadiuscover1Flash,
                      hightcover1Flash, startAnglecover1Flash,
                      spanningAnglecover1Flash);

  G4LogicalVolume *logcover1 =
      new G4LogicalVolume(cover1, cover1material, "cover1", 0, 0, 0);

  physicover1Flash = new G4PVPlacement(
      G4Transform3D(rm5, G4ThreeVector((cover1XPositionFlash), 0., 0.)),
      "cover1", logcover1, motherPhys, false, 0);

  logcover1->SetVisAttributes(gray);

  // cover disk between monitor piece 4 and first applicator

  const G4double outRadiuscover3Flash = 50.0 * mm;
  const G4double innRadiuscover3Flash = 24.75 * mm;
  const G4double hightcover3Flash = 0.01 * mm;
  const G4double startAnglecover3Flash = 0. * deg;
  const G4double spanningAnglecover3Flash = 360. * deg;
  const G4double cover3XPositionFlash = -899.44 * mm;

  cover3 = new G4Tubs("cover3", innRadiuscover3Flash, outRadiuscover3Flash,
                      hightcover3Flash, startAnglecover3Flash,
                      spanningAnglecover3Flash);

  G4LogicalVolume *logcover3 =
      new G4LogicalVolume(cover3, cover3material, "cover3", 0, 0, 0);

  physicover3Flash = new G4PVPlacement(
      G4Transform3D(rm5, G4ThreeVector((cover3XPositionFlash), 0., 0.)),
      "cover3", logcover3, motherPhys, false, 0);

  logcover3->SetVisAttributes(gray);
}

void Applicator80BeamLine::FlashBeamLineFinalApplicator() {
  // -----------------------//
  // FINAL APPLICATOR Flash  //
  //------------------------//

  // const G4double outRadiusFinalApplicatorFlash = 45. *mm;
  // const G4double innRadiusFinalApplicatorFlash = 40. *mm;
  const G4double hightFinalApplicatorFlash = 250.0 * mm;
  const G4double startAngleFinalApplicatorFlash = 0. * deg;
  const G4double spanningAngleFinalApplicatorFlash = 360. * deg;
  const G4double finalApplicatorXPositionFlash = -449.44 * mm;

  G4double phi6 = 90. * deg;

  G4RotationMatrix rm6;
  rm6.rotateY(phi6);

  solidFinalApplicatorFlash = new G4Tubs(
      "FinalApplicatorFlash", innerRadiusFinalApplicatorFlash,
      OuterRadiusFinalApplicatorFlash, hightFinalApplicatorFlash,
      startAngleFinalApplicatorFlash, spanningAngleFinalApplicatorFlash);

  G4LogicalVolume *logFinalApplicatorFlash = new G4LogicalVolume(
      solidFinalApplicatorFlash, finalApplicatorMaterialFlash,
      "FinalApplicatorFlash", 0, 0, 0);

  physiFinalApplicatorFlash = new G4PVPlacement(
      G4Transform3D(rm6,
                    G4ThreeVector((finalApplicatorXPositionFlash), 0., 0.)),
      "FinalApplicatorFlash", logFinalApplicatorFlash, motherPhys, false, 0);

  //  logFinalApplicatorFlash ->
  //  SetVisAttributes(G4VisAttributes::GetInvisible());
  logFinalApplicatorFlash->SetVisAttributes(magenta);
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// MESSENGER ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void Applicator80BeamLine::SetInnerRadiusFinalApplicatorFlash(
    G4double value) { /// NB here it isn't enough to set the radius of
                      /// applicator, but also of all the junctions
  solidFinalApplicatorFlash->SetInnerRadius(value);
  solidFirstApplicatorFlash->SetInnerRadius(value);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "Inner Radius of the applicator Flash is (mm):"
         << solidFinalApplicatorFlash->GetInnerRadius() / mm << G4endl;
}

/////////////////////////////////////////////////////////////////////////

void Applicator80BeamLine::SetOuterRadiusFinalApplicatorFlash(G4double value) {
  solidFinalApplicatorFlash->SetOuterRadius(value);
  solidFirstApplicatorFlash->SetOuterRadius(value);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "Outer Radius of the final applicator Flash is (mm):"
         << solidFinalApplicatorFlash->GetOuterRadius() / mm << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
