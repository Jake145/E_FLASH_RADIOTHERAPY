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
/// \file FlashDetectorConstruction.cc
/// \brief Implementation of the FlashDetectorConstruction class

#include "FlashDetectorConstruction.hh"
//#include "FlashHit.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Region.hh"
#include "G4SDManager.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"

#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"

#include "G4UserLimits.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SystemOfUnits.hh"

#include "Applicator80BeamLine.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4SystemOfUnits.hh"
//#include "FSensitiveDetector.hh"
//#include "OpticFiberSD.hh"
//#include "PhotoDiodeSD.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PhysicalConstants.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VisAttributes.hh"
//#include "G4VSensitiveDetector.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashDetectorConstruction::FlashDetectorConstruction()
    : G4VUserDetectorConstruction(), physicalTreatmentRoom(0), Collimator(0),
      fCheckOverlaps(true) {

  DefineMaterials();

  // Construct();
  // ConstructSDandField();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashDetectorConstruction::~FlashDetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FlashDetectorConstruction::DefineMaterials() {
  G4NistManager *man = G4NistManager::Instance();
  G4int ncomponents;
  G4bool isotopes = false;
  G4double prelude_density = 7.4 * g / cm3;
  G4Material *prelude =
      new G4Material("prelude", prelude_density, ncomponents = 4);
  prelude->AddElement(man->FindOrBuildElement("Lu"), 71 * perCent);
  prelude->AddElement(man->FindOrBuildElement("Si"), 7 * perCent);
  prelude->AddElement(man->FindOrBuildElement("O"), 18 * perCent);
  prelude->AddElement(man->FindOrBuildElement("Y"), 4 * perCent);

  G4Material *scintillator =
      new G4Material("scintillator", prelude_density, ncomponents = 2);
  scintillator->AddMaterial(prelude, 99.81 * perCent);
  scintillator->AddElement(man->FindOrBuildElement("Ce"), 0.19 * perCent);

  G4MaterialPropertiesTable *mpt = new G4MaterialPropertiesTable();
  const G4int num = 20;
  G4double ene[num] = {1.79 * eV, 1.85 * eV, 1.91 * eV, 1.97 * eV,

                       2.04 * eV, 2.11 * eV, 2.19 * eV, 2.27 * eV,

                       2.36 * eV, 2.45 * eV, 2.56 * eV, 2.67 * eV,

                       2.80 * eV, 2.94 * eV, 3.09 * eV, 3.25 * eV,

                       3.44 * eV, 3.65 * eV, 3.89 * eV, 4.16 * eV};

  G4double fast[num] = {0.01, 0.10, 0.20, 0.50,

                        0.90, 1.70, 2.90, 5.00,

                        8.30, 12.5, 17.0, 22.9,

                        26.4, 25.6, 16.8, 4.20,

                        0.30, 0.20, 0.10, 0.01};

  G4double rLyso[num] = {1.81, 1.81, 1.81, 1.81,

                         1.81, 1.81, 1.81, 1.81,

                         1.81, 1.81, 1.81, 1.81,

                         1.81, 1.81, 1.81, 1.81,

                         1.81, 1.81, 1.81, 1.81};

  G4double abs[num] = {3.5 * m, 3.5 * m, 3.5 * m, 3.5 * m,

                       3.5 * m, 3.5 * m, 3.5 * m, 3.5 * m,

                       3.5 * m, 3.5 * m, 3.5 * m, 3.5 * m,

                       3.5 * m, 3.5 * m, 3.5 * m, 3.5 * m,

                       3.5 * m, 3.5 * m, 3.5 * m, 3.5 * m};

  mpt->AddProperty("FASTCOMPONENT", ene, fast, num);

  mpt->AddProperty("RINDEX", ene, rLyso, num);

  mpt->AddProperty("ABSLENGTH", ene, abs, num);

  mpt->AddConstProperty("SCINTILLATIONYIELD", 2700 / MeV);

  mpt->AddConstProperty("RESOLUTIONSCALE", 1);

  mpt->AddConstProperty("FASTTIMECONSTANT", 41 * ns);
  scintillator->SetMaterialPropertiesTable(mpt);

  // fluorinated polymer
  G4int polyeth = 1;
  G4int nC_eth = 2 * polyeth;
  G4int nH_eth = 4 * polyeth;
  G4double z, a, density;
  fH = new G4Element("H", "H", z = 1., a = 1.01 * g / mole);
  fC = new G4Element("C", "C", z = 6., a = 12.01 * g / mole);
  fN = new G4Element("N", "N", z = 7., a = 14.01 * g / mole);
  fO = new G4Element("O", "O", z = 8., a = 16.00 * g / mole);
  fPethylene2 = new G4Material("Pethylene2", 1400 * kg / m3, 2);
  fPethylene2->AddElement(fH, nH_eth);
  fPethylene2->AddElement(fC, nC_eth);

  //--------------------------------------------------
  // Fluorinated Polyethylene
  //--------------------------------------------------
  G4double PhotonEnergy[] = {
      2.00 * eV, 2.03 * eV, 2.06 * eV, 2.09 * eV, 2.12 * eV, 2.15 * eV,
      2.18 * eV, 2.21 * eV, 2.24 * eV, 2.27 * eV, 2.30 * eV, 2.33 * eV,
      2.36 * eV, 2.39 * eV, 2.42 * eV, 2.45 * eV, 2.48 * eV, 2.51 * eV,
      2.54 * eV, 2.57 * eV, 2.60 * eV, 2.63 * eV, 2.66 * eV, 2.69 * eV,
      2.72 * eV, 2.75 * eV, 2.78 * eV, 2.81 * eV, 2.84 * eV, 2.87 * eV,
      2.90 * eV, 2.93 * eV, 2.96 * eV, 2.99 * eV, 3.02 * eV, 3.05 * eV,
      3.08 * eV, 3.11 * eV, 3.14 * eV, 3.17 * eV, 3.20 * eV, 3.23 * eV,
      3.26 * eV, 3.29 * eV, 3.32 * eV, 3.35 * eV, 3.38 * eV, 3.41 * eV,
      3.44 * eV, 3.47 * eV};

  const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

  G4double refractiveIndexClad2[] = {
      1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
      1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
      1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
      1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
      1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42};

  assert(sizeof(refractiveIndexClad2) == sizeof(PhotonEnergy));
  G4double absClad[] = {
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m};

  assert(sizeof(absClad) == sizeof(PhotonEnergy));

  // Add entries into properties table
  G4MaterialPropertiesTable *mptClad2 = new G4MaterialPropertiesTable();
  mptClad2->AddProperty("RINDEX", PhotonEnergy, refractiveIndexClad2, nEntries);
  mptClad2->AddProperty("ABSLENGTH", PhotonEnergy, absClad, nEntries);

  fPethylene2->SetMaterialPropertiesTable(mptClad2);



  G4double fractionmass;
  std::vector<G4int> natoms;
  std::vector<G4double> fractionMass;
  std::vector<G4String> elements;
  
  
  nist = G4NistManager::Instance();
  cryst_mat = nist->FindOrBuildMaterial("scintillator");
  //PMMA = nist->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);

  elements.push_back("C");     natoms.push_back(5);
  elements.push_back("H");     natoms.push_back(8);
  elements.push_back("O");     natoms.push_back(2);

  density = 1.190*g/cm3;

  PMMA = nist->
          ConstructNewMaterial("PMMA", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Cladding (polyethylene)
  //--------------------------------------------------

  elements.push_back("C");     natoms.push_back(2);
  elements.push_back("H");     natoms.push_back(4);

  density = 1.200*g/cm3;

  PE = nist->
          ConstructNewMaterial("Pethylene", elements, natoms, density);

  elements.clear();
  natoms.clear();
  //--------------------------------------------------
  //  PMMA for fibers
  //--------------------------------------------------
  G4double photonEnergy_fib[] = {
      2.00 * eV, 2.03 * eV, 2.06 * eV, 2.09 * eV, 2.12 * eV, 2.15 * eV,
      2.18 * eV, 2.21 * eV, 2.24 * eV, 2.27 * eV, 2.30 * eV, 2.33 * eV,
      2.36 * eV, 2.39 * eV, 2.42 * eV, 2.45 * eV, 2.48 * eV, 2.51 * eV,
      2.54 * eV, 2.57 * eV, 2.60 * eV, 2.63 * eV, 2.66 * eV, 2.69 * eV,
      2.72 * eV, 2.75 * eV, 2.78 * eV, 2.81 * eV, 2.84 * eV, 2.87 * eV,
      2.90 * eV, 2.93 * eV, 2.96 * eV, 2.99 * eV, 3.02 * eV, 3.05 * eV,
      3.08 * eV, 3.11 * eV, 3.14 * eV, 3.17 * eV, 3.20 * eV, 3.23 * eV,
      3.26 * eV, 3.29 * eV, 3.32 * eV, 3.35 * eV, 3.38 * eV, 3.41 * eV,
      3.44 * eV, 3.47 * eV};

  const G4int nEntries_fib = sizeof(photonEnergy_fib) / sizeof(G4double);

  G4double refractiveIndex_fiber[] = {
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49};

  assert(sizeof(refractiveIndex_fiber) == sizeof(photonEnergy_fib));

  G4double abs_fiber[] = {
      5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m,
      5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m,
      5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m,
      5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m,
      5.40 * m, 1.10 * m, 1.10 * m, 1.10 * m, 1.10 * m, 1.10 * m, 1.10 * m,
      1.10 * m, 1. * mm,  1. * mm,  1. * mm,  1. * mm,  1. * mm,  1. * mm,
      1. * mm,  1. * mm,  1. * mm,  1. * mm,  1. * mm,  1. * mm,  1. * mm,
      1. * mm};

  assert(sizeof(abs_fiber) == sizeof(photonEnergy_fib));

  /* G4double emission_fib[] =
   {0.05, 0.10, 0.30, 0.50, 0.75, 1.00, 1.50, 1.85, 2.30, 2.75,
    3.25, 3.80, 4.50, 5.20, 6.00, 7.00, 8.50, 9.50, 11.1, 12.4,
    12.9, 13.0, 12.8, 12.3, 11.1, 11.0, 12.0, 11.0, 17.0, 16.9,
    15.0, 9.00, 2.50, 1.00, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

   assert(sizeof(emission_fib) == sizeof(photonEnergy_fib)); */

  // Add entries into properties table
  G4MaterialPropertiesTable *mpt_fiber = new G4MaterialPropertiesTable();
  mpt_fiber->AddProperty("RINDEX", photonEnergy_fib, refractiveIndex_fiber,
                         nEntries_fib);

  mpt_fiber->AddProperty("ABSLENGTH", photonEnergy_fib, abs_fiber,
                         nEntries_fib);
  // mptWLSfiber->AddProperty("WLSCOMPONENT",photonEnergy_fib,emissionFib,nEntries_fib);

  PMMA->SetMaterialPropertiesTable(mpt_fiber);
  G4cout << "PMMA G4MaterialPropertiesTable" << G4endl;
  mpt_fiber->DumpTable();

  //PE = nist->FindOrBuildMaterial("G4_POLYETHYLENE", isotopes);

  G4double refractiveIndex_Clad[] = {
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49};

  assert(sizeof(refractiveIndex_Clad) == sizeof(photonEnergy_fib));

  G4double abs_Clad[] = {
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m};

  assert(sizeof(abs_Clad) == sizeof(photonEnergy_fib));

  // Add entries into properties table
  G4MaterialPropertiesTable *mptClad = new G4MaterialPropertiesTable();
  mptClad->AddProperty("RINDEX", photonEnergy_fib, refractiveIndex_Clad,
                       nEntries_fib);
  mptClad->AddProperty("ABSLENGTH", photonEnergy_fib, abs_Clad, nEntries_fib);

  PE->SetMaterialPropertiesTable(mptClad);

    TEFLON = nist->FindOrBuildMaterial("G4_TEFLON",isotopes);

  G4double photonEnergy_teflon[] =
              { 7.897*eV,7.208*eV, 6.702*eV,  4.999*eV};
    G4double refractiveIndex3[] =
              { 1.432, 1.308, 1.364, 1.329};
  const G4int nEntries_teflon = sizeof(photonEnergy_teflon)/sizeof(G4double);
    G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();
    myMPT3->AddProperty("RINDEX", photonEnergy_teflon, refractiveIndex3,
  nEntries_teflon);

    G4cout << "Teflon G4MaterialPropertiesTable" << G4endl;
    myMPT3->DumpTable();

    TEFLON->SetMaterialPropertiesTable(myMPT3);
    
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume * FlashDetectorConstruction::ConstructPhantom(G4double Cx,G4double Cy,G4double Cz,G4double d,G4double Oz) {

  G4Material *phantomMaterial = nist->FindOrBuildMaterial("G4_WATER");

  // ------------ Generate & Add Material Properties Table ------------
  //

  G4double photonEnergy[] = {
      2.034 * eV, 2.068 * eV, 2.103 * eV, 2.139 * eV, 2.177 * eV, 2.216 * eV,
      2.256 * eV, 2.298 * eV, 2.341 * eV, 2.386 * eV, 2.433 * eV, 2.481 * eV,
      2.532 * eV, 2.585 * eV, 2.640 * eV, 2.697 * eV, 2.757 * eV, 2.820 * eV,
      2.885 * eV, 2.954 * eV, 3.026 * eV, 3.102 * eV, 3.181 * eV, 3.265 * eV,
      3.353 * eV, 3.446 * eV, 3.545 * eV, 3.649 * eV, 3.760 * eV, 3.877 * eV,
      4.002 * eV, 4.136 * eV};

  const G4int nEntries = sizeof(photonEnergy) / sizeof(G4double);

  //
  // Water
  //
  G4double refractiveIndex1[] = {
      1.3435, 1.344,  1.3445, 1.345,  1.3455, 1.346,  1.3465, 1.347,
      1.3475, 1.348,  1.3485, 1.3492, 1.35,   1.3505, 1.351,  1.3518,
      1.3522, 1.3530, 1.3535, 1.354,  1.3545, 1.355,  1.3555, 1.356,
      1.3568, 1.3572, 1.358,  1.3585, 1.359,  1.3595, 1.36,   1.3608};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

 /* G4double absorption[] = {
      3.448 * m,  4.082 * m,  6.329 * m,  9.174 * m,  12.346 * m, 13.889 * m,
      15.152 * m, 17.241 * m, 18.868 * m, 20.000 * m, 26.316 * m, 35.714 * m,
      45.455 * m, 47.619 * m, 52.632 * m, 52.632 * m, 55.556 * m, 52.632 * m,
      52.632 * m, 47.619 * m, 45.455 * m, 41.667 * m, 37.037 * m, 33.333 * m,
      30.000 * m, 28.500 * m, 27.000 * m, 24.500 * m, 22.000 * m, 19.500 * m,
      17.500 * m, 14.500 * m};

  assert(sizeof(absorption) == sizeof(photonEnergy));

  G4double scintilFast[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  assert(sizeof(scintilFast) == sizeof(photonEnergy));

  G4double scintilSlow[] =
            { 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
              7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00,
              3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00,
              4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
              7.00, 6.00, 5.00, 4.00 };

  assert(sizeof(scintilSlow) == sizeof(photonEnergy));*/
  G4MaterialPropertiesTable *myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX", photonEnergy, refractiveIndex1, nEntries)
      ->SetSpline(true);
  /*myMPT1->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries)
      ->SetSpline(true);
  myMPT1->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
        ->SetSpline(true);
  myMPT1->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
        ->SetSpline(true);

  myMPT1->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  myMPT1->AddConstProperty("YIELDRATIO",0.8);

  G4double energy_water[] = {
      1.56962 * eV, 1.58974 * eV, 1.61039 * eV, 1.63157 * eV, 1.65333 * eV,
      1.67567 * eV, 1.69863 * eV, 1.72222 * eV, 1.74647 * eV, 1.77142 * eV,
      1.7971 * eV,  1.82352 * eV, 1.85074 * eV, 1.87878 * eV, 1.90769 * eV,
      1.93749 * eV, 1.96825 * eV, 1.99999 * eV, 2.03278 * eV, 2.06666 * eV,
      2.10169 * eV, 2.13793 * eV, 2.17543 * eV, 2.21428 * eV, 2.25454 * eV,
      2.29629 * eV, 2.33962 * eV, 2.38461 * eV, 2.43137 * eV, 2.47999 * eV,
      2.53061 * eV, 2.58333 * eV, 2.63829 * eV, 2.69565 * eV, 2.75555 * eV,
      2.81817 * eV, 2.88371 * eV, 2.95237 * eV, 3.02438 * eV, 3.09999 * eV,
      3.17948 * eV, 3.26315 * eV, 3.35134 * eV, 3.44444 * eV, 3.54285 * eV,
      3.64705 * eV, 3.75757 * eV, 3.87499 * eV, 3.99999 * eV, 4.13332 * eV,
      4.27585 * eV, 4.42856 * eV, 4.59258 * eV, 4.76922 * eV, 4.95999 * eV,
      5.16665 * eV, 5.39129 * eV, 5.63635 * eV, 5.90475 * eV, 6.19998 * eV};

  const G4int numentries_water = sizeof(energy_water) / sizeof(G4double);

  // assume 100 times larger than the rayleigh scattering for now.
  G4double mie_water[] = {
      167024.4 * m, 158726.7 * m, 150742 * m,   143062.5 * m, 135680.2 * m,
      128587.4 * m, 121776.3 * m, 115239.5 * m, 108969.5 * m, 102958.8 * m,
      97200.35 * m, 91686.86 * m, 86411.33 * m, 81366.79 * m, 76546.42 * m,
      71943.46 * m, 67551.29 * m, 63363.36 * m, 59373.25 * m, 55574.61 * m,
      51961.24 * m, 48527.00 * m, 45265.87 * m, 42171.94 * m, 39239.39 * m,
      36462.50 * m, 33835.68 * m, 31353.41 * m, 29010.30 * m, 26801.03 * m,
      24720.42 * m, 22763.36 * m, 20924.88 * m, 19200.07 * m, 17584.16 * m,
      16072.45 * m, 14660.38 * m, 13343.46 * m, 12117.33 * m, 10977.70 * m,
      9920.416 * m, 8941.407 * m, 8036.711 * m, 7202.470 * m, 6434.927 * m,
      5730.429 * m, 5085.425 * m, 4496.467 * m, 3960.210 * m, 3473.413 * m,
      3032.937 * m, 2635.746 * m, 2278.907 * m, 1959.588 * m, 1675.064 * m,
      1422.710 * m, 1200.004 * m, 1004.528 * m, 833.9666 * m, 686.1063 * m};

  assert(sizeof(mie_water) == sizeof(energy_water));
  // gforward, gbackward, forward backward ratio
  G4double mie_water_const[3] = {0.99, 0.99, 0.8};

  myMPT1->AddProperty("MIEHG", energy_water, mie_water, numentries_water)
      ->SetSpline(true);
  myMPT1->AddConstProperty("MIEHG_FORWARD", mie_water_const[0]);
  myMPT1->AddConstProperty("MIEHG_BACKWARD", mie_water_const[1]);
  myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO", mie_water_const[2]);

  G4cout << "Water G4MaterialPropertiesTable" << G4endl;
  myMPT1->DumpTable();
phantomMaterial->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);*/
  phantomMaterial->SetMaterialPropertiesTable(myMPT1);
  

  G4double phantomSizeX = 10.0 * cm, phantomSizeY = 30.0 * cm,
           phantomSizeZ = 30.0 * cm;
  G4ThreeVector phantomPosition = G4ThreeVector(-(199.4 * mm - phantomSizeX/2) , 0. * mm, 0. * mm);
  // Definition of the solid volume of the Phantom
  phantom = new G4Box("Phantom", phantomSizeX / 2, phantomSizeY / 2,
                      phantomSizeZ / 2);

  // Definition of the logical volume of the Phantom
  phantomLogicalVolume =
      new G4LogicalVolume(phantom, phantomMaterial, "phantomLog", 0, 0, 0);

  // Definition of the physics volume of the Phantom
  phant_phys =
      new G4PVPlacement(0, phantomPosition, "phantomPhys", phantomLogicalVolume,
                        physicalTreatmentRoom, false, 0);
  /*
  G4OpticalSurface* opWaterSurface = new G4OpticalSurface("WaterSurface");
    opWaterSurface->SetType(dielectric_LUTDAVIS);
    opWaterSurface->SetFinish(Rough_LUT);
    opWaterSurface->SetModel(DAVIS);

    G4LogicalBorderSurface* waterSurface =
            new G4LogicalBorderSurface("WaterSurface",
                                   phant_phys,physicalTreatmentRoom,opWaterSurface);

    G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
          (waterSurface->GetSurface(phant_phys,physicalTreatmentRoom)->
                                                         GetSurfaceProperty());
    if (opticalSurface) opticalSurface->DumpInfo();
    */
    
    
  //============================DETECTOR_SUPPORT=====================================//
  
  //G4double support_x=1*cm, wedge_X=1*mm,wedge_Y=1*mm,wedge_Z=1*cm;
  G4double support_x=1*cm, wedge_X=Cz+2*d,wedge_Y=Cy+2*d,wedge_Z=Cx+Oz+2*d; 
  G4ThreeVector xTrans(-(support_x/2-wedge_X/2), 0, -wedge_Z/2);
  G4Box* support_whole= new G4Box("Support_w",support_x/2,10*cm,10*cm);
  G4Box* wedge = new G4Box("Wedge",wedge_X/2,wedge_Y/2,wedge_Z/2);
  G4SubtractionSolid* support_wedged = new G4SubtractionSolid("Wedged_Support",support_whole,wedge,0,xTrans);
  
  DetectorSupport =
      new G4LogicalVolume(support_wedged, PMMA, "SupportLog");

G4RotationMatrix rotmp = G4RotationMatrix();
  rotmp.rotateY(0 * deg);

  G4ThreeVector positionp = G4ThreeVector((phantomSizeX/2+support_x/2),0,Cx/2);
  G4Transform3D transformp = G4Transform3D(rotmp, positionp);
  
  G4VPhysicalVolume *DTsupp = new G4PVPlacement(
      transformp, DetectorSupport, "supportphys", phantomLogicalVolume, false, 0);
  
  //=================================================================================//
  //================Second Piece of Phantom==================================//
  
  G4ThreeVector phantomPosition_2 = G4ThreeVector((15 * cm -(phantomSizeX / 2 +support_x/2))+support_x/2 , 0. * mm, -Cx/2);
  // Definition of the solid volume of the Phantom
  G4Box * phantom_2 = new G4Box("Phantom_2", 15 * cm -(phantomSizeX / 2 +support_x/2) , phantomSizeY / 2,
                      phantomSizeZ / 2);

  // Definition of the logical volume of the Phantom
  G4LogicalVolume* phantomLogicalVolume_2 =
      new G4LogicalVolume(phantom_2, phantomMaterial, "phantomLog_2", 0, 0, 0);

  // Definition of the physics volume of the Phantom
  G4VPhysicalVolume * phant_phys_2 =
      new G4PVPlacement(0, phantomPosition_2, "phantomPhys_2", phantomLogicalVolume_2,
                        DTsupp, false, 0);
  
  
  
  //==========================================================================//
  
  G4Region *PhantomRegion = new G4Region("Phantom_reg");
  phantomLogicalVolume->SetRegion(PhantomRegion);
  PhantomRegion->AddRootLogicalVolume(phantomLogicalVolume);

  // Visualisation attributes of the phantom
  red = new G4VisAttributes(G4Colour(0 / 255., 255 / 255., 0 / 255.));
  red->SetVisibility(true);

blue = new G4VisAttributes(G4Colour(0 / 255., 0./ 255., 255. / 255.));
  blue->SetVisibility(true);

  phantomLogicalVolume->SetVisAttributes(red);
 DetectorSupport->SetVisAttributes(blue);
   phantomLogicalVolume_2->SetVisAttributes(red);
  G4double maxStep = 0.1 * mm;
  fStepLimit = new G4UserLimits(maxStep);
  phantomLogicalVolume->SetUserLimits(fStepLimit);

  return phant_phys;
}

G4VPhysicalVolume *FlashDetectorConstruction::Construct() {
  // -----------------------------
  // Treatment room - World volume
  //------------------------------
  // Treatment room sizes
  const G4double worldX = 400.0 * cm;
  const G4double worldY = 400.0 * cm;
  const G4double worldZ = 400.0 * cm;
  G4bool isotopes = false;

  G4Material *airNist =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  // Air
  //
  /*
  G4double photonEnergy_air[] =
              { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
                2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
                2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
                2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
                2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
                3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
                3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
                3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
    G4double refractiveIndex2[] =
              { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00 };
  const G4int nEntries_air = sizeof(photonEnergy_air)/sizeof(G4double);
    G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
    myMPT2->AddProperty("RINDEX", photonEnergy_air, refractiveIndex2,
  nEntries_air);

    G4cout << "Air G4MaterialPropertiesTable" << G4endl;
    myMPT2->DumpTable();

    airNist->SetMaterialPropertiesTable(myMPT2); */

  G4Box *treatmentRoom = new G4Box("TreatmentRoom", worldX, worldY, worldZ);
  G4LogicalVolume *logicTreatmentRoom = new G4LogicalVolume(
      treatmentRoom, airNist, "logicTreatmentRoom", 0, 0, 0);
  physicalTreatmentRoom =
      new G4PVPlacement(0, G4ThreeVector(), "physicalTreatmentRoom",
                        logicTreatmentRoom, 0, false, 0);

  // The treatment room is invisible in the Visualisation
  logicTreatmentRoom->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  // -----------------------------
  // detector + phantom +Default dimensions
  //------------------------------

  red = new G4VisAttributes(G4Colour(0 / 255., 255 / 255., 0 / 255.));
  red->SetVisibility(true);
  green = new G4VisAttributes(G4Colour(255 / 255., 0. / 255., 0 / 255.));
  green->SetVisibility(true);
//crystal dims
  G4double cryst_dX = 1 * cm, cryst_dY = 1 * mm, cryst_dZ = 1 * mm;
  G4double gap = 0 * mm; // a gap for wrapping, change to add gap
  G4double dX = cryst_dX - gap, dY = cryst_dY - gap, dZ = cryst_dZ - gap;
 //fiber dim
 G4double opticfiber_core_dx = 5 * cm;
  G4double opticfiber_core_diameter = 0.98 * mm;
  G4double optic_fiber_clad_diameter = 2.2 * mm;
  G4double optic_fiber_cladding_diameter = 1. * mm;
  //Teflon dim
  G4double fPTFEThickness = 0.3 * mm;
  //construct collimatore
    Collimator = new Applicator80BeamLine(physicalTreatmentRoom);
  // constuct phantom
  phantom_physical=ConstructPhantom(dX,dY,dZ,fPTFEThickness,opticfiber_core_dx);
  //constuct detector
  G4Box *solidCryst = new G4Box("crystal", dX / 2, dY / 2, dZ / 2);
  //
  logicCryst = new G4LogicalVolume(solidCryst,   // its solid
                                   cryst_mat,    // its material
                                   "CrystalLV"); // its name
  G4RotationMatrix rotm = G4RotationMatrix();
  rotm.rotateY(90 * deg);

  G4ThreeVector position = G4ThreeVector(-(1*cm/2-dZ/2-fPTFEThickness), 0, -(dX/2+fPTFEThickness));
  G4Transform3D transform = G4Transform3D(rotm, position);

  G4VPhysicalVolume *phys_cryst = new G4PVPlacement(
      transform, logicCryst, "crystalphys", DetectorSupport, false, 0);
     
//OOOOOOOOOOOOOOOOOOoooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOOoo
  G4OpticalSurface *opteflonSurface_up =
      new G4OpticalSurface("teflonSurface");
  opteflonSurface_up->SetType(dielectric_LUTDAVIS);

  opteflonSurface_up->SetModel(DAVIS);
  opteflonSurface_up->SetFinish(PolishedTeflon_LUT);


 G4LogicalBorderSurface *teflonSurface_up =
      new G4LogicalBorderSurface("teflonSurface_up", phys_cryst,
                                 phantom_physical, opteflonSurface_up);
  G4LogicalBorderSurface *teflonSurface_down =
      new G4LogicalBorderSurface("teflonSurface_down", phantom_physical,
                                 phys_cryst, opteflonSurface_up);

  G4OpticalSurface *opticalSurface_teflon = dynamic_cast<G4OpticalSurface *>(
      teflonSurface_up->GetSurface(phys_cryst, phantom_physical)
          ->GetSurfaceProperty());
  if (opticalSurface_teflon)
    opticalSurface_teflon->DumpInfo();
//OOOOOOOOOOOOOOOOOOoooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOOoo
  // optic fiber
  //
  
  ////////////////////

  G4Tubs *opticfiber_cladding = new G4Tubs(
      "OF_clad", opticfiber_core_diameter / 2,
      optic_fiber_cladding_diameter / 2, opticfiber_core_dx / 2, 0., twopi);

  ////////////////////

  G4Tubs *opticfiber_core =
      new G4Tubs("OF_core", 0. * cm, opticfiber_core_diameter / 2,
                 opticfiber_core_dx / 2, 0., CLHEP::twopi);

  opticfiber_core_log = new G4LogicalVolume(opticfiber_core, // its solid
                                            PMMA,            // its material
                                            "OF_core_LV");   // its name
  G4Tubs *opticfiber_clad = new G4Tubs(
      "OF_clad", opticfiber_core_diameter / 2,
      (optic_fiber_clad_diameter - optic_fiber_cladding_diameter) / 2,
      opticfiber_core_dx / 2, 0., twopi);

  G4LogicalVolume *opticfiber_cladding_log =
      new G4LogicalVolume(opticfiber_cladding, // its solid
                          fPethylene2,         // its material
                          "OF_cladding_LV");   // its name

  G4LogicalVolume *opticfiber_clad_log =
      new G4LogicalVolume(opticfiber_clad, // its solid
                          PE,              // its material
                          "OF_clad_LV");   // its name
  G4ThreeVector position_opt = G4ThreeVector(0. * mm, 0. * mm, 0. * mm);
  G4ThreeVector position_clad = G4ThreeVector(
      0. * mm + (dX / 2 + opticfiber_core_dx / 2), 0. * mm, 0. * mm);

  G4RotationMatrix rotm_opt = G4RotationMatrix();
  rotm_opt.rotateY(0 * deg);
  G4Transform3D transform_opt = G4Transform3D(rotm_opt, position_opt);

  G4RotationMatrix rotm_clad = G4RotationMatrix();
  rotm_clad.rotateY(90 * deg);
  G4Transform3D transform_clad = G4Transform3D(rotm_clad, position_clad);

  //
  G4VPhysicalVolume *physclad = new G4PVPlacement(
      transform_clad,
      // position
      opticfiber_clad_log, "outerfiber", logicCryst, false, 0);

  G4VPhysicalVolume *physcore = new G4PVPlacement(transform_clad,
                                                  // position
                                                  opticfiber_core_log,

                                                  "OF_core_phys",

                                                  logicCryst, false, 0);
  G4VPhysicalVolume *claddingcore = new G4PVPlacement(transform_clad,
                                                      // position
                                                      opticfiber_cladding_log,

                                                      "OF_cladding_phys",

                                                      logicCryst, false, 0);

  G4OpticalSurface *opcore_scint =
      new G4OpticalSurface("OpticFiberandScintillator");
  opcore_scint->SetType(dielectric_LUTDAVIS);
  opcore_scint->SetFinish(PolishedESRGrease_LUT);
  opcore_scint->SetModel(DAVIS);
  G4LogicalBorderSurface *core_scint_up = new G4LogicalBorderSurface(
      "teflonSurface_of_up", phys_cryst, physcore, opcore_scint);
  G4LogicalBorderSurface *core_scint_down = new G4LogicalBorderSurface(
      "teflonSurface_of_down", physcore, phys_cryst, opcore_scint);
  // OOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooo

  G4VSolid * t1= new G4Box("t1",dX/2+fPTFEThickness, dY/2+fPTFEThickness, dZ/2+fPTFEThickness);
  G4VSolid *t2 = new G4Box("t2",dX, dY, dZ);

  G4ThreeVector zTrans(0, 0,0);

  G4SubtractionSolid* solidFullTeflon = new G4SubtractionSolid("hollowteflon", t1, t2, 0, zTrans);


  G4VSolid *tHole = new G4Tubs("sTeflonHole",
			       0.0*cm,
			       (optic_fiber_clad_diameter - optic_fiber_cladding_diameter) / 2,
			       fPTFEThickness/2,
			       0.*deg,
			       360.*deg);
  G4RotationMatrix*yRot =new G4RotationMatrix;


  //yRot->rotateX(90*deg);
  yRot->rotateY(90*deg);
  G4ThreeVector zHole((dX/2 + fPTFEThickness/2), 0,0  );

  G4SubtractionSolid* solidHoleTeflon = new G4SubtractionSolid("hollowHoleteflon", solidFullTeflon, tHole, yRot, zHole);
 G4LogicalVolume *logicTeflon = new G4LogicalVolume(solidHoleTeflon, TEFLON,"lTeflon",0,0,0);
 G4VPhysicalVolume *physiTeflon = new G4PVPlacement(0, G4ThreeVector(0,0,0), logicTeflon, "pTeflon", logicCryst, false, 0);
  
  
  //OOOOOOOOOOOOOOOOOOoooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOOoo
  G4OpticalSurface *scintteflonSurface_up =
      new G4OpticalSurface("steflonSurface");
  scintteflonSurface_up->SetType(dielectric_LUTDAVIS);

  scintteflonSurface_up->SetModel(DAVIS);
  scintteflonSurface_up->SetFinish(PolishedTeflon_LUT);


 G4LogicalBorderSurface *ssteflonSurface_up =
      new G4LogicalBorderSurface("steflonSurface_up", phys_cryst,
                                 physiTeflon, scintteflonSurface_up);
  G4LogicalBorderSurface *ssteflonSurface_down =
      new G4LogicalBorderSurface("steflonSurface_down", physiTeflon,
                                 phys_cryst, scintteflonSurface_up);

  G4OpticalSurface *scintteflonSurface_up_0 = dynamic_cast<G4OpticalSurface *>(
      ssteflonSurface_up->GetSurface(phys_cryst, physiTeflon)
          ->GetSurfaceProperty());
  if (scintteflonSurface_up_0)
    scintteflonSurface_up_0->DumpInfo();
//OOOOOOOOOOOOOOOOOOoooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOOoo
  //OOOOOOOOOOOOOOOOOOoooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOOoo
  G4OpticalSurface *phantteflonSurface_up =
      new G4OpticalSurface("pteflonSurface");
  phantteflonSurface_up->SetType(dielectric_LUTDAVIS);

  phantteflonSurface_up->SetModel(DAVIS);
  phantteflonSurface_up->SetFinish(PolishedTeflon_LUT);


 G4LogicalBorderSurface *pteflonSurface_up =
      new G4LogicalBorderSurface("pteflonSurface_up", phant_phys,
                                 physiTeflon, phantteflonSurface_up);
  G4LogicalBorderSurface *pteflonSurface_down =
      new G4LogicalBorderSurface("pteflonSurface_down", physiTeflon,
                                 phant_phys, phantteflonSurface_up);

  G4OpticalSurface *pteflonSurface_up_0 = dynamic_cast<G4OpticalSurface *>(
      pteflonSurface_up->GetSurface(phant_phys, physiTeflon)
          ->GetSurfaceProperty());
  if (pteflonSurface_up_0)
    pteflonSurface_up_0->DumpInfo();
//OOOOOOOOOOOOOOOOOOoooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOOoo

G4OpticalSurface *cladteflonSurface_up =
      new G4OpticalSurface("cteflonSurface");
  cladteflonSurface_up->SetType(dielectric_LUTDAVIS);

  cladteflonSurface_up->SetModel(DAVIS);
  cladteflonSurface_up->SetFinish(PolishedTeflon_LUT);


 G4LogicalBorderSurface *cteflonSurface_up =
      new G4LogicalBorderSurface("cteflonSurface_up", physclad,
                                 physiTeflon, cladteflonSurface_up);
  G4LogicalBorderSurface *cteflonSurface_down =
      new G4LogicalBorderSurface("cteflonSurface_down", physiTeflon,
                                 physclad, cladteflonSurface_up);

  G4OpticalSurface *cteflonSurface_up_0 = dynamic_cast<G4OpticalSurface *>(
      cteflonSurface_up->GetSurface(physclad, physiTeflon)
          ->GetSurfaceProperty());
  if (cteflonSurface_up_0)
    cteflonSurface_up_0->DumpInfo();
 //0000000000000000000000000000000000000000000ooooooooooooooooooooooooooo000000000000000000000000ooo

 
  // Definiamo ora le ultime tre superfici ottiche interessanti

  /*
  G4OpticalSurface* opcore_scint = new
  G4OpticalSurface("OpticFiberandScintillator");
    opcore_scint->SetType(dielectric_LUTDAVIS);
    opcore_scint->SetFinish(Polished_LUT);
    opcore_scint->SetModel(DAVIS);

    G4LogicalBorderSurface* core_scint=
            new G4LogicalBorderSurface("opfibscint",
                                   physcore,phys_cryst,opcore_scint);

    G4OpticalSurface* opticalSurface_7 = dynamic_cast <G4OpticalSurface*>
          (core_scint->GetSurface(physcore,phys_cryst)->
                                                         GetSurfaceProperty());
    if (opticalSurface_7) opticalSurface_7->DumpInfo();

     */
  G4OpticalSurface *opcore_clad = new G4OpticalSurface("OpticFiberandClad");
  opcore_clad->SetType(dielectric_LUTDAVIS);
  opcore_clad->SetFinish(Polished_LUT);
  opcore_clad->SetModel(DAVIS);

  G4LogicalBorderSurface *core_clad_up = new G4LogicalBorderSurface(
      "opfibclad_up", physcore, physclad, opcore_clad);
  G4LogicalBorderSurface *core_clad_down = new G4LogicalBorderSurface(
      "opfibclad_down", physclad, physcore, opcore_clad);

  G4OpticalSurface *opticalSurface_8 = dynamic_cast<G4OpticalSurface *>(
      core_clad_up->GetSurface(physcore, physclad)->GetSurfaceProperty());
  if (opticalSurface_8)
    opticalSurface_8->DumpInfo();
/*
  // Ora creaiamo il photodetector

  G4double ephoton[] = {7.0 * eV, 7.14 * eV};
  const G4int num = sizeof(ephoton) / sizeof(G4double);
  G4double photocath_EFF[] = {1., 1.}; // Enables 'detection' of photons
  assert(sizeof(photocath_EFF) == sizeof(ephoton));
  G4double photocath_ReR[] = {1.92, 1.92};
  assert(sizeof(photocath_ReR) == sizeof(ephoton));
  G4double photocath_ImR[] = {1.69, 1.69};
  assert(sizeof(photocath_ImR) == sizeof(ephoton));
  G4MaterialPropertiesTable *photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("EFFICIENCY", ephoton, photocath_EFF, num);
  photocath_mt->AddProperty("REALRINDEX", ephoton, photocath_ReR, num);
  photocath_mt->AddProperty("IMAGINARYRINDEX", ephoton, photocath_ImR, num);

  /*G4OpticalSurface* photocath_opsurf=
    new G4OpticalSurface("photocath_opsurf",glisur,polished,
                         dielectric_dielectric);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

  //****************** Build Photodiode coupling (Air)
  G4double innerRadius_air = 0. * cm;
  G4double startAngle_air = 0. * deg;
  G4double spanningAngle_air = 360. * deg;
  G4double height_air = 0.01 * cm;

  // the "photocathode" is a metal slab at the back of the glass that
  // is only a very rough approximation of the real thing since it only
  // absorbs or detects the photons based on the efficiency set below
  G4Tubs *air_coupling = new G4Tubs(
      "photocath_tube", innerRadius_air, opticfiber_core_diameter / 2,
      height_air / 2, startAngle_air, spanningAngle_air);

  G4RotationMatrix rotm_air = G4RotationMatrix();
  rotm_air.rotateY(0 * deg);
  G4Transform3D transform_air = G4Transform3D(
      rotm_air, G4ThreeVector(0, 0, (opticfiber_core_dx / 2 + height_air / 2)));

  G4LogicalVolume *air_coupling_log =
      new G4LogicalVolume(air_coupling, airNist, "Air_coupling_LV");

  G4VPhysicalVolume *couplingphys =
      new G4PVPlacement(transform_air, air_coupling_log, "couplingphys",
                        opticfiber_core_log, false, 0);

  //****************** Build Photodiode
  G4double innerRadius_pmt = 0. * cm;
  G4double startAngle_pmt = 0. * deg;
  G4double spanningAngle_pmt = 360. * deg;
  G4double height_pmt = 0.1 * cm;

  // the "photocathode" is a metal slab at the back of the glass that
  // is only a very rough approximation of the real thing since it only
  // absorbs or detects the photons based on the efficiency set below
  fPhotocath = new G4Tubs("photocath_tube", innerRadius_pmt,
                          opticfiber_core_diameter / 2, height_pmt / 2,
                          startAngle_pmt, spanningAngle_pmt);

  G4RotationMatrix rotm_pd = G4RotationMatrix();
  rotm_pd.rotateY(0 * deg);
  G4Transform3D transform_pd = G4Transform3D(
      rotm_pd, G4ThreeVector(0, 0, (height_air / 2 + height_pmt / 2)));

  fPhotocath_log = new G4LogicalVolume(fPhotocath, PMMA, "PhotoDiode_LV");

  G4VPhysicalVolume *photophys = new G4PVPlacement(
      transform_pd, fPhotocath_log, "photodiode", air_coupling_log, false, 0);

  //**Create logical skin surfaces
  // new G4LogicalSkinSurface("photocath_surf",fPhotocath_log,photocath_opsurf);

  G4OpticalSurface *opcore_pd = new G4OpticalSurface("OpticFiberandp");
  opcore_pd->SetType(dielectric_LUTDAVIS);
  opcore_pd->SetFinish(Polished_LUT);
  opcore_pd->SetModel(DAVIS);

  G4LogicalBorderSurface *core_pd =
      new G4LogicalBorderSurface("opfibpd", physcore, photophys, opcore_pd);

  G4OpticalSurface *opticalSurface_9 = dynamic_cast<G4OpticalSurface *>(
      core_pd->GetSurface(physcore, photophys)->GetSurfaceProperty());
  if (opticalSurface_9)
    opticalSurface_9->DumpInfo();*/

  //

  

 G4OpticalSurface* opphantom_clad = new G4OpticalSurface("PhantomandClad");
   opphantom_clad->SetType(dielectric_LUTDAVIS);
   opphantom_clad->SetFinish(Rough_LUT);
   opphantom_clad->SetModel(DAVIS);

   G4LogicalBorderSurface* phantom_clad=
           new G4LogicalBorderSurface("opfibclad",
                                  phantom_physical,physclad,opphantom_clad);

   G4OpticalSurface* opticalSurface_9 = dynamic_cast <G4OpticalSurface*>
         (phantom_clad->GetSurface(phantom_physical,physclad)->
                                                        GetSurfaceProperty());
   if (opticalSurface_9) opticalSurface_9->DumpInfo();

         
  G4VisAttributes *white = new G4VisAttributes(G4Colour());
  white->SetVisibility(true);
  // white -> SetForceSolid(true);

  G4VisAttributes *blue = new G4VisAttributes(G4Colour(0., 0., 1.));
  blue->SetVisibility(true);
  // blue -> SetForceSolid(true);

  G4VisAttributes *grey = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  grey->SetVisibility(true);
  // grey-> SetForceSolid(true);

  G4VisAttributes *yellow = new G4VisAttributes(G4Colour(1., 1., 0.));
  yellow->SetVisibility(true);
  // yellow-> SetForceSolid(true);

  G4VisAttributes *skyBlue1 =
      new G4VisAttributes(G4Colour(135 / 255., 206 / 255., 235 / 255.));
  skyBlue1->SetVisibility(true);
  logicCryst->SetVisAttributes(red);
  /* air_coupling_log->SetVisAttributes(yellow);
  fPhotocath_log->SetVisAttributes(blue);
   logicwrapper_long->SetVisAttributes(skyBlue1);
        logicwrapper_side->SetVisAttributes(skyBlue1);
            logicwrapper_little->SetVisAttributes(skyBlue1);*/
  opticfiber_core_log->SetVisAttributes(red);
  opticfiber_clad_log->SetVisAttributes(skyBlue1);
  opticfiber_cladding_log->SetVisAttributes(green);

  G4double maxStep_det = 0.1 * mm;
  fStepLimit = new G4UserLimits(maxStep_det);
  logicCryst->SetUserLimits(fStepLimit);

  logicCryst->SetVisAttributes(green);
  
  logicTeflon->SetVisAttributes(red);



  G4Region *CrystalRegion = new G4Region("crystal_reg");
  logicCryst->SetRegion(CrystalRegion);
  CrystalRegion->AddRootLogicalVolume(logicCryst);

  G4Region *OFcoreRegion = new G4Region("OF_core_reg");
  opticfiber_core_log->SetRegion(OFcoreRegion);
  OFcoreRegion->AddRootLogicalVolume(opticfiber_core_log);

  G4Region *OFcladRegion = new G4Region("OF_clad_reg");
  opticfiber_clad_log->SetRegion(OFcladRegion);
  OFcladRegion->AddRootLogicalVolume(opticfiber_clad_log);

  return physicalTreatmentRoom;
}

////////////////////////////////////////////////////////////////////////////////////////////

void FlashDetectorConstruction::ConstructSDandField() {
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  G4MultiFunctionalDetector *cryst = new G4MultiFunctionalDetector("crystalSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
  G4VPrimitiveScorer *primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
  SetSensitiveDetector("CrystalLV", cryst);

  /*G4String SDname_of_cryst="CrystSD";
   G4SDManager* SDman_of_cryst = G4SDManager::GetSDMpointer();
   sd_of_cryst = new FSensitiveDetector(SDname_of_cryst);


  SDman_of_cryst->AddNewDetector( sd_of_cryst );
  SetSensitiveDetector("CrystalLV",sd_of_cryst);


  G4String SDname_of_of="OpticFiberSD";
   G4SDManager* SDman_of_of = G4SDManager::GetSDMpointer();
   sd_of_of = new OpticFiberSD(SDname_of_of);


  SDman_of_of->AddNewDetector( sd_of_of );
  SetSensitiveDetector("OF_core_LV",sd_of_of);

  G4String SDname_of_pd="PhotoDiodeSD";
   G4SDManager* SDman_of_pd = G4SDManager::GetSDMpointer();
   sd_of_pd = new PhotoDiodeSD(SDname_of_pd);


  SDman_of_pd->AddNewDetector( sd_of_pd );
  SetSensitiveDetector("PhotoDiode_LV",sd_of_pd);*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
