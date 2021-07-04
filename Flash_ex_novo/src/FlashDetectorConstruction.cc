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
  G4double prelude_density = 7.25 * g / cm3;
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
  G4double PhotonEnergy[] = {1.79 * eV, 1.85 * eV, 1.91 * eV, 1.97 * eV,

                       2.04 * eV, 2.11 * eV, 2.19 * eV, 2.27 * eV,

                       2.36 * eV, 2.45 * eV, 2.56 * eV, 2.67 * eV,

                       2.80 * eV, 2.94 * eV, 3.09 * eV, 3.25 * eV,

                       3.44 * eV, 3.65 * eV, 3.89 * eV, 4.16 * eV};

  const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

  G4double refractiveIndexClad2[] = {
      1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
      1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42
      };

  assert(sizeof(refractiveIndexClad2) == sizeof(PhotonEnergy));
  G4double absClad[] = {
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m};

  assert(sizeof(absClad) == sizeof(PhotonEnergy));

  // Add entries into properties table
  G4MaterialPropertiesTable *mptClad2 = new G4MaterialPropertiesTable();
  mptClad2->AddProperty("RINDEX", PhotonEnergy, refractiveIndexClad2, nEntries);
  mptClad2->AddProperty("ABSLENGTH", PhotonEnergy, absClad, nEntries);

  fPethylene2->SetMaterialPropertiesTable(mptClad2);




  std::vector<G4int> natoms;

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
  G4double photonEnergy_fib[] = {1.79 * eV, 1.85 * eV, 1.91 * eV, 1.97 * eV,

                       2.04 * eV, 2.11 * eV, 2.19 * eV, 2.27 * eV,

                       2.36 * eV, 2.45 * eV, 2.56 * eV, 2.67 * eV,

                       2.80 * eV, 2.94 * eV, 3.09 * eV, 3.25 * eV,

                       3.44 * eV, 3.65 * eV, 3.89 * eV, 4.16 * eV};

  const G4int nEntries_fib = sizeof(photonEnergy_fib) / sizeof(G4double);

  G4double refractiveIndex_fiber[] = {
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49
      };

  assert(sizeof(refractiveIndex_fiber) == sizeof(photonEnergy_fib));

  

  // Add entries into properties table
  G4MaterialPropertiesTable *mpt_fiber = new G4MaterialPropertiesTable();
  mpt_fiber->AddProperty("RINDEX", photonEnergy_fib, refractiveIndex_fiber,
                         nEntries_fib);

  //mpt_fiber->AddProperty("ABSLENGTH", photonEnergy_fib, abs_fiber,nEntries_fib);
  // mptWLSfiber->AddProperty("WLSCOMPONENT",photonEnergy_fib,emissionFib,nEntries_fib);

  PMMA->SetMaterialPropertiesTable(mpt_fiber);
  G4cout << "PMMA G4MaterialPropertiesTable" << G4endl;
  mpt_fiber->DumpTable();

  //PE = nist->FindOrBuildMaterial("G4_POLYETHYLENE", isotopes);

  G4double refractiveIndex_Clad[] = {
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
      1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49
      };

  assert(sizeof(refractiveIndex_Clad) == sizeof(photonEnergy_fib));

  G4double abs_Clad[] = {
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m,
      20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m, 20.0 * m
      };

  assert(sizeof(abs_Clad) == sizeof(photonEnergy_fib));

  // Add entries into properties table
  G4MaterialPropertiesTable *mptClad = new G4MaterialPropertiesTable();
  mptClad->AddProperty("RINDEX", photonEnergy_fib, refractiveIndex_Clad,
                       nEntries_fib);
  mptClad->AddProperty("ABSLENGTH", photonEnergy_fib, abs_Clad, nEntries_fib);

  PE->SetMaterialPropertiesTable(mptClad);

    TEFLON = nist->FindOrBuildMaterial("G4_TEFLON",isotopes);

  G4double photonEnergy_teflon[] =
              {1.79 * eV, 1.85 * eV, 1.91 * eV, 1.97 * eV,

                       2.04 * eV, 2.11 * eV, 2.19 * eV, 2.27 * eV,

                       2.36 * eV, 2.45 * eV, 2.56 * eV, 2.67 * eV,

                       2.80 * eV, 2.94 * eV, 3.09 * eV, 3.25 * eV,

                       3.44 * eV, 3.65 * eV, 3.89 * eV, 4.16 * eV};
    G4double refractiveIndex3[] =
              {1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.35,
               1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.35};
  const G4int nEntries_teflon = sizeof(photonEnergy_teflon)/sizeof(G4double);
    G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();
    myMPT3->AddProperty("RINDEX", photonEnergy_teflon, refractiveIndex3,
  nEntries_teflon);

    G4cout << "Teflon G4MaterialPropertiesTable" << G4endl;
    myMPT3->DumpTable();

    TEFLON->SetMaterialPropertiesTable(myMPT3);
    // EJ-212 Scintillator material
    // Values from http://www.eljentechnology.com/index.php/joomla-overview/what-is-new-in-1-5/64-ej-212
    // and from http://www.eljentechnology.com/images/stories/Technical_Information/Physical%20Constant%20of%20Plastic%20Scintillators.pdf
    // Define the material and its composition
    EJ212 = new G4Material("EJ-212", density=1.023*g/cm3, ncomponents=2);
    EJ212->AddElement(fH, 0.5243407708);
    EJ212->AddElement(fC, 0.4756592292);
    
    
    // Add the material's optical properties
    const G4int n_ej212_entries = 26;
    
    // All the properties are stored as key/value pairs. The keys are photon
    // energies (not wavelengths!) the values are the property for that energy
    G4double ej212_photon_energies[n_ej212_entries] =
    {2.38*eV, 2.41*eV, 2.43*eV, 2.45*eV, 2.48*eV, 2.50*eV, 2.53*eV, 2.55*eV,
        2.58*eV, 2.61*eV, 2.64*eV, 2.66*eV, 2.69*eV, 2.72*eV, 2.75*eV, 2.78*eV,
        2.82*eV, 2.85*eV, 2.88*eV, 2.92*eV, 2.93*eV, 2.95*eV, 2.99*eV, 3.02*eV,
        3.06*eV, 3.10*eV};
    
    // For scintillation you define the relative amplitudes of each photon
    // energy in the emission spectrum. This is normalised to 1.
    G4double ej212_emission_spectrum [n_ej212_entries] =
    {0.005, 0.010, 0.020, 0.050, 0.095, 0.100, 0.100, 0.080, 0.070, 0.070,
        0.067, 0.055, 0.045, 0.036, 0.031, 0.027, 0.023, 0.020, 0.018, 0.016,
        0.014, 0.012, 0.010, 0.009, 0.008, 0.007};
    
    // Refractive index; this is static for most photon energies.
    G4double ej212_refractive_index [n_ej212_entries] =
        {1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58,
         1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58,
         1.58, 1.58, 1.58, 1.58};
    
    // Absorption/Attenuation length (distance for signal to drop to 1/e intial)
    G4double ej212_absorption_length [n_ej212_entries] =
       {250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 250*cm,
        250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 250*cm,
        250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 250*cm};
    
    G4MaterialPropertiesTable* ej212_MPT = new G4MaterialPropertiesTable();
    ej212_MPT->AddProperty("RINDEX",                // The property name
                           ej212_photon_energies,   // photon energies (keys)
                           ej212_refractive_index,  // value for that photon energy
                           n_ej212_entries);        // number of key/value pairs
    
    // Add the scintillation emission spectrum (normalised to 1)
    ej212_MPT->AddProperty("FASTCOMPONENT",
                           ej212_photon_energies,
                           ej212_emission_spectrum,
                           n_ej212_entries);
    
    // Specify the absorption lengths
    ej212_MPT->AddProperty("ABSLENGTH",
                           ej212_photon_energies,
                           ej212_absorption_length,
                           n_ej212_entries);
    
    // How many photons / MeV / Step, this is assumed to be linear
    ej212_MPT->AddConstProperty("SCINTILLATIONYIELD", 10000.0/MeV);
    // Scale the variance of the number of photons / MeV / Step (e.g. due to
    // doping of the material and/or Fano factor)
    ej212_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
    // If EJ-212 had a
    ej212_MPT->AddConstProperty("YIELDRATIO", 1.0);
    // The time constant for the time evolution of the number of scintillation
    // photons (i.e. tau in the formula 'A*exp(-t/tau)' )
    ej212_MPT->AddConstProperty("FASTTIMECONSTANT", 2.4*ns);
    // There is also a 'FASTSCINTILLATIONRISETIME' that accounts for the rise
    // time of the scintillator but given that for EJ-212 this is 0.9ns
    // it will be ignored.
    
    // Set the material table
    EJ212->SetMaterialPropertiesTable(ej212_MPT);
}


G4VPhysicalVolume * FlashDetectorConstruction::ConstructPhantom_Support(G4double CollPos, G4double Cx,G4double Cy,G4double Cz,G4double d,G4double Oz,G4bool plastic_bool) {

  G4Material *phantomMaterial = nist->FindOrBuildMaterial("G4_WATER");

  // ------------ Generate & Add Material Properties Table ------------
  //

  G4double photonEnergy[] = {1.79 * eV, 1.85 * eV, 1.91 * eV, 1.97 * eV,

                       2.04 * eV, 2.11 * eV, 2.19 * eV, 2.27 * eV,

                       2.36 * eV, 2.45 * eV, 2.56 * eV, 2.67 * eV,

                       2.80 * eV, 2.94 * eV, 3.09 * eV, 3.25 * eV,

                       3.44 * eV, 3.65 * eV, 3.89 * eV, 4.16 * eV};

  const G4int nEntries = sizeof(photonEnergy) / sizeof(G4double);

  //
  // Water
  //
  G4double refractiveIndex1[] = {
      1.35,   1.3505, 1.351,  1.3518,
      1.3522, 1.3530, 1.3535, 1.354,  1.3545, 1.355,  1.3555, 1.356,
      1.3568, 1.3572, 1.358,  1.3585, 1.359,  1.3595, 1.36,   1.3608};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

 
  G4MaterialPropertiesTable *myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX", photonEnergy, refractiveIndex1, nEntries)
      ->SetSpline(true);
 
  phantomMaterial->SetMaterialPropertiesTable(myMPT1);
  
 G4double Position_coefficient=  CollPos;
  G4double phantomSizeX = 4 * mm, phantomSizeY = 30.0 * cm,
           phantomSizeZ = 30.0 * cm, phantom_coordinateX=(Position_coefficient * mm + phantomSizeX/2);
  
  
    G4ThreeVector phantomPosition = G4ThreeVector(phantom_coordinateX , 0. * mm, 0. * mm);
  // Definition of the solid volume of the Phantom
  phantom = new G4Box("Phantom", phantomSizeX / 2, phantomSizeY / 2,
                      phantomSizeZ / 2);

  // Definition of the logical volume of the Phantom
  phantomLogicalVolume =
      new G4LogicalVolume(phantom, PMMA, "phantomLog", 0, 0, 0);

  // Definition of the physics volume of the Phantom
  phant_phys =
      new G4PVPlacement(0, phantomPosition, "phantomPhys", phantomLogicalVolume,physicalTreatmentRoom, false, 0);
  
    
    
  //============================DETECTOR_SUPPORT=====================================//
  
  //G4double support_x=1*cm, wedge_X=1*mm,wedge_Y=1*mm,wedge_Z=1*cm;
  if (plastic_bool == false){
   support_x=2*cm;
   G4double wedge_X=Cz+2*d,wedge_Y=Cy+2*d,wedge_Z=Cx+Oz+2*d; 
  G4ThreeVector xTrans(-(support_x/2-wedge_X/2), 0, -wedge_Z/2);
  G4Box* support_whole= new G4Box("Support_w",support_x/2,10*cm,10*cm);
  G4Box* wedge = new G4Box("Wedge",wedge_X/2,wedge_Y/2,wedge_Z/2);
  //G4SubtractionSolid* support_wedged = new G4SubtractionSolid("Wedged_Support",support_whole,wedge,0,xTrans);
  
  DetectorSupport =
      new G4LogicalVolume(support_whole, PMMA, "SupportLog");

G4RotationMatrix rotmp = G4RotationMatrix();
  rotmp.rotateY(0 * deg);
   supp_coordinateX=phantom_coordinateX+(phantomSizeX/2+support_x/2);

  G4ThreeVector positionp = G4ThreeVector(supp_coordinateX,0,Cx/2);
  G4Transform3D transformp = G4Transform3D(rotmp, positionp);
  
  new G4PVPlacement(
      transformp, DetectorSupport, "supportphys", logicTreatmentRoom, false, 0);
      
  AirBox = new G4LogicalVolume(wedge, airNist, "filler");
  
    G4Transform3D transformab = G4Transform3D(rotmp, xTrans);
    
  new G4PVPlacement(transformab, AirBox, "f", DetectorSupport, false, 0);}
  
  else if (plastic_bool == true){
  
   support_x=2*cm;
  G4Box* support_whole= new G4Box("Support_w",support_x/2,10*cm,10*cm);
  
  DetectorSupport =
      new G4LogicalVolume(support_whole, PMMA, "SupportLog");

G4RotationMatrix rotmp = G4RotationMatrix();
  rotmp.rotateY(0 * deg);
   supp_coordinateX=phantom_coordinateX+(phantomSizeX/2+support_x/2);

  G4ThreeVector positionp = G4ThreeVector(supp_coordinateX,0,0);
  G4Transform3D transformp = G4Transform3D(rotmp, positionp);
  
  new G4PVPlacement(
      transformp, DetectorSupport, "supportphys", logicTreatmentRoom, false, 0);
      
  
  }
  //=================================================================================//
  //================Second Piece of Phantom==================================//
  G4double phantom2_coordinateX= supp_coordinateX +(15 * cm -(phantomSizeX / 2 +support_x/2))+support_x/2;
  G4ThreeVector phantomPosition_2 = G4ThreeVector(phantom2_coordinateX , 0. * mm, 0. * mm);
  // Definition of the solid volume of the Phantom
  G4Box * phantom_2 = new G4Box("Phantom_2", 15 * cm -(phantomSizeX / 2 +support_x/2) , phantomSizeY / 2,
                      phantomSizeZ / 2);

  // Definition of the logical volume of the Phantom
  G4LogicalVolume* phantomLogicalVolume_2 =
      new G4LogicalVolume(phantom_2, PMMA, "phantomLog_2", 0, 0, 0);

  // Definition of the physics volume of the Phantom
 
      new G4PVPlacement(0, phantomPosition_2, "phantomPhys_2", phantomLogicalVolume_2,
                        physicalTreatmentRoom, false, 0);
  
  
  
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
G4VPhysicalVolume *FlashDetectorConstruction::ConstructPhantom(G4double CollPos){


G4Material *phantomMaterial = nist->FindOrBuildMaterial("G4_WATER");

  // ------------ Generate & Add Material Properties Table ------------
  //

  G4double photonEnergy[] = {1.79 * eV, 1.85 * eV, 1.91 * eV, 1.97 * eV,

                       2.04 * eV, 2.11 * eV, 2.19 * eV, 2.27 * eV,

                       2.36 * eV, 2.45 * eV, 2.56 * eV, 2.67 * eV,

                       2.80 * eV, 2.94 * eV, 3.09 * eV, 3.25 * eV,

                       3.44 * eV, 3.65 * eV, 3.89 * eV, 4.16 * eV};

  const G4int nEntries = sizeof(photonEnergy) / sizeof(G4double);

  //
  // Water
  //
  G4double refractiveIndex1[] = {
      1.35,   1.3505, 1.351,  1.3518,
      1.3522, 1.3530, 1.3535, 1.354,  1.3545, 1.355,  1.3555, 1.356,
      1.3568, 1.3572, 1.358,  1.3585, 1.359,  1.3595, 1.36,   1.3608};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

 
  G4MaterialPropertiesTable *myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX", photonEnergy, refractiveIndex1, nEntries)
      ->SetSpline(true);
  
  phantomMaterial->SetMaterialPropertiesTable(myMPT1);
  
G4double Position_coefficient= CollPos;
  G4double phantomSizeX = 30.0 * cm, phantomSizeY = 30.0 * cm,
           phantomSizeZ = 30.0 * cm, phantom_coordinateX=(Position_coefficient * mm + phantomSizeX/2);
  //G4ThreeVector phantomPosition = G4ThreeVector(-(199.4 * mm - phantomSizeX/2) , 0. * mm, 0. * mm);
  
    G4ThreeVector phantomPosition = G4ThreeVector(phantom_coordinateX , 0. * mm, 0. * mm);
  // Definition of the solid volume of the Phantom
  phantom = new G4Box("Phantom", phantomSizeX / 2, phantomSizeY / 2,
                      phantomSizeZ / 2);

  // Definition of the logical volume of the Phantom
  phantomLogicalVolume =
      new G4LogicalVolume(phantom, PMMA, "phantomLog", 0, 0, 0);

  // Definition of the physics volume of the Phantom
  phant_phys =
      new G4PVPlacement(0, phantomPosition, "phantomPhys", phantomLogicalVolume,
                        physicalTreatmentRoom, false, 0);
  
    
    G4Region *PhantomRegion = new G4Region("Phantom_reg");
  phantomLogicalVolume->SetRegion(PhantomRegion);
  PhantomRegion->AddRootLogicalVolume(phantomLogicalVolume);

  // Visualisation attributes of the phantom
  red = new G4VisAttributes(G4Colour(0 / 255., 255 / 255., 0 / 255.));
  red->SetVisibility(true);

blue = new G4VisAttributes(G4Colour(0 / 255., 0./ 255., 255. / 255.));
  blue->SetVisibility(true);

  phantomLogicalVolume->SetVisAttributes(red);


  G4double maxStep = 0.1 * mm;
  fStepLimit = new G4UserLimits(maxStep);
  phantomLogicalVolume->SetUserLimits(fStepLimit);

  return phant_phys;


}
G4VPhysicalVolume *FlashDetectorConstruction::BuildDetector(G4double dX,G4double dY,G4double dZ,G4double fPTFEThickness,G4double opticfiber_core_dx,G4bool plastic_bool){

fCheckOverlaps = true;

G4double maxStep_det = 0.1 * mm;
  fStepLimit = new G4UserLimits(maxStep_det);
  
 green = new G4VisAttributes(G4Colour(0 / 255., 255 / 255., 0 / 255.));
  green->SetVisibility(true);
  red = new G4VisAttributes(G4Colour(255 / 255., 0. / 255., 0 / 255.));
  red->SetVisibility(true);
  G4VisAttributes *white = new G4VisAttributes(G4Colour());
  white->SetVisibility(true);
  // white -> SetForceSolid(true);

  blue = new G4VisAttributes(G4Colour(0., 0., 1.));
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
 //fiber dim
if( plastic_bool == false){
  G4double opticfiber_core_diameter = 0.98 * mm;
  G4double optic_fiber_clad_diameter = 2.2 * mm;
  G4double optic_fiber_cladding_diameter = 1. * mm;
  //Teflon dim

   //constuct detector
  G4Box *solidCryst = new G4Box("crystal", dX / 2, dY / 2, dZ / 2);
  //
  logicCryst = new G4LogicalVolume(solidCryst,   // its solid
                                   cryst_mat,    // its material
                                   "CrystalLV"); // its name
  G4RotationMatrix rotm = G4RotationMatrix();
  rotm.rotateY(90 * deg);

  G4ThreeVector position = G4ThreeVector(0, 0, (opticfiber_core_dx/2));
  G4Transform3D transform = G4Transform3D(rotm, position);

  phys_cryst = new G4PVPlacement(
      transform, logicCryst, "crystalphys", AirBox, false, 0,fCheckOverlaps);
     

  
  ////////////////////

  G4Tubs *opticfiber_cladding = new G4Tubs(
      "OF_cladding", opticfiber_core_diameter / 2,
      optic_fiber_cladding_diameter / 2, opticfiber_core_dx / 2, 0., twopi);

  ////////////////////

  G4Tubs *opticfiber_core =
      new G4Tubs("OF_core", 0. * cm, opticfiber_core_diameter / 2,
                 opticfiber_core_dx / 2, 0., CLHEP::twopi);

  opticfiber_core_log = new G4LogicalVolume(opticfiber_core, // its solid
                                            PMMA,            // its material
                                            "OF_core_LV");   // its name
  G4Tubs *opticfiber_clad = new G4Tubs(
      "OF_clad", optic_fiber_cladding_diameter / 2,
      optic_fiber_clad_diameter / 2,
      
      opticfiber_core_dx / 2, 0., twopi);

  G4LogicalVolume *opticfiber_cladding_log =
      new G4LogicalVolume(opticfiber_cladding, // its solid
                          fPethylene2,         // its material
                          "OF_cladding_LV");   // its name

  G4LogicalVolume *opticfiber_clad_log =
      new G4LogicalVolume(opticfiber_clad, // its solid
                          PE,              // its material
                          "OF_clad_LV");   // its name
  G4ThreeVector position_opt = G4ThreeVector(0. * mm, 0. * mm, -dX/2);
  //G4ThreeVector position_clad = G4ThreeVector(
  //    0. * mm + (dX / 2 + opticfiber_core_dx / 2), 0. * mm, 0. * mm);
	G4ThreeVector position_clad=position_opt;
  G4RotationMatrix rotm_opt = G4RotationMatrix();
  rotm_opt.rotateX(0 * deg);
  G4Transform3D transform_opt = G4Transform3D(rotm_opt, position_opt);

  G4RotationMatrix rotm_clad = G4RotationMatrix();
  rotm_clad.rotateX(0 * deg);
  G4Transform3D transform_clad = G4Transform3D(rotm_clad, position_clad);

  //
  new G4PVPlacement(
      transform_clad,
      // position
      opticfiber_clad_log, "outerfiber", AirBox, false, 0,fCheckOverlaps);

  new G4PVPlacement(transform_clad,
                                                  // position
                                                  opticfiber_core_log,

                                                  "OF_core_phys",

                                                  AirBox, false, 0,fCheckOverlaps);
  new G4PVPlacement(transform_clad,
                                                      // position
                                                      opticfiber_cladding_log,

                                                      "OF_cladding_phys",

                                                      AirBox, false, 0,fCheckOverlaps);

  G4VSolid * t1= new G4Box("t1",dX/2+fPTFEThickness, dY/2+fPTFEThickness, dZ/2+fPTFEThickness);
  G4VSolid *t2 = new G4Box("t2",dX/2 , dY/2 , dZ/2 );
G4RotationMatrix rotm_t2 = G4RotationMatrix();
  rotm_t2.rotateX(0 * deg);
  G4ThreeVector zTrans(0, 0,0);

  G4SubtractionSolid* solidFullTeflon = new G4SubtractionSolid("hollowteflon", t1, t2, 0, zTrans);


  G4VSolid *tHole = new G4Tubs("sTeflonHole",
			       0.0*cm,
			       optic_fiber_clad_diameter  / 2 + 0.01*mm,
			       fPTFEThickness/2,
			       0.*deg,
			       360.*deg);
  G4RotationMatrix*yRot =new G4RotationMatrix;


  //yRot->rotateX(90*deg);
  yRot->rotateY(270*deg);
  G4ThreeVector zHole((dX/2 + fPTFEThickness/2), 0,0  );

  G4SubtractionSolid* solidHoleTeflon = new G4SubtractionSolid("hollowHoleteflon", solidFullTeflon, tHole, yRot, zHole);
 G4LogicalVolume *logicTeflon = new G4LogicalVolume(solidHoleTeflon, TEFLON,"lTeflon",0,0,0);
 new G4PVPlacement(yRot, G4ThreeVector(0, 0, (opticfiber_core_dx/2)), logicTeflon, "pTeflon", AirBox, false, 0,fCheckOverlaps);
  
  

        

 
  opticfiber_core_log->SetVisAttributes(red);
  opticfiber_clad_log->SetVisAttributes(skyBlue1);
  opticfiber_cladding_log->SetVisAttributes(green);
  logicTeflon->SetVisAttributes(red);
   
  G4Region *OFcoreRegion = new G4Region("OF_core_reg");
  opticfiber_core_log->SetRegion(OFcoreRegion);
  OFcoreRegion->AddRootLogicalVolume(opticfiber_core_log);

  G4Region *OFcladRegion = new G4Region("OF_clad_reg");
  opticfiber_clad_log->SetRegion(OFcladRegion);
  OFcladRegion->AddRootLogicalVolume(opticfiber_clad_log);
  }
  
  else if (plastic_bool==true){
  
  G4Box *solidCryst = new G4Box("crystal", dX / 2, dY / 2, dZ / 2);
  //
  logicCryst = new G4LogicalVolume(solidCryst,   // its solid
                                   EJ212,    // its material
                                   "CrystalLV"); // its name
  G4RotationMatrix rotm = G4RotationMatrix();
  rotm.rotateY(90 * deg);

  G4ThreeVector position = G4ThreeVector(-5*mm, 0, 0);
  G4Transform3D transform = G4Transform3D(rotm, position);

  phys_cryst = new G4PVPlacement(
      transform, logicCryst, "crystalphys", DetectorSupport, false, 0,fCheckOverlaps);
     

  
  
  
  }
  
   logicCryst->SetVisAttributes(red);
  logicCryst->SetUserLimits(fStepLimit);

  logicCryst->SetVisAttributes(green);
  




  G4Region *CrystalRegion = new G4Region("crystal_reg");
  logicCryst->SetRegion(CrystalRegion);
  CrystalRegion->AddRootLogicalVolume(logicCryst);
  return phys_cryst;
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

  airNist =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  // Air
  //
  

  G4Box *treatmentRoom = new G4Box("TreatmentRoom", worldX, worldY, worldZ);
  logicTreatmentRoom = new G4LogicalVolume(
      treatmentRoom, airNist, "logicTreatmentRoom", 0, 0, 0);
  physicalTreatmentRoom =
      new G4PVPlacement(0, G4ThreeVector(), "physicalTreatmentRoom",
                        logicTreatmentRoom, 0, false, 0);

  // The treatment room is invisible in the Visualisation
  logicTreatmentRoom->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  // -----------------------------
  // detector + phantom +Default dimensions
  //------------------------------
  select_EJ212 = true;
//LYSO_dimensions/////////////////////
if (select_EJ212 == false){
 G4double cryst_dX = 1 * cm, cryst_dY = 2 * mm, cryst_dZ = 2 * mm;
   opticfiber_core_dx_ = 5 * cm;
      fPTFEThickness_ = 0.5 * mm;
  G4double gap = 0 * mm; // a gap for wrapping, change to add gap
    dX_ = cryst_dX - gap;
     dY_ = cryst_dY - gap;
      dZ_ = cryst_dZ - gap;
//EJ212_dimensions/////////////////////////

}
else if (select_EJ212==true){
 dX_=2 * cm;
 dZ_ = 2 * mm;
 dY_= 1*cm;
   opticfiber_core_dx_ = 0 *mm;
    fPTFEThickness_ = 0*mm;
}

  //construct collimatore
    Collimator = new Applicator80BeamLine(physicalTreatmentRoom);
  // constuct phantom//////
  
  //phantom_physical=ConstructPhantom(Collimator->finalApplicatorXPositionFlash + Collimator->hightFinalApplicatorFlash);
  phantom_physical=ConstructPhantom_Support(Collimator->finalApplicatorXPositionFlash + Collimator->hightFinalApplicatorFlash,dX_,dY_,dZ_,fPTFEThickness_,opticfiber_core_dx_,select_EJ212);
///construct Detector//////
  detector_physical=BuildDetector(dX_,dY_,dZ_,fPTFEThickness_,opticfiber_core_dx_,select_EJ212);
  

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

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
