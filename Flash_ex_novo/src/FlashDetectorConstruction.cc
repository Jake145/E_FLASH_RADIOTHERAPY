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

#include "G4RunManager.hh"

#include "G4Region.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Applicator80BeamLine.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashDetectorConstruction::FlashDetectorConstruction()
: G4VUserDetectorConstruction(),physicalTreatmentRoom(0),Collimator(0),
  fCheckOverlaps(true)
{

  DefineMaterials();  
  
  //Construct();
  //ConstructSDandField();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashDetectorConstruction::~FlashDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FlashDetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();
  
  G4bool isotopes = false;
  
  G4Element*  O = man->FindOrBuildElement("O" , isotopes); 
  G4Element* Si = man->FindOrBuildElement("Si", isotopes);
  G4Element* Lu = man->FindOrBuildElement("Lu", isotopes); 
  G4Element* Y = man->FindOrBuildElement("Y", isotopes); 
  
  G4Material* LYSO = new G4Material("Lu2Y2SiO5", 7.3*g/cm3, 4);
  LYSO->AddElement(Lu, 2);
  LYSO->AddElement(Si, 1);
  LYSO->AddElement(O , 5);
  LYSO->AddElement(Y,2);
  
G4cout << "The material of Detector has been changed is LYSO" << G4endl;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FlashDetectorConstruction::ConstructPhantom(){


G4NistManager* nist = G4NistManager::Instance();
G4Material* phantomMaterial = nist->FindOrBuildMaterial("G4_WATER");
G4double phantomSizeX=20.0*cm, phantomSizeY=20.0*cm, phantomSizeZ=20.0*cm;
G4ThreeVector phantomPosition = G4ThreeVector(-99.4*mm,  0.*mm,0.*mm);
// Definition of the solid volume of the Phantom
    phantom = new G4Box("Phantom", 
			phantomSizeX/2, 
			phantomSizeY/2, 
			phantomSizeZ/2);
    
// Definition of the logical volume of the Phantom
     phantomLogicalVolume = new G4LogicalVolume(phantom,	
					     phantomMaterial, 
					     "phantomLog", 0, 0, 0);
  
    // Definition of the physics volume of the Phantom
    new G4PVPlacement(0,
	                                    phantomPosition,
					    "phantomPhys",
					    phantomLogicalVolume,
					    physicalTreatmentRoom,
					    false,
					    0);

G4Region* PhantomRegion = new G4Region("Phantom_reg");
phantomLogicalVolume->SetRegion(PhantomRegion);
PhantomRegion->AddRootLogicalVolume(phantomLogicalVolume);

// Visualisation attributes of the phantom
    red = new G4VisAttributes(G4Colour(0/255., 255/255. ,0/255.));
    red -> SetVisibility(true);
   
    phantomLogicalVolume -> SetVisAttributes(red); 

    
    G4double maxStep = 0.1*mm;
  fStepLimit = new G4UserLimits(maxStep);
  phantomLogicalVolume->SetUserLimits(fStepLimit);
  
  //return physicalTreatmentRoom;
}




G4VPhysicalVolume* FlashDetectorConstruction::Construct()
{
// -----------------------------
  // Treatment room - World volume
  //------------------------------
  // Treatment room sizes
  const G4double worldX = 400.0 *cm;
  const G4double worldY = 400.0 *cm;
  const G4double worldZ = 400.0 *cm;
  G4bool isotopes = false;
 
  G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  G4Box* treatmentRoom = new G4Box("TreatmentRoom",worldX,worldY,worldZ);
  G4LogicalVolume* logicTreatmentRoom = new G4LogicalVolume(treatmentRoom, 
                                                            airNist, 
                                                            "logicTreatmentRoom", 
							    0,0,0);
  physicalTreatmentRoom = new G4PVPlacement(0,
					    G4ThreeVector(),
					    "physicalTreatmentRoom", 
					    logicTreatmentRoom, 
					    0,false,0);
 

  // The treatment room is invisible in the Visualisation
  logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::GetInvisible());
    ConstructPhantom();
 // -----------------------------
  // detector
  //------------------------------

red = new G4VisAttributes(G4Colour(0/255., 255/255. ,0/255.));
    red -> SetVisibility(true);
    green = new G4VisAttributes(G4Colour(255/255., 0./255. ,0/255.));
    green -> SetVisibility(true);
   
    

G4double cryst_dX = 1*cm, cryst_dY = 1*mm, cryst_dZ = 1*mm;
  G4double gap = 0*mm;        //a gap for wrapping, change to add gap
  G4double dX = cryst_dX - gap, dY = cryst_dY - gap, dZ = cryst_dZ - gap ;
// optic fiber
  //
  G4double opticfiber_core_dx= 5*cm;
   G4double opticfiber_core_radius=0.98*mm;
   G4double optic_fiber_clad_radius=1*mm;

G4NistManager* nist = G4NistManager::Instance();

  G4Material* cryst_mat   = nist->FindOrBuildMaterial("Lu2Y2SiO5");
  G4Material* PMMA = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", 
  isotopes);
  G4Material* PE = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYETHYLENE", 
  isotopes);    

G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, dZ/2);
                     
  G4LogicalVolume* logicCryst = 
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name
G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(0*deg); 
    
    G4ThreeVector position = G4ThreeVector(-50.0*mm,  0.*mm,0.*mm);     
    G4Transform3D transform = G4Transform3D(rotm,position);
                                    
    new G4PVPlacement(transform,logicCryst,            
                      "crystalphys",         
                      phantomLogicalVolume,
                      false,0);                  
  G4Tubs* opticfiber_core =
    new G4Tubs("OF_core", 0.*cm, opticfiber_core_radius, opticfiber_core_dx/2, 0., CLHEP::twopi);
      
  G4LogicalVolume* opticfiber_core_log =                         
    new G4LogicalVolume(opticfiber_core,       //its solid
                        PMMA,         //its material
                        "OF_core_LV");         //its name
  G4Tubs* opticfiber_clad =
    new G4Tubs("OF_clad", opticfiber_core_radius,optic_fiber_clad_radius, opticfiber_core_dx/2, 0., twopi);
      
  G4LogicalVolume* opticfiber_clad_log =                         
    new G4LogicalVolume(opticfiber_clad,       //its solid
                        PE,         //its material
                        "OF_clad_LV");         //its name
  G4ThreeVector position_opt = G4ThreeVector(0.*mm ,  0.*mm,0.*mm);  
  G4ThreeVector position_clad = G4ThreeVector(0.*mm + (dX/2 + opticfiber_core_dx/2),  0.*mm,0.*mm);  
   
    
    G4RotationMatrix rotm_opt  = G4RotationMatrix();
    rotm_opt.rotateY(0*deg); 
    G4Transform3D transform_opt = G4Transform3D(rotm_opt,position_opt);
    
    
    G4RotationMatrix rotm_clad  = G4RotationMatrix();
    rotm_clad.rotateY(90*deg); 
    G4Transform3D transform_clad = G4Transform3D(rotm_clad,position_clad);
    
    
  // place rings within detector 
  //
    new G4PVPlacement(transform_clad,                     // rotation
                       //position
                      opticfiber_clad_log,             //its logical volume
                      "outerfiber",                //its name
                      logicCryst,false,0);                 //no boolean operation
                             // checking overlaps 
    new G4PVPlacement(transform_opt,                     // rotation
                           //position
                          opticfiber_core_log,             //its logical volume
                          "OF_clad_LV",                //its name
                          opticfiber_clad_log,false,0);                 //no boolean operation
                                 // checking overlaps 
    
   
    G4VisAttributes * skyBlue1 = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
     skyBlue1->SetVisibility(true);
    logicCryst -> SetVisAttributes(red);
    opticfiber_core_log -> SetVisAttributes(skyBlue1);
    opticfiber_clad_log -> SetVisAttributes(red);

    
   G4double maxStep_det = 0.1*mm;
  fStepLimit = new G4UserLimits(maxStep_det);
  logicCryst->SetUserLimits(fStepLimit);
  
  logicCryst -> SetVisAttributes(green); 

  Collimator = new Applicator80BeamLine(physicalTreatmentRoom);
  
  G4Region* CrystalRegion = new G4Region("crystal_reg");
logicCryst->SetRegion(CrystalRegion);
CrystalRegion->AddRootLogicalVolume(logicCryst);

G4Region* OFcoreRegion = new G4Region("OF_core_reg");
opticfiber_core_log->SetRegion(OFcoreRegion);
OFcoreRegion->AddRootLogicalVolume(opticfiber_core_log);

G4Region* OFcladRegion = new G4Region("OF_clad_reg");
opticfiber_clad_log->SetRegion(OFcladRegion);
OFcladRegion->AddRootLogicalVolume(opticfiber_clad_log);
  return physicalTreatmentRoom;
   }

////////////////////////////////////////////////////////////////////////////////////////////


void FlashDetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  
  // declare crystal as a MultiFunctionalDetector scorer
  //  
  G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystalSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
  SetSensitiveDetector("CrystalLV",cryst);
  
  // declare patient as a MultiFunctionalDetector scorer
  //  
  //G4MultiFunctionalDetector* patient = new G4MultiFunctionalDetector("Phantom");
  //G4SDManager::GetSDMpointer()->AddNewDetector(phantom);
  //G4VPrimitiveScorer* primitiv2 = new G4PSDoseDeposit("dose");
  //patient->RegisterPrimitive(primitiv2);
  //SetSensitiveDetector("phantomLog",phantom);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
