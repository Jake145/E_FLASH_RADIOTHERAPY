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

#include "VHEE_collimator.hh"
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

VHEE_collimator::VHEE_collimator(G4VPhysicalVolume *physicalTreatmentRoom)
    : motherPhys(physicalTreatmentRoom), X_center(0.), y_center(0.),z_center(0.),radious_field(0.),number_of_leafs(0.)
       {
  ConstructCollimator();
}

VHEE_collimator::~VHEE_collimator() {}

void VHEE_collimator::ConstructCollimator() {
  // Sets default geometry and materials
  SetMaterial();

  // Construct the whole collimator Beam Line
  ConstructColl();
}

void VHEE_collimator::SetMaterial() {
  
  X_center = -25 *cm; 
  y_center = 0.*cm;
  z_center = 0.*cm;
  radious_field = 3*cm; 
  number_of_leafs = 10;

  blue = new G4VisAttributes(G4Colour(0., 0., 1.));
  blue->SetVisibility(true);
  G4NistManager *man = G4NistManager::Instance();
  G4bool isotopes = false;
  G4int ncomponents;
  G4double CSI_density = 3.16 * g / cm3;
  Bromuro_SI = new G4Material("Silicon_carbide", CSI_density, ncomponents = 2);
  Bromuro_SI->AddElement(man->FindOrBuildElement("C"), 50 * perCent);
  Bromuro_SI->AddElement(man->FindOrBuildElement("Si"), 50 * perCent);
  

 
}

void VHEE_collimator::ConstructColl() {

  
  // Components of the Collimator
  G4bool checkOverlaps = true;
  G4double X_leaf = 7.*cm, Y_leaf = 15.*cm;
  G4double h = radious_field/ number_of_leafs;
  G4Box *wall = new G4Box("wall", X_leaf /2,Y_leaf  , Y_leaf/2);
  G4Box *leaf = new G4Box("leaf", X_leaf / 2, Y_leaf / 2, h / 2);

    LeafLV = new G4LogicalVolume(leaf, Bromuro_SI, "Leaf_LV");
    
    
    
    G4double phi1 = 0. * deg;

  	G4RotationMatrix rm1;
  	rm1.rotateY(phi1);
  
    G4int j = 0;
    
     for (G4int i=0;i<number_of_leafs+1;i++)
  {
    
    G4double Z_ =  (z_center + radious_field - i * h);
    G4double Y_ = y_center + std::sqrt(radious_field * radious_field - (Z_ - z_center)*(Z_ - z_center));
    
  	
    new G4PVPlacement(G4Transform3D(rm1, 
    G4ThreeVector(X_center, Y_ + Y_leaf/2, Z_)),"leaf_phys",LeafLV, motherPhys,false,j,checkOverlaps);
    j++;
    
    new G4PVPlacement(G4Transform3D(rm1, 
    G4ThreeVector(X_center, Y_ + Y_leaf/2, - Z_)),"leaf_phys",LeafLV, motherPhys,false,j,checkOverlaps);
    j++;
   
   new G4PVPlacement(G4Transform3D(rm1, 
    G4ThreeVector(X_center, - (Y_ + Y_leaf/2), - Z_)),"leaf_phys",LeafLV, motherPhys,false,j,checkOverlaps);
    j++;
    
    new G4PVPlacement(G4Transform3D(rm1, 
    G4ThreeVector(X_center, -(Y_ + Y_leaf/2), Z_)),"leaf_phys",LeafLV, motherPhys,false,j,checkOverlaps);
    j++;
  } 

LeafLV->SetVisAttributes(blue);

G4LogicalVolume *Wall_lv = new G4LogicalVolume(wall, Bromuro_SI, "wall_LV");

new G4PVPlacement(G4Transform3D(rm1, 
    G4ThreeVector(X_center, y_center, z_center + radious_field + Y_leaf/2 + h/2)),"wall_phys",Wall_lv, motherPhys,false,0,checkOverlaps);

new G4PVPlacement(G4Transform3D(rm1, 
    G4ThreeVector(X_center, y_center, z_center - radious_field - Y_leaf/2 - h/2)),"wall_phys",Wall_lv, motherPhys,false,1,checkOverlaps);
    
    Wall_lv->SetVisAttributes(blue);

}



