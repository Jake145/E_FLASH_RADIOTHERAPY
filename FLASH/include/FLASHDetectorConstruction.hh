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
// This is the first version of FLASH, a Geant4-based application
// 
//////////////////////////////////////////////////////////////////////////////////////////////


#ifndef FLASHDetectorConstruction_H
#define FLASHDetectorConstruction_H 1

#include "G4Box.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4Tubs.hh"
#include "G4VUserDetectorConstruction.hh"
class G4VPhysicalVolume;
class G4LogicalVolume;
class FLASHDetectorMessenger;

class FLASHDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  FLASHDetectorConstruction(G4VPhysicalVolume*);

  virtual ~FLASHDetectorConstruction();

  //  G4VPhysicalVolume *detectorPhysicalVolume;  aggiunto

private: 

  void ConstructPhantom();
  void ConstructDetector();
  
  void ConstructSensitiveDetector();
  //void ParametersCheck();

    void DefineMaterials();

public: 
    //virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
               
  

// Get detector position relative to WORLD
inline G4ThreeVector GetDetectorToWorldPosition()
  {
    return phantomPosition + position;
  }

/////////////////////////////////////////////////////////////////////////////
// Get displacement between phantom and detector by detector position (center of), phantom (center of) and detector sizes
inline G4ThreeVector GetDetectorToPhantomPosition()
{
    return G4ThreeVector(phantomSizeX/2 - dX/2 + position.getX(),
                         phantomSizeY/2 - dY/2 + position.getY(),
                         phantomSizeZ/2 - dZ/2 + position.getZ()
		          );
}

/////////////////////////////////////////////////////////////////////////////
// Calculate (and set) detector position by displacement, phantom and detector sizes
inline void SetDetectorPosition()
  {
	  // Adjust detector position
	  position.setX(position.getX());
	  position.setY(position.getY());
	  position.setZ(position.getZ());
     
    //G4cout << "*************** DetectorToPhantomPosition " << detectorToPhantomPosition/cm << "\n";
    //G4cout << "*************** DetectorPosition " << detectorPosition/cm << "\n";
  }
/////////////////////////////////////////////////////////////////////////////
/*
inline bool IsInside(G4double detectorX,
		     G4double detectorY,
		     G4double detectorZ,
		     G4double phantomX,
		     G4double phantomY,
		     G4double phantomZ,
		     G4ThreeVector detToPhantomPosition)
{

// Dimensions check... X Y and Z
// Firstly check what dimension we are modifying
	{
	    if (detectorX > phantomX) 
		 {
		    G4cout << "Error: Detector X dimension must be smaller or equal to the correspondent of the phantom" << G4endl;
		    return false;
		 }
	    if ( (phantomX - detectorX) < detToPhantomPosition.getX()) 
	         {
		    G4cout << "Error: X dimension doesn't fit with detector to phantom relative position" << G4endl;
		    return false;
	         }
	}

	{
	    if (detectorY > phantomY) 
		 {
		    G4cout << "Error: Detector Y dimension must be smaller or equal to the correspondent of the phantom" << G4endl;
		    return false;
		 }
	    if ( (phantomY - detectorY) < detToPhantomPosition.getY()) 
	     {
		   G4cout << "Error: Y dimension doesn't fit with detector to phantom relative position" << G4endl;
		   return false;
	     }
	}			 

	{
	    if (detectorZ > phantomZ) 
		 {
		    G4cout << "Error: Detector Z dimension must be smaller or equal to the correspondent of the phantom" << G4endl;
		    return false;
		 }
	    if ( (phantomZ - detectorZ) < detToPhantomPosition.getZ()) 
	     {
		   G4cout << "Error: Z dimension doesn't fit with detector to phantom relative position" << G4endl;
		   return false;
	     }
	}

	return true;
}*/
/////////////////////////////////////////////////////////////////////////////

  G4bool  SetPhantomMaterial(G4String material);
  void SetVoxelSize(G4double sizeX, G4double sizeY, G4double sizeZ);
  void SetDetectorSize(G4double sizeX, G4double sizeY, G4double sizeZ);
  void SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ);
  void SetPhantomPosition(G4ThreeVector);
  void SetDetectorPosition(G4ThreeVector position);
  //void UpdateGeometry();
  

  void PrintParameters();
  //G4LogicalVolume* GetDetectorLogicalVolume(){ return detectorLogicalVolume;}

  

  

private:

  FLASHDetectorMessenger* detectorMessenger; 

  G4VisAttributes* red;

  G4VPhysicalVolume* motherPhys;

  G4Box *phantom ;
  G4Tubs *solidCryst, *opticfiber_core, *opticfiber_clad;
  G4LogicalVolume *phantomLogicalVolume, *logicCryst,*opticfiber_core_log,*opticfiber_clad_log; 
  G4VPhysicalVolume *phantomPhysicalVolume,   *physcryst, *physclad, *physcore;
  
  G4double phantomSizeX; 
  G4double phantomSizeY; 
  G4double phantomSizeZ;

  G4double dX; 
  G4double dY; 
  G4double dZ;

  G4double cryst_dX; 
  G4double cryst_dY; 
  G4double cryst_dZ;

G4double gap;
G4double opticfiber_core_dx;
G4double opticfiber_core_radius; 
G4double optic_fiber_clad_radius;

  G4ThreeVector phantomPosition, position, position_opt; //  phantom center, detector center, detector to phantom relative position

  G4Material *phantomMaterial;
  G4Region* aRegion;
  
  
  
};
#endif



