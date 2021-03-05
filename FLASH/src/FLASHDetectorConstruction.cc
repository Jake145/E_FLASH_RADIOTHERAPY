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
// * technical work of the GEANT4 appaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This is the first version of FLASH, a Geant4-based application
//
// 
//////////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "FLASHDetectorConstruction.hh"
#include "G4Tubs.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSDoseDeposit3D.hh"
#include "FLASHDetectorMessenger.hh"

#include "G4SDManager.hh"

#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
/////////////////////////////////////////////////////////////////////////////
FLASHDetectorConstruction::FLASHDetectorConstruction(G4VPhysicalVolume* physicalTreatmentRoom)
:G4VUserDetectorConstruction(),
    motherPhys(physicalTreatmentRoom), // pointer to WORLD volume 
    phantom(0), 
    phantomLogicalVolume(0), 
    phantomPhysicalVolume(0), 
    aRegion(0)
{

 
  


  // Messenger to change parameters of the phantom/detector geometry
  detectorMessenger = new FLASHDetectorMessenger(this);

  // Define here the material of the water phantom and of the detector
  SetPhantomMaterial("G4_WATER"); 

  // Construct geometry (messenger commands)

  // Detector
  // Default crystal sizes
  G4double cryst_dX = 1*cm, cryst_dY = 1*mm, cryst_dZ = 1*mm;
  G4double gap = 0*mm;        //a gap for wrapping, change to add gap
  G4double dX = cryst_dX - gap, dY = cryst_dY - gap, dZ = cryst_dZ - gap ;
// optic fiber
  //
  G4double opticfiber_core_dx= 5*cm;
   G4double opticfiber_core_radius=980*um;
   G4double optic_fiber_clad_radius=1000*um;

 
  // Phantom 
  SetPhantomSize(20. *cm, 20. *cm, 20. *cm);   
  SetPhantomPosition(G4ThreeVector(-99.4 *mm, 0. *mm, 0. *mm)); 
  SetDetectorPosition(G4ThreeVector(-99.4 *mm, 0. *mm, 0. *mm));  

 

  // Write virtual parameters to the real ones and check for consistency      
  //UpdateGeometry(); 
DefineMaterials();
}


/////////////////////////////////////////////////////////////////////////////
FLASHDetectorConstruction::~FLASHDetectorConstruction()
{ 
    delete detectorMessenger;
}

/////////////////////////////////////////////////////////////////////////////
// ConstructPhantom() is the method that reconstuct a water box (called phantom 
// (or water phantom)). 
// A water phantom can be considered a good
// approximation of a an human body. 
////////////////////////////////////////////////////////////////////////////

void FLASHDetectorConstruction::ConstructPhantom()
{

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
    phantomPhysicalVolume = new G4PVPlacement(0,
	                                    phantomPosition,
					    "phantomPhys",
					    phantomLogicalVolume,
					    motherPhys,
					    false,
					    0);

// Visualisation attributes of the phantom
    red = new G4VisAttributes(G4Colour(0/255., 255/255. ,0/255.));
    red -> SetVisibility(true);
    //red -> SetForceSolid(true);
    //red -> SetForceWireframe(true);
    phantomLogicalVolume -> SetVisAttributes(red); 
    //phantomLogicalVolume -> SetVisAttributes(G4VisAttributes::GetInvisible());
}

/////////////////////////////////////////////////////////////////////////////
void FLASHDetectorConstruction::DefineMaterials()
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
  

    
}

void FLASHDetectorConstruction::ConstructDetector()
{
G4bool isotopes = false;
G4NistManager* nist = G4NistManager::Instance();
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* cryst_mat   = nist->FindOrBuildMaterial("Lu2Y2SiO5");
  G4Material* PMMA = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", 
  isotopes);
  G4Material* PE = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYETHYLENE", 
  isotopes);    

G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, cryst_dZ/2);
                     
  G4LogicalVolume* logicCryst = 
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name
G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg); 
    
    G4ThreeVector position = G4ThreeVector(-99.4*mm,  0.*mm,0.*mm);     
    G4Transform3D transform = G4Transform3D(rotm,position);
                                    
    G4VPhysicalVolume* physcryst=new G4PVPlacement(transform,logicCryst,            
                      "crystal",         
                      phantomLogicalVolume,
                      false,0);                  
  G4Tubs* opticfiber_core =
    new G4Tubs("OF_core", 0.*cm, opticfiber_core_radius, opticfiber_core_dx, 0., CLHEP::twopi);
      
  G4LogicalVolume* opticfiber_core_log =                         
    new G4LogicalVolume(opticfiber_core,       //its solid
                        PMMA,         //its material
                        "OF_core_LV");         //its name
  G4Tubs* opticfiber_clad =
    new G4Tubs("OF_clad", opticfiber_core_radius,optic_fiber_clad_radius, opticfiber_core_dx, 0., twopi);
      
  G4LogicalVolume* opticfiber_clad_log =                         
    new G4LogicalVolume(opticfiber_clad,       //its solid
                        PE,         //its material
                        "OF_clad_LV");         //its name
  G4ThreeVector position_opt = G4ThreeVector(-99.4*mm + dX/2,  0.*mm,0.*mm);     
    G4Transform3D transform_opt = G4Transform3D(rotm,position_opt); 
  // place rings within detector 
  //
    G4VPhysicalVolume* physclad=new G4PVPlacement(transform,                     // rotation
                       //position
                      opticfiber_clad_log,             //its logical volume
                      "outerfiber",                //its name
                      phantomLogicalVolume,false,0);                 //no boolean operation
                             // checking overlaps 
    G4VPhysicalVolume* physcore=new G4PVPlacement(transform,                     // rotation
                           //position
                          opticfiber_core_log,             //its logical volume
                          "OF_clad_LV",                //its name
                          opticfiber_clad_log,false,0);                 //no boolean operation
                                 // checking overlaps 
    
   
    G4VisAttributes * skyBlue1 = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
    logicCryst -> SetVisAttributes(skyBlue1);
    opticfiber_core_log -> SetVisAttributes(skyBlue1);
    opticfiber_clad_log -> SetVisAttributes(skyBlue1);
   
  // **************
  // Cut per Region    
  // **************
  
  // A smaller cut is fixed in the phantom to calculate the energy deposit with the
  // required accuracy 
    if (!aRegion)
    {
	aRegion = new G4Region("crystal");
	logicCryst -> SetRegion(aRegion);
	aRegion -> AddRootLogicalVolume(logicCryst);
    }
 G4cout << "The Crystal has been built --- Add a scoring mesh for it  in the GUI if appropriate (similar to the phantom one)" << G4endl;
 
}


/////////////////////////////////////////////////////////////////////////////
/*
void  FLASHDetectorConstruction::ParametersCheck()
{
    // Check phantom/detector sizes & relative position
    if (!IsInside(detectorSizeX, 
		detectorSizeY, 
		detectorSizeZ,
		phantomSizeX,
		phantomSizeY,
		phantomSizeZ,
		detectorToPhantomPosition
		))
      G4Exception("FLASHDetectorConstruction::ParametersCheck()", "FLASH0001", FatalException, "Error: Detector is not fully inside Phantom!");
}*/

/////////////////
// MESSENGERS //
////////////////

G4bool FLASHDetectorConstruction::SetPhantomMaterial(G4String material)
{

    if (G4Material* pMat = G4NistManager::Instance()->FindOrBuildMaterial(material, false) )
    {
	phantomMaterial  = pMat;
	if (phantomLogicalVolume) 
	{
	    phantomLogicalVolume ->  SetMaterial(pMat);

	    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
	    G4cout << "The material of Phantom/Detector has been changed to " << material << G4endl;
	}
    }
    else
    {
	G4cout << "WARNING: material \"" << material << "\" doesn't exist in NIST elements/materials"
	    " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl; 
	G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl; 
	return false;
    }

    return true;
}




/////////////////////////////////////////////////////////////////////////////
void FLASHDetectorConstruction::SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) phantomSizeX = sizeX;
    if (sizeY > 0.) phantomSizeY = sizeY;
    if (sizeZ > 0.) phantomSizeZ = sizeZ;
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void FLASHDetectorConstruction::SetDetectorSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) {cryst_dX = sizeX;}
    if (sizeY > 0.) {cryst_dY = sizeY;}
    if (sizeZ > 0.) {cryst_dZ = sizeZ;}
}
/////////////////////////////////////////////////////////////////////////////

void FLASHDetectorConstruction::SetVoxelSize(G4double , G4double , G4double)
{
    G4cout<< "SetVoxelSize method is not needed anymore " << G4endl;
}
void FLASHDetectorConstruction::SetPhantomPosition(G4ThreeVector pos)
{
    phantomPosition = pos;
}

/////////////////////////////////////////////////////////////////////////////
void FLASHDetectorConstruction::SetDetectorPosition(G4ThreeVector displ)
{
    position = displ;
    position_opt=displ+G4ThreeVector(dX/2,dY/2,dZ/2);
}


////////////////////////////////////////////////////////////////////////////////
/*
void FLASHDetectorConstruction::UpdateGeometry()
{
    ParametersCheck();

    G4GeometryManager::GetInstance() -> OpenGeometry();
    if (phantom)
    {
	phantom -> SetXHalfLength(phantomSizeX/2);
	phantom -> SetYHalfLength(phantomSizeY/2);
	phantom -> SetZHalfLength(phantomSizeZ/2);
	phantomPhysicalVolume -> SetTranslation(phantomPosition);
    }
    else   ConstructPhantom();

    // Get the center of the detector 
    SetDetectorPosition();
    if (solidCryst)
    {
	
	solidCryst -> SetZHalfLength(dZ/2);
	physcryst -> SetTranslation(position);      
    }
    else    ConstructDetector();

   

    // Inform the kernel about the new geometry
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();

    PrintParameters();
}

*/

void FLASHDetectorConstruction::PrintParameters()
{

    G4cout << "The (X,Y,Z) dimensions of the phantom are : (" << 
	G4BestUnit( phantom -> GetXHalfLength()*2., "Length") << ',' << 
	G4BestUnit( phantom -> GetYHalfLength()*2., "Length") << ',' << 
	G4BestUnit( phantom -> GetZHalfLength()*2., "Length") << ')' << G4endl; 
    
    G4cout << "The Z dimension of the crystal are : (" << 
	
	G4BestUnit( solidCryst -> GetZHalfLength()*2., "Length") << ')' << G4endl; 

    G4cout << "Displacement between Phantom and World is: "; 
    G4cout << "DX= "<< G4BestUnit(phantomPosition.getX(),"Length") << 
	"DY= "<< G4BestUnit(phantomPosition.getY(),"Length") << 
	"DZ= "<< G4BestUnit(phantomPosition.getZ(),"Length") << G4endl;
}

void FLASHDetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  
  // declare crystal as a MultiFunctionalDetector scorer
  //  
  G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
  SetSensitiveDetector("CrystalLV",cryst);
  
  
}
