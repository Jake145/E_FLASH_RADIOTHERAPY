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
// This is the first version of FLASH, a Geant4-based application
//
// 
//////////////////////////////////////////////////////////////////////////////////////////////



#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"  
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "G4NistElementBuilder.hh"
#include "G4SubtractionSolid.hh"   
#include "FLASHDetectorConstruction.hh" 
#include "Applicator80BeamLine.hh"
#include "Applicator80BeamLineMessenger.hh" 

Applicator80BeamLine::Applicator80BeamLine():
  physicalTreatmentRoom(0),flashDetectorConstruction(0),
  
  
  solidFirstApplicatorFLASH(0),
  physiFirstApplicatorFLASH(0),
  solidFinalApplicatorFLASH(0),
  physiFinalApplicatorFLASH(0),

  solidGiunz1FinalAppFLASH(0),
  physiGiunz1FinalAppFLASH(0),

  solidGiunz2FinalAppFLASH(0),  
  physiGiunz2FinalAppFLASH(0),

  solidGiunz3FinalAppFLASH(0),
  physiGiunz3FinalAppFLASH(0),  

  solidGiunz3FinalAppIntFLASH(0),
  physiGiunz3FinalAppIntFLASH(0), 
  
  solidGiunz4FinalAppFLASH(0),
  physiGiunz4FinalAppFLASH(0),

  solidGiunz5FinalAppFLASH(0),
  physiGiunz5FinalAppFLASH(0),

  protector1(0),
  physiprotector1(0),
  
  protector2(0),
  physiprotector2(0),
  
  protector3(0),
  physiprotector3(0),
  
  protector4(0),
  physiprotector4(0),
  
  cover1(0),
  physicover1FLASH(0),
  
  cover2(0),
  physicover2FLASH(0),
  
  cover3(0),
  physicover3FLASH(0),
  
  
  /*
  solidBlocco1FLASH(0),
  physiBlocco1FLASH(0),

  solidBlocco2FLASH(0),
  physiBlocco2FLASH(0),

  solidBlocco3FLASH(0),
  physiBlocco3FLASH(0),

  solidBlocco20mmFLASH(0),
  physiBlocco20mmFLASH(0),  
  
  solidCM1_1_2FLASH(0),
  physiCM1_1_2FLASH(0),
  
  solidCM1_2_2FLASH(0),
  physiCM1_2_2FLASH(0),
  
  solidCM2_1_2FLASH(0),
  physiCM2_1_2FLASH(0),

  solidCM2_2_2FLASH(0),
  physiCM2_2_2FLASH(0),  

  solidCCMFLASH(0),
  physiCCMFLASH(0),

  solidPFS1FLASH(0),
  physiPFS1FLASH(0),

  solidPFS2FLASH(0),
  physiPFS2FLASH(0),

  solidPFS3FLASH(0),
  physiPFS3FLASH(0),
  */
  solidFTFLASH(0),
  physiFTFLASH(0)


{
  // Messenger to change parameters of the applicator80BeamLine geometry
  //applicatorMessenger = new Applicator80BeamLineMessenger(this);

}

Applicator80BeamLine::~Applicator80BeamLine()
{
  
  delete flashDetectorConstruction;
}




G4VPhysicalVolume* Applicator80BeamLine::Construct()
{ 
  // Sets default geometry and materials
  SetDefaultDimensions();
  
  // Construct the whole Applicator Beam Line 
  ConstructApplicator80BeamLine();

  
  // FLASHDetectorConstruction builds ONLY the phantom and the detector with its associated ROGeometry
  flashDetectorConstruction = new FLASHDetectorConstruction(physicalTreatmentRoom); 
  
  return physicalTreatmentRoom;
}


// DEFAULT MATERIAL ARE ALSO PROVIDED	
// and COLOURS ARE ALSO DEFINED
// ----------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////
void Applicator80BeamLine::SetDefaultDimensions()
{

   // Set of coulors that can be used
  white = new G4VisAttributes( G4Colour());
  white -> SetVisibility(true);
  //white -> SetForceSolid(true);
	
  blue = new G4VisAttributes(G4Colour(0. ,0. ,1.));
  blue -> SetVisibility(true);
  //blue -> SetForceSolid(true);
	
  gray = new G4VisAttributes( G4Colour(0.5, 0.5, 0.5 ));
  gray-> SetVisibility(true);
  //gray-> SetForceSolid(true);
	
  red = new G4VisAttributes(G4Colour(1. ,0. ,0.));
  red-> SetVisibility(true);
  //red-> SetForceSolid(true);
	
  yellow = new G4VisAttributes(G4Colour(1., 1., 0. ));
  yellow-> SetVisibility(true);
  //yellow-> SetForceSolid(true);
	
  green = new G4VisAttributes( G4Colour(25/255. , 255/255. ,  25/255. ));
  green -> SetVisibility(true);
  //green -> SetForceSolid(true);
	
  darkGreen = new G4VisAttributes( G4Colour(0/255. , 100/255. ,  0/255. ));
  darkGreen -> SetVisibility(true);
  //darkGreen -> SetForceSolid(true);
		
  darkOrange3 = new G4VisAttributes( G4Colour(205/255. , 102/255. ,  000/255. ));
  darkOrange3 -> SetVisibility(true);
  //darkOrange3 -> SetForceSolid(true);
	
  skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
  skyBlue -> SetVisibility(true);
  //skyBlue -> SetForceSolid(true);
  magenta = new G4VisAttributes( G4Colour(255/255. , 0/255. ,  255/255. ));
  magenta -> SetVisibility(true);

  // Geometry FIRST APPLICATOR DEFAULTS

  G4double defaultOuterRadiusFirstApplicatorFLASH = 55. * mm;
  OuterRadiusFirstApplicatorFLASH = defaultOuterRadiusFirstApplicatorFLASH;

  G4double defaultinnerRadiusFirstApplicatorFLASH = 50. * mm;
  innerRadiusFirstApplicatorFLASH = defaultinnerRadiusFirstApplicatorFLASH;

  // Geometry FINAL APPLICATOR DEFAULTS

  G4double defaultOuterRadiusFinalApplicatorFLASH = 55. *mm;
  OuterRadiusFinalApplicatorFLASH = defaultOuterRadiusFinalApplicatorFLASH;

  G4double defaultinnerRadiusFinalApplicatorFLASH = 50. *mm;
  innerRadiusFinalApplicatorFLASH = defaultinnerRadiusFinalApplicatorFLASH;

  // DEFAULT DEFINITION OF THE MATERIALS
  // All elements and compound definition follows the NIST database
 
  // ELEMENTS
  G4bool isotopes = false;
  G4Material* aluminumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al", isotopes);
  //G4Material* tantalumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ta", isotopes); 
  //G4Material* copperNistAsMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu", isotopes);
  G4Element* zincNist = G4NistManager::Instance()->FindOrBuildElement("Zn");
  G4Element* copperNist = G4NistManager::Instance()->FindOrBuildElement("Cu");

  // COMPOUND
  //G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  //G4Material* kaptonNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON", isotopes);
  G4Material* galacticNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic", isotopes);
  G4Material* PMMANist = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);
  //G4Material* mylarNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR", isotopes);
  G4Material* titanioNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ti", isotopes); 
  
  
  G4double d; // Density
  G4int nComponents;// Number of components 
  G4double fractionmass; // Fraction in mass of an element in a material

  d = 8.40*g/cm3;   // brass 
  nComponents = 2;
  G4Material* brass = new G4Material("Brass", d, nComponents);  
  brass -> AddElement(zincNist, fractionmass = 30 *perCent);
  brass -> AddElement(copperNist, fractionmass = 70 *perCent);

 
  // MATERIAL ASSIGNMENT

  // Material of the FINAL APPLICATOR FLASH
  firstApplicatorMaterialFLASH = PMMANist;

// Material of the FINAL APPLICATOR FLASH
  finalApplicatorMaterialFLASH = PMMANist;

 // monitor 1 FINAL APPLICATOR FLASH
  Giunz1FinalAppMaterialFLASH = aluminumNist;

 // monitor 2 FINAL APPLICATOR FLASH
  Giunz2FinalAppMaterialFLASH = aluminumNist;
 
 // part of monitor
 cover1material = aluminumNist;
 
 //part of monitor
 cover2material = aluminumNist;
 
 //part of applicator
 cover3material = PMMANist;
 
 //all protectors encase the monitors and are made of pmma
 protector1material=PMMANist;
 protector2material=PMMANist;
 protector3material=PMMANist;
 protector4material=PMMANist;
 
 // monitor 3 FINAL APPLICATOR Int FLASH
  Giunz3FinalAppMaterialFLASH = aluminumNist;

 // monitor 4 FINAL APPLICATOR FLASH
  Giunz4FinalAppMaterialFLASH = aluminumNist;

 // Junction 5 FINAL APPLICATOR FLASH
  Giunz5FinalAppMaterialFLASH = PMMANist;
  
 // Superior Final Part Monitor Chambers Material
    FTFLASHMaterialFLASH = titanioNist;
	
 // Vacuum Source
    VSFLASHMaterialFLASH = galacticNist;

}

/////////////////////////////////////////////////////////////////////////////
void Applicator80BeamLine::ConstructApplicator80BeamLine()
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
 
  // Components of the Applicator Beam Line

  FlashBeamLineVacuumSource();
  FlashBeamLineTitaniumWindows();
  FlashBeamLineFirstApplicator();
  
  FlashBeamLineJunctions(); 
  FlashBeamLineFinalApplicator();
    
}


void Applicator80BeamLine::FlashBeamLineVacuumSource()
{
 // ---------------------------------------------------------------//
  //                     Vacuum Source                             //
  // ---------------------------------------------------------------//

  
  G4double phi1 = 90. *deg;     

            
   G4RotationMatrix rm1;               
   rm1.rotateY(phi1);

  const G4double outRadiusVSFLASH = 44.75 *mm;
  const G4double innRadiusVSFLASH = 0.*mm;
  const G4double hightVSFLASH = 0.5 *mm;
  const G4double startAngleVSFLASH = 0.*deg;
  const G4double spanningAngleVSFLASH = 360.*deg;
  const G4double XPositionVSFLASH = -1000.0 *mm;
    
  solidVSFLASH = new G4Tubs("VSFLASH", innRadiusVSFLASH, 
				    outRadiusVSFLASH,
				    hightVSFLASH, 
				    startAngleVSFLASH, 
				    spanningAngleVSFLASH);

  G4LogicalVolume* logVSFLASH = new G4LogicalVolume(solidVSFLASH, 
							      VSFLASHMaterialFLASH, "VSFLASH", 0, 0, 0);

  physiVSFLASH = new G4PVPlacement(G4Transform3D(rm1, G4ThreeVector((XPositionVSFLASH),0.,0.)),
					   "VSFLASH", logVSFLASH, physicalTreatmentRoom, false, 0); 

  logVSFLASH -> SetVisAttributes(green);
}
/////////
void Applicator80BeamLine::FlashBeamLineTitaniumWindows()
{
// ---------------------------------------------------------------//
  //                     Titanium Window                        //
  // ---------------------------------------------------------------//
	//with just this piece ssd=1.6cm
  G4double phi2 = 90. *deg;     

            
   G4RotationMatrix rm2;               
   rm2.rotateY(phi2);


  const G4double outRadiusFTFLASH = 44.75 *mm;
  const G4double innRadiusFTFLASH = 8.5 *mm;
  const G4double hightFTFLASH = 0.03 *mm;
  const G4double startAngleFTFLASH = 0.*deg;
  const G4double spanningAngleFTFLASH = 360.*deg;
  const G4double XPositionFTFLASH = -999.47 *mm;

  solidFTFLASH = new G4Tubs("FTFLASH", innRadiusFTFLASH, 
				    outRadiusFTFLASH,
				    hightFTFLASH, 
				    startAngleFTFLASH, 
				    spanningAngleFTFLASH);

  G4LogicalVolume* logFTFLASH = new G4LogicalVolume(solidFTFLASH, 
							      FTFLASHMaterialFLASH, "FTFLASH", 0, 0, 0);

  physiFTFLASH = new G4PVPlacement(G4Transform3D(rm2, G4ThreeVector((XPositionFTFLASH),0.,0.)),
					   "FTFLASH", logFTFLASH, physicalTreatmentRoom, false, 0); 

  logFTFLASH -> SetVisAttributes(yellow);
}
void Applicator80BeamLine::FlashBeamLineFirstApplicator()
{


	// -----------------------//
	  // FIRST APPLICATOR FLASH  //
	  //------------------------//

	 // const G4double outRadiusFirstApplicatorFLASH = 45. *mm;
	 // const G4double innRadiusFirstApplicatorFLASH = 40. *mm;
	const G4double hightFirstApplicatorFLASH = 100. * mm;
	const G4double startAngleFirstApplicatorFLASH = 0. * deg;
	const G4double spanningAngleFirstApplicatorFLASH = 360. * deg;
	const G4double firstApplicatorXPositionFLASH = -799.44 * mm;




	G4double phi6 = 90. * deg;


	G4RotationMatrix rm6;
	rm6.rotateY(phi6);


	solidFirstApplicatorFLASH = new G4Tubs("FirstApplicatorFLASH", innerRadiusFirstApplicatorFLASH,
		OuterRadiusFirstApplicatorFLASH,
		hightFirstApplicatorFLASH,
		startAngleFirstApplicatorFLASH,
		spanningAngleFirstApplicatorFLASH);

	G4LogicalVolume* logFirstApplicatorFLASH = new G4LogicalVolume(solidFirstApplicatorFLASH,
		firstApplicatorMaterialFLASH, "FirstApplicatorFLASH", 0, 0, 0);

	physiFirstApplicatorFLASH = new G4PVPlacement(G4Transform3D(rm6, G4ThreeVector((firstApplicatorXPositionFLASH), 0., 0.)),
		"FirstApplicatorFLASH", logFirstApplicatorFLASH, physicalTreatmentRoom, false, 0);

	//  logFirstApplicatorFLASH -> SetVisAttributes(G4VisAttributes::GetInvisible()); 
	logFirstApplicatorFLASH->SetVisAttributes(magenta);
}


void Applicator80BeamLine::FlashBeamLineJunctions()
{


	G4double phi5 = 90. * deg;


	G4RotationMatrix rm5;
	rm5.rotateY(phi5);
	// --------------------------------- //
	  // Junction 5 FIRST APPLICATOR FLASH //
	  // --------------------------------- //

	const G4double outRadiusGiunz5FinalAppFLASH = 60 * mm;
	const G4double innRadiusGiunz5FinalAppFLASH = 55 * mm;
	const G4double hightGiunz5FinalAppFLASH = 10. * mm;
	const G4double startAngleGiunz5FinalAppFLASH = 0. * deg;
	const G4double spanningAngleGiunz5FinalAppFLASH = 360. * deg;
	const G4double Giunz5FinalAppXPositionFLASH = -689.44 * mm;

	solidGiunz5FinalAppFLASH = new G4Tubs("Giunz5FinalAppFLASH", innRadiusGiunz5FinalAppFLASH,
		outRadiusGiunz5FinalAppFLASH,
		hightGiunz5FinalAppFLASH,
		startAngleGiunz5FinalAppFLASH,
		spanningAngleGiunz5FinalAppFLASH);

	G4LogicalVolume* logGiunz5FinalAppFLASH = new G4LogicalVolume(solidGiunz5FinalAppFLASH,
		Giunz5FinalAppMaterialFLASH, "Giunz5FinalAppFLASH", 0, 0, 0);

	physiGiunz5FinalAppFLASH = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((Giunz5FinalAppXPositionFLASH), 0., 0.)),
		"Giunz5FinalAppFLASH", logGiunz5FinalAppFLASH, physicalTreatmentRoom, false, 0);

	logGiunz5FinalAppFLASH->SetVisAttributes(yellow);

	// --------------------------------- //
	  // Junction 4 FINAL APPLICATOR FLASH //
	  // --------------------------------- //

	const G4double innRadiusGiunz4FinalAppFLASH = 22.25 * mm;
	const G4double outRadiusGiunz4FinalAppFLASH = 27.25 * mm;
	const G4double hightGiunz4FinalAppFLASH = 9.35 * mm;
	const G4double startAngleGiunz4FinalAppFLASH = 0. * deg;
	const G4double spanningAngleGiunz4FinalAppFLASH = 360. * deg;
	const G4double Giunz4FinalAppXPositionFLASH = -908.79 * mm;

	solidGiunz4FinalAppFLASH = new G4Tubs("Giunz4FinalAppFLASH", innRadiusGiunz4FinalAppFLASH,
		outRadiusGiunz4FinalAppFLASH,
		hightGiunz4FinalAppFLASH,
		startAngleGiunz4FinalAppFLASH,
		spanningAngleGiunz4FinalAppFLASH);

	G4LogicalVolume* logGiunz4FinalAppFLASH = new G4LogicalVolume(solidGiunz4FinalAppFLASH,
		Giunz4FinalAppMaterialFLASH, "Giunz4FinalAppFLASH", 0, 0, 0);

	physiGiunz4FinalAppFLASH = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((Giunz4FinalAppXPositionFLASH), 0., 0.)),
		"Giunz4FinalAppFLASH", logGiunz4FinalAppFLASH, physicalTreatmentRoom, false, 0);

	logGiunz4FinalAppFLASH->SetVisAttributes(blue);
//protector4
const G4double outRadiusprotector4FLASH = 32 * mm;
	const G4double innRadiusprotector4FLASH = 27.25 * mm;
	const G4double hightprotector4FLASH = 9.35 * mm;
	const G4double startAngleprotector4FLASH = 0. * deg;
	const G4double spanningAngleprotector4FLASH = 360. * deg;
	const G4double protector4XPositionFLASH = -908.79 * mm;

	protector4 = new G4Tubs("protector4", innRadiusprotector4FLASH,
		outRadiusprotector4FLASH,
		hightprotector4FLASH,
		startAngleprotector4FLASH,
		spanningAngleprotector4FLASH);

	G4LogicalVolume* logprotector4FLASH = new G4LogicalVolume(protector4,
		protector4material, "protector4", 0, 0, 0);

	physiprotector4 = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((protector4XPositionFLASH), 0., 0.)),
		"protector4", logprotector4FLASH, physicalTreatmentRoom, false, 0);



	// --------------------------------- //
	 // monitor 3 FINAL APPLICATOR FLASH //
	 // --------------------------------- //

	const G4double outRadiusGiunz3FinalAppFLASH = 18.25 * mm;
	const G4double innRadiusGiunz3FinalAppFLASH =  13.75 * mm;
	const G4double hightGiunz3FinalAppFLASH = 4.65 * mm;
	const G4double startAngleGiunz3FinalAppFLASH = 0. * deg;
	const G4double spanningAngleGiunz3FinalAppFLASH = 360. * deg;
	const G4double Giunz3FinalAppXPositionFLASH = -922.79 * mm;

	solidGiunz3FinalAppFLASH = new G4Tubs("Giunz3FinalAppFLASH", innRadiusGiunz3FinalAppFLASH,
		outRadiusGiunz3FinalAppFLASH,
		hightGiunz3FinalAppFLASH,
		startAngleGiunz3FinalAppFLASH,
		spanningAngleGiunz3FinalAppFLASH);

	G4LogicalVolume* logicsolidGiunz3FinalAppFLASH = new G4LogicalVolume(solidGiunz3FinalAppFLASH,
		Giunz3FinalAppMaterialFLASH, "Giunz3FinalAppFLASH", 0, 0, 0);

	physiGiunz3FinalAppFLASH = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((Giunz3FinalAppXPositionFLASH), 0., 0.)),
		"Giunz3FinalAppFLASH", logicsolidGiunz3FinalAppFLASH, physicalTreatmentRoom, false, 0);

	logicsolidGiunz3FinalAppFLASH->SetVisAttributes(yellow);
	
//protector3
const G4double outRadiusprotector3FLASH = 25 * mm;
	const G4double innRadiusprotector3FLASH = 18.25 * mm;
	const G4double hightprotector3FLASH = 4.65 * mm;
	const G4double startAngleprotector3FLASH = 0. * deg;
	const G4double spanningAngleprotector3FLASH = 360. * deg;
	const G4double protector3XPositionFLASH = -922.79 * mm;

	protector3 = new G4Tubs("protector3", innRadiusprotector3FLASH,
		outRadiusprotector3FLASH,
		hightprotector3FLASH,
		startAngleprotector3FLASH,
		spanningAngleprotector3FLASH);

	G4LogicalVolume* logprotector3FLASH = new G4LogicalVolume(protector3,
		protector3material, "protector3", 0, 0, 0);

	physiprotector3 = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((protector3XPositionFLASH), 0., 0.)),
		"protector3", logprotector3FLASH, physicalTreatmentRoom, false, 0);

	logprotector3FLASH->SetVisAttributes(gray);
	//cover disk between monitor piece 3 and 4

const G4double innRadiuscover2FLASH = 18.25 * mm;
	const G4double outRadiuscover2FLASH = 22.25 * mm;
	const G4double hightcover2FLASH = 0.01 * mm;
	const G4double startAnglecover2FLASH= 0. * deg;
	const G4double spanningAnglecover2FLASH = 360. * deg;
	const G4double cover2XPositionFLASH = -918.14 * mm;

	cover2 = new G4Tubs("cover2", innRadiuscover2FLASH,
		outRadiuscover2FLASH,
		hightcover2FLASH,
		startAnglecover2FLASH,
		spanningAnglecover2FLASH);

	G4LogicalVolume* logcover2 = new G4LogicalVolume(cover2,
		cover2material, "cover2", 0, 0, 0);

	physicover2FLASH = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((cover2XPositionFLASH), 0., 0.)),
		"cover2", logcover2, physicalTreatmentRoom, false, 0);

	logcover2->SetVisAttributes(gray);
	// --------------------------------- //
	  // monitor 2 FINAL APPLICATOR FLASH //
	  // --------------------------------- //

	const G4double innRadiusGiunz2FinalAppFLASH = 15.0 * mm;
	const G4double outRadiusGiunz2FinalAppFLASH = 20.0 * mm;
	const G4double hightGiunz2FinalAppFLASH = 30 * mm;
	const G4double startAngleGiunz2FinalAppFLASH = 0. * deg;
	const G4double spanningAngleGiunz2FinalAppFLASH = 360. * deg;
	const G4double Giunz2FinalAppXPositionFLASH = -957.44 * mm;

	solidGiunz2FinalAppFLASH = new G4Tubs("Giunz2FinalAppFLASH", innRadiusGiunz2FinalAppFLASH,
		outRadiusGiunz2FinalAppFLASH,
		hightGiunz2FinalAppFLASH,
		startAngleGiunz2FinalAppFLASH,
		spanningAngleGiunz2FinalAppFLASH);

	G4LogicalVolume* logGiunz2FinalAppFLASH = new G4LogicalVolume(solidGiunz2FinalAppFLASH,
		Giunz2FinalAppMaterialFLASH, "Giunz2FinalAppFLASH", 0, 0, 0);

	physiGiunz2FinalAppFLASH = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((Giunz2FinalAppXPositionFLASH), 0., 0.)),
		"Giunz2FinalAppFLASH", logGiunz2FinalAppFLASH, physicalTreatmentRoom, false, 0);

	logGiunz2FinalAppFLASH->SetVisAttributes(white);

//protector2
const G4double outRadiusprotector2FLASH = 25 * mm;
	const G4double innRadiusprotector2FLASH = 20.0 * mm;
	const G4double hightprotector2FLASH = 30 * mm;
	const G4double startAngleprotector2FLASH = 0. * deg;
	const G4double spanningAngleprotector2FLASH = 360. * deg;
	const G4double protector2XPositionFLASH = -957.44 * mm;

	protector2 = new G4Tubs("protector2", innRadiusprotector2FLASH,
		outRadiusprotector2FLASH,
		hightprotector2FLASH,
		startAngleprotector2FLASH,
		spanningAngleprotector2FLASH);

	G4LogicalVolume* logprotector2FLASH = new G4LogicalVolume(protector2,
		protector2material, "protector2", 0, 0, 0);

	physiprotector2 = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((protector2XPositionFLASH), 0., 0.)),
		"protector2", logprotector2FLASH, physicalTreatmentRoom, false, 0);

	logprotector2FLASH->SetVisAttributes(gray);

	// --------------------------------- //
	  // monitor 1 FINAL APPLICATOR FLASH //
	  // --------------------------------- //

	const G4double innRadiusGiunz1FinalAppFLASH = 10.0 * mm;
	const G4double outRadiusGiunz1FinalAppFLASH = 15.0 * mm;
	const G4double hightGiunz1FinalAppFLASH = 6.0 * mm;
	const G4double startAngleGiunz1FinalAppFLASH = 0. * deg;
	const G4double spanningAngleGiunz1FinalAppFLASH = 360. * deg;
	const G4double Giunz1FinalAppXPositionFLASH = -993.44 * mm;

	solidGiunz1FinalAppFLASH = new G4Tubs("Giunz1FinalAppFLASH", innRadiusGiunz1FinalAppFLASH,
		outRadiusGiunz1FinalAppFLASH,
		hightGiunz1FinalAppFLASH,
		startAngleGiunz1FinalAppFLASH,
		spanningAngleGiunz1FinalAppFLASH);

	G4LogicalVolume* logGiunz1FinalAppFLASH = new G4LogicalVolume(solidGiunz1FinalAppFLASH,
		Giunz1FinalAppMaterialFLASH, "Giunz1FinalAppFLASH", 0, 0, 0);

	physiGiunz1FinalAppFLASH = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((Giunz1FinalAppXPositionFLASH), 0., 0.)),
		"Giunz1FinalAppFLASH", logGiunz1FinalAppFLASH, physicalTreatmentRoom, false, 0);

	logGiunz1FinalAppFLASH->SetVisAttributes(gray);
	
//protector1
const G4double outRadiusprotector1FLASH = 20 * mm;
	const G4double innRadiusprotector1FLASH = 15 * mm;
	const G4double hightprotector1FLASH = 6.0 * mm;
	const G4double startAngleprotector1FLASH = 0. * deg;
	const G4double spanningAngleprotector1FLASH = 360. * deg;
	const G4double protector1XPositionFLASH = -993.44 * mm;

	protector1 = new G4Tubs("protector1", innRadiusprotector1FLASH,
		outRadiusprotector1FLASH,
		hightprotector1FLASH,
		startAngleprotector1FLASH,
		spanningAngleprotector1FLASH);

	G4LogicalVolume* logprotector1FLASH = new G4LogicalVolume(protector1,
		protector1material, "protector1", 0, 0, 0);

	physiprotector1 = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((protector1XPositionFLASH), 0., 0.)),
		"protector1", logprotector1FLASH, physicalTreatmentRoom, false, 0);

	logprotector1FLASH->SetVisAttributes(gray);
	
	
//cover disk between monitor piece 1 and 2

const G4double outRadiuscover1FLASH = 20. * mm;
	const G4double innRadiuscover1FLASH = 15 * mm;
	const G4double hightcover1FLASH = 0.01 * mm;
	const G4double startAnglecover1FLASH= 0. * deg;
	const G4double spanningAnglecover1FLASH = 360. * deg;
	const G4double cover1XPositionFLASH = -987.44 * mm;

	cover1 = new G4Tubs("cover1", innRadiuscover1FLASH,
		outRadiuscover1FLASH,
		hightcover1FLASH,
		startAnglecover1FLASH,
		spanningAnglecover1FLASH);

	G4LogicalVolume* logcover1 = new G4LogicalVolume(cover1,
		cover1material, "cover1", 0, 0, 0);

	physicover1FLASH = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((cover1XPositionFLASH), 0., 0.)),
		"cover1", logcover1, physicalTreatmentRoom, false, 0);

	logcover1->SetVisAttributes(gray);
	
	//cover disk between monitor piece 4 and first applicator

const G4double outRadiuscover3FLASH = 50.0 * mm;
	const G4double innRadiuscover3FLASH = 24.75 * mm;
	const G4double hightcover3FLASH = 0.01 * mm;
	const G4double startAnglecover3FLASH= 0. * deg;
	const G4double spanningAnglecover3FLASH = 360. * deg;
	const G4double cover3XPositionFLASH = -899.44 * mm;

	cover3 = new G4Tubs("cover3", innRadiuscover3FLASH,
		outRadiuscover3FLASH,
		hightcover3FLASH,
		startAnglecover3FLASH,
		spanningAnglecover3FLASH);

	G4LogicalVolume* logcover3 = new G4LogicalVolume(cover3,
		cover3material, "cover3", 0, 0, 0);

	physicover3FLASH = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((cover3XPositionFLASH), 0., 0.)),
		"cover3", logcover3, physicalTreatmentRoom, false, 0);

	logcover3->SetVisAttributes(gray);
}


void Applicator80BeamLine::FlashBeamLineFinalApplicator()
{
// -----------------------//
  // FINAL APPLICATOR FLASH  //
  //------------------------//

 // const G4double outRadiusFinalApplicatorFLASH = 45. *mm;
 // const G4double innRadiusFinalApplicatorFLASH = 40. *mm;
  const G4double hightFinalApplicatorFLASH = 250.0 *mm;
  const G4double startAngleFinalApplicatorFLASH = 0.*deg;
  const G4double spanningAngleFinalApplicatorFLASH = 360.*deg;
  const G4double finalApplicatorXPositionFLASH = -449.44*mm;

  
 
  
  G4double phi6 = 90. *deg;     

         
   G4RotationMatrix rm6;               
   rm6.rotateY(phi6);

    
  solidFinalApplicatorFLASH = new G4Tubs("FinalApplicatorFLASH", innerRadiusFinalApplicatorFLASH, 
				    OuterRadiusFinalApplicatorFLASH,
				    hightFinalApplicatorFLASH, 
				    startAngleFinalApplicatorFLASH, 
				    spanningAngleFinalApplicatorFLASH);

  G4LogicalVolume* logFinalApplicatorFLASH = new G4LogicalVolume(solidFinalApplicatorFLASH, 
							      finalApplicatorMaterialFLASH, "FinalApplicatorFLASH", 0, 0, 0);

  physiFinalApplicatorFLASH = new G4PVPlacement(G4Transform3D(rm6, G4ThreeVector((finalApplicatorXPositionFLASH),0.,0.)),
					   "FinalApplicatorFLASH", logFinalApplicatorFLASH, physicalTreatmentRoom, false, 0); 

  //  logFinalApplicatorFLASH -> SetVisAttributes(G4VisAttributes::GetInvisible()); 
  logFinalApplicatorFLASH -> SetVisAttributes(magenta); 
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////// MESSENGER ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


void Applicator80BeamLine::SetInnerRadiusFinalApplicatorFLASH(G4double value)
{///NB here it isn't enough to set the radius of applicator, but also of all the junctions
  solidFinalApplicatorFLASH -> SetInnerRadius(value);
  solidFirstApplicatorFLASH->SetInnerRadius(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout<<"Inner Radius of the applicator FLASH is (mm):"
	<< solidFinalApplicatorFLASH -> GetInnerRadius()/mm
	<< G4endl; 
}

/////////////////////////////////////////////////////////////////////////

void Applicator80BeamLine::SetOuterRadiusFinalApplicatorFLASH(G4double value)
{
  solidFinalApplicatorFLASH -> SetOuterRadius(value);
  solidFirstApplicatorFLASH -> SetOuterRadius(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout<<"Outer Radius of the final applicator FLASH is (mm):"
	<< solidFinalApplicatorFLASH -> GetOuterRadius()/mm
	<< G4endl; 
}

/////////////////////////////////////////////////////////////////////////////



