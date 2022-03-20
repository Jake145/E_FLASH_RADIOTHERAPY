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
#include "doiPETGlobalParameters.hh"
#include "doiPETAnalysis.hh"
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


#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashDetectorConstruction::FlashDetectorConstruction()
    : G4VUserDetectorConstruction(), physicalTreatmentRoom(0),
      fCheckOverlaps(true) {

  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashDetectorConstruction::~FlashDetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FlashDetectorConstruction::DefineMaterials() {


 
 

  G4NistManager *man = G4NistManager::Instance();
  //Defile Aluminum material for the detetor cover
	Aluminum = man->FindOrBuildMaterial("G4_Al", isotopes);

	//Define elements for the GSO  crystal (scintillator) 
	O = man->FindOrBuildElement("O" , isotopes); 
	Si = man->FindOrBuildElement("Si", isotopes);
	Gd = man->FindOrBuildElement("Gd", isotopes);  


	//define GSO crystal for PET detector
	GSO = new G4Material("GSO", 6.7*g/cm3, 3);
	GSO->AddElement(Gd, 2);
	GSO->AddElement(Si, 1);
	GSO->AddElement(O,  5); 
	crystalMaterial   = man->FindOrBuildMaterial("GSO");
  G4int ncomponents;
  G4bool isotopes = false;
  G4double prelude_density = 7.25 * g / cm3;
  G4Material *prelude =
      new G4Material("prelude", prelude_density, ncomponents = 4);
  prelude->AddElement(man->FindOrBuildElement("Lu"), 71 * perCent);
  prelude->AddElement(man->FindOrBuildElement("Si"), 7 * perCent);
  prelude->AddElement(man->FindOrBuildElement("O"), 18 * perCent);
  prelude->AddElement(man->FindOrBuildElement("Y"), 4 * perCent);



  G4Material* LSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
  LSO->AddElement(man->FindOrBuildElement("Lu"), 2);
  LSO->AddElement(man->FindOrBuildElement("Si"), 1);
  LSO->AddElement(man->FindOrBuildElement("O") , 5);  
  
  
  G4Material *scintillator =
      new G4Material("scintillator", prelude_density, ncomponents = 2);
  scintillator->AddMaterial(prelude, 99.81 * perCent);
  scintillator->AddElement(man->FindOrBuildElement("Ce"), 0.19 * perCent);

  

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
 

  std::vector<G4int> natoms;

  std::vector<G4String> elements;

  nist = G4NistManager::Instance();
  cryst_mat = nist->FindOrBuildMaterial("scintillator");

  elements.push_back("C");
  natoms.push_back(5);
  elements.push_back("H");
  natoms.push_back(8);
  elements.push_back("O");
  natoms.push_back(2);

  density = 1.190 * g / cm3;

  PMMA = nist->ConstructNewMaterial("PMMA", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Cladding (polyethylene)
  //--------------------------------------------------

  elements.push_back("C");
  natoms.push_back(2);
  elements.push_back("H");
  natoms.push_back(4);

  density = 1.200 * g / cm3;

  PE = nist->ConstructNewMaterial("Pethylene", elements, natoms, density);

  elements.clear();
  natoms.clear();
 

  TEFLON = nist->FindOrBuildMaterial("G4_TEFLON", isotopes);

  
  ej212 = new G4Material("EJ-212", density = 1.023 * g / cm3, ncomponents = 2);
  ej212->AddElement(fH, 0.5243407708);
  ej212->AddElement(fC, 0.4756592292);

  
  EJ212 = nist->FindOrBuildMaterial("EJ-212");
}

void FlashDetectorConstruction::Construct_PET()
{  
  
	//Each block detector is identified with its unique number, blockIndex.
	blockIndex = 0;

	//Each crystal is identified with its unique number, crystalIndex
	crystalIndex = 0;



	//Define air volume (box) to fill the detector block. Crystal elements (scintillators) is then placed.
	sizeOfAirBox_DOI = (numberOfCrystal_DOI * sizeOfCrystal_DOI) + (numberOfCrystal_DOI - 1)*crystalGap_DOI;
	sizeOfAirBox_axial = (numberOfCrystal_axial * sizeOfCrystal_axial) + (numberOfCrystal_axial - 1)*crystalGap_axial;
	sizeOfAirBox_tangential = (numberOfCrystal_tangential * sizeOfCrystal_tangential) + (numberOfCrystal_tangential - 1)*crystalGap_tangential;

	
	pAnalysis = doiPETAnalysis::GetInstance();
	pAnalysis->GetSizeOfDetector(sizeOfAirBox_DOI,sizeOfAirBox_tangential, sizeOfAirBox_axial);

	G4cout<<"size of crytal element: "<<sizeOfCrystal_tangential<<" "<<sizeOfCrystal_axial<<" "<<sizeOfCrystal_DOI<<G4endl;
	G4cout<<"Size of detector block (without Al cover): "<<sizeOfAirBox_tangential<<" "<<sizeOfAirBox_axial<<" "<<sizeOfAirBox_DOI<<G4endl;



	//Define the size of the detector block. 
	sizeOfBlockDetector_DOI = sizeOfAirBox_DOI + AluminumCoverThickness;
	sizeOfBlockDetector_axial = sizeOfAirBox_axial + AluminumCoverThickness;
	sizeOfBlockDetector_tangential = sizeOfAirBox_tangential + AluminumCoverThickness;



	//Define solid shape for the detector block
	G4Box* blockDetector = new G4Box("blockDetector",sizeOfBlockDetector_DOI/2,sizeOfBlockDetector_tangential/2,sizeOfBlockDetector_axial/2);

	//Define the logical volume for the detector block
	blockDetector_logicalV = new G4LogicalVolume(blockDetector,Aluminum,"blockDetector_logicalV", 0,0,0);



	//Define air (box) inside the detector block. Crystal elements will be placed in it.
	G4Box* airBox = new G4Box("airBox", sizeOfAirBox_DOI/2, sizeOfAirBox_tangential/2,sizeOfAirBox_axial/2);

	//Define the logical volume
	airBox_logicalV = new G4LogicalVolume(airBox,airNist,"airBox_logicalV", 0,0,0);

	//Define its physical volume and place it inside the detector block
	airBox_physicalV = new G4PVPlacement (0,G4ThreeVector(0,0,0),airBox_logicalV,"airBox_physicalV", blockDetector_logicalV,false,0,fCheckOverlaps);


	///////////////////////////////////////// Arrange the PET ring and place the PET detectors in the ring(s) ////////////////////////////////////

	for(G4int Ring = 0; Ring< numberOfRings; Ring++)
	{
		//place the detectors in a ring along the axial direction. Note that the ring gap between two adjcent rings is measured from scintillator to scintillator.  It does not include the Aluminum thickness cover.

		detectorPositionZ = 30*cm+(Ring-((G4double)numberOfRings)/2 + 0.5)*(sizeOfBlockDetector_axial + ringGap - AluminumCoverThickness);

		for(G4int i = 0; i<numberOfDetector_perRing; i++)
		{
			//The azimuthal angle to arrange the detectors in a ring
			thetaDetector = (double)(i*twopi/numberOfDetector_perRing);

			//The radius of the scanner is measured from opposing crystal (scintillator) faces. It does not include the Aluminum thickness cover.
			detectorPositionX = (scannerRadius + sizeOfBlockDetector_DOI/2 - AluminumCoverThickness/2)*std::cos(thetaDetector);
			detectorPositionY = (scannerRadius + sizeOfBlockDetector_DOI/2 - AluminumCoverThickness/2)*std::sin(thetaDetector);

			//Define the rotation matrix for correct placement of detetors
			G4RotationMatrix rotm_PET = G4RotationMatrix();
			rotm_PET.rotateZ(thetaDetector);
			G4ThreeVector uz_PET = G4ThreeVector(detectorPositionX,detectorPositionY,detectorPositionZ);
			G4Transform3D transform = G4Transform3D(rotm_PET,uz_PET);

			//Define the physical volume of the detectors.
			blockDetector_physicalV = new G4PVPlacement (transform,blockDetector_logicalV,"blockDetector_physicalV", logicTreatmentRoom,false,blockIndex,fCheckOverlaps);
			blockIndex++;
			//G4cout<<Ring<<" "<<detectorPositionX- ((sizeOfBlockDetector_DOI - AluminumCoverThickness)/2)*cos(thetaDetector)<<" "<<detectorPositionY- ((sizeOfBlockDetector_DOI- AluminumCoverThickness)/2)*sin(thetaDetector)<<" "<<detectorPositionZ<<G4endl;

		}
	}

	//Define the solid crystal
	G4VSolid* CrystalSolid = new G4Box("Crystal", sizeOfCrystal_DOI/2., sizeOfCrystal_tangential/2., sizeOfCrystal_axial/2.);

	//Define the local volume of the crystal
	crystal_logicalV = new G4LogicalVolume(CrystalSolid,crystalMaterial,"Crystal_logicalV", 0,0,0);

	//Place the crystals inside the detectors and give them a unique number with crystalIndex.
	for(G4int i_DOI = 0; i_DOI<numberOfCrystal_DOI; i_DOI++){
		crystalPositionX=(i_DOI-((G4double)numberOfCrystal_DOI)/2 + 0.5)*(sizeOfCrystal_DOI + crystalGap_DOI);
		for(G4int i_axial=0; i_axial< numberOfCrystal_axial;i_axial++){
			crystalPositionZ = (i_axial-((G4double)numberOfCrystal_axial)/2 + 0.5)*(sizeOfCrystal_axial + crystalGap_axial);
			for(G4int i_tan=0; i_tan<numberOfCrystal_tangential;i_tan++){
				crystalPositionY=(i_tan-((G4double)numberOfCrystal_tangential)/2 + 0.5)*(sizeOfCrystal_tangential + crystalGap_tangential);

				//G4cout<<crystalIndex<<" "<<crystalPositionX<<" "<<crystalPositionY<<" "<<crystalPositionZ<<G4endl;
				//place the crystal inside the block detector. 
				crystal_physicalV = new G4PVPlacement (0, G4ThreeVector (crystalPositionX,crystalPositionY,crystalPositionZ), crystal_logicalV, "Crystal_physicalV", airBox_logicalV,false,crystalIndex/*,fCheckOverlaps*/);
				crystalIndex++;
			}
		}
	}

	//******************  Visualization *****************************//

	//visualization for the block detector
	G4VisAttributes* blockDetectorVisAtt;
	blockDetectorVisAtt = new G4VisAttributes(G4Colour(1,1.0,1.0));
	blockDetectorVisAtt->SetVisibility (true);
	//blockDetectorVisAtt->SetForceWireframe (true);
	blockDetector_logicalV->SetVisAttributes (blockDetectorVisAtt);
	//blockDetector_logicalV->SetVisAttributes (G4VisAttributes::Invisible);

	//visualization for the the box filled with air
	G4VisAttributes* airBoxVisAtt;
	airBoxVisAtt = new G4VisAttributes(G4Colour(1,1.0,1.0));
	airBoxVisAtt->SetVisibility (true);
	airBoxVisAtt->SetForceWireframe (true);
	airBox_logicalV->SetVisAttributes (airBoxVisAtt);
	airBox_logicalV->SetVisAttributes (G4VisAttributes::Invisible);

	//visualization for the crystal
	G4VisAttributes* crystalVisAtt;
	crystalVisAtt = new G4VisAttributes(G4Colour(0.88,0.55,1.0));
	//crystalVisAtt->SetVisibility (true);
	crystalVisAtt->SetForceWireframe (true);
	crystal_logicalV->SetVisAttributes (crystalVisAtt);
	crystal_logicalV->SetVisAttributes (G4VisAttributes::Invisible);
 


}

G4VPhysicalVolume *
FlashDetectorConstruction::ConstructPhantom(G4double CollPos, G4bool dishomo) {

G4bool fCheckOverlaps=true;
red = new G4VisAttributes(G4Colour(0 / 255., 255 / 255., 0 / 255.));
  red->SetVisibility(true);

  blue = new G4VisAttributes(G4Colour(0 / 255., 0. / 255., 255. / 255.));
  blue->SetVisibility(true);

  G4Material *phantomMaterial = nist->FindOrBuildMaterial("G4_WATER");

  
    G4double Position_coefficient = CollPos;
if (dishomo == false){

  G4double phantomSizeX = 20.0 * cm, phantomSizeY = 20.0 * cm,
           phantomSizeZ = 60.0 * cm,
           phantom_coordinateZ = (Position_coefficient * mm + phantomSizeZ / 2);

  G4ThreeVector phantomPosition =
      G4ThreeVector(0*mm, 0. * mm, phantom_coordinateZ);
  // Definition of the solid volume of the Phantom
  phantom = new G4Box("Phantom", phantomSizeX / 2, phantomSizeY / 2,
                      phantomSizeZ / 2);

  // Definition of the logical volume of the Phantom
  phantomLogicalVolume =
      new G4LogicalVolume(phantom, PMMA, "phantomLog", 0, 0, 0);

  // Definition of the physics volume of the Phantom
  phant_phys =
      new G4PVPlacement(0, phantomPosition, "phantomPhys", phantomLogicalVolume,
                        physicalTreatmentRoom, false, 0,fCheckOverlaps);

  G4Region *PhantomRegion = new G4Region("Phantom_reg");
  phantomLogicalVolume->SetRegion(PhantomRegion);
  PhantomRegion->AddRootLogicalVolume(phantomLogicalVolume);

  // Visualisation attributes of the phantom
  red = new G4VisAttributes(G4Colour(0 / 255., 255 / 255., 0 / 255.));
  red->SetVisibility(true);

  blue = new G4VisAttributes(G4Colour(0 / 255., 0. / 255., 255. / 255.));
  blue->SetVisibility(true);

  phantomLogicalVolume->SetVisAttributes(red);

  G4double maxStep = 0.1 * mm;
  fStepLimit = new G4UserLimits(maxStep);
  phantomLogicalVolume->SetUserLimits(fStepLimit);}
  
  else{
  
  G4Material* bone = G4NistManager::Instance()->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU", false);
  //=====================================first piece of phantom====================//
  G4double Depth_firstpiece = 40 * cm;
  G4double phantomSizeZ = Depth_firstpiece,
           phantomSizeY = 20.0 * cm, phantomSizeX = 20.0 * cm,
           phantom_coordinateZ = (Position_coefficient * mm + phantomSizeX / 2);

  if (phantomSizeX != 0) {
    G4ThreeVector phantomPosition =
        G4ThreeVector( 0. * mm, 0. * mm,phantom_coordinateZ);
    // Definition of the solid volume of the Phantom
    phantom = new G4Box("Phantom", phantomSizeX / 2, phantomSizeY / 2,
                        phantomSizeZ / 2);

    // Definition of the logical volume of the Phantom
    phantomLogicalVolume =
        new G4LogicalVolume(phantom, phantomMaterial, "phantomLog", 0, 0, 0);

    // Definition of the physics volume of the Phantom
    phant_phys = new G4PVPlacement(0, phantomPosition, "phantomPhys",
                                   phantomLogicalVolume, physicalTreatmentRoom,
                                   false, 0,fCheckOverlaps);
  }

  //============================Dishomogeneity=====================================//

 

    support_z = 10 * cm;
    G4Box *support_whole =
        new G4Box("Support_w", phantomSizeX / 2, phantomSizeY / 2, support_z / 2);

    DetectorSupport = new G4LogicalVolume(support_whole, airNist, "SupportLog");

    G4RotationMatrix rotmp = G4RotationMatrix();
    rotmp.rotateY(0 * deg);
    supp_coordinateZ = phantom_coordinateZ + (phantomSizeZ / 2 + support_z / 2);

    G4ThreeVector positionp = G4ThreeVector( 0, 0,phantomSizeY);
    G4Transform3D transformp = G4Transform3D(rotmp, positionp);

    new G4PVPlacement(transformp, DetectorSupport, "supportphys",
                      logicTreatmentRoom, false, 0,fCheckOverlaps);
  
  //=================================================================================//
  //================Second Piece of Phantom==================================//
  G4double phantom2_coordinateZ =
      supp_coordinateZ + (30 * cm - (phantomSizeZ / 2 + support_z / 2)) +
      support_z / 2;
  G4ThreeVector phantomPosition_2 =
      G4ThreeVector( 0. * mm, 0. * mm,phantom2_coordinateZ);
  // Definition of the solid volume of the Phantom
  G4Box *phantom_2 =
      new G4Box("Phantom_2", 
                phantomSizeX / 2, phantomSizeY / 2,30 * cm - (phantomSizeZ / 2 + support_z / 2));

  // Definition of the logical volume of the Phantom
  G4LogicalVolume *phantomLogicalVolume_2 =
      new G4LogicalVolume(phantom_2, phantomMaterial, "phantomLog_2", 0, 0, 0);

  // Definition of the physics volume of the Phantom

  new G4PVPlacement(0, phantomPosition_2, "phantomPhys_2",
                    phantomLogicalVolume_2, physicalTreatmentRoom, false, 0,fCheckOverlaps);

  //==========================================================================//

  G4Region *PhantomRegion = new G4Region("Phantom_reg");

  if (phantomSizeZ != 0) {
    phantomLogicalVolume->SetRegion(PhantomRegion);
    PhantomRegion->AddRootLogicalVolume(phantomLogicalVolume);
    phantomLogicalVolume->SetVisAttributes(red);
    phantomLogicalVolume->SetUserLimits(fStepLimit);
  }
  DetectorSupport->SetVisAttributes(blue);
  phantomLogicalVolume_2->SetVisAttributes(red);

  }
  
  
  
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

  airNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  // Air
  //

  G4Box *treatmentRoom = new G4Box("TreatmentRoom", worldX, worldY, worldZ);
  logicTreatmentRoom = new G4LogicalVolume(treatmentRoom, airNist,
                                           "logicTreatmentRoom", 0, 0, 0);
  physicalTreatmentRoom =
      new G4PVPlacement(0, G4ThreeVector(), "physicalTreatmentRoom",
                        logicTreatmentRoom, 0, false, 0);

  // The treatment room is invisible in the Visualisation
  logicTreatmentRoom->SetVisAttributes(G4VisAttributes::GetInvisible());

  // -----------------------------
  // detector + phantom +Default dimensions
  //------------------------------
 
  // construct collimator
  
  
  

    phantom_physical =
        ConstructPhantom(0*cm,false); //put true for inhomogeneous phantom
   



Construct_PET();





//////////////////////////============================///////////////////
  return physicalTreatmentRoom;
}







