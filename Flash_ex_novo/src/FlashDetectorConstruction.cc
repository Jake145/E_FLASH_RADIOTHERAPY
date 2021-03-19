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

#include "G4Region.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

#include "G4SystemOfUnits.hh"
#include "Applicator80BeamLine.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "FSensitiveDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PhysicalConstants.hh"
//#include "G4VSensitiveDetector.hh"
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
  
  G4Material* LYSO = new G4Material("Lu2Y2SiO5", 7.1*g/cm3, 4);
  LYSO->AddElement(Lu, 2);
  LYSO->AddElement(Si, 1);
  LYSO->AddElement(O , 5);
  LYSO->AddElement(Y,2);
  
 G4double energy[]    = {3.061*eV, 2.952*eV, 2.844*eV,2.689*eV,2.551*eV,2.403*eV,2.271*eV};
  G4double rindex[]    = {1.833, 1.827, 1.822,1.818,1.813,1.810,1.806};
 G4double absorption[] = {136.2*nm, 142.8*nm, 149.84*nm,160.84*nm,160.84*nm,171.84*nm,185.04*nm,198.24*nm}; 
   const G4int nEntries = sizeof(energy)/sizeof(G4double);
  G4MaterialPropertiesTable*MPT =new G4MaterialPropertiesTable();
  // property independent of energy
  MPT->AddConstProperty("SCINTILLATIONYIELD", 27600./MeV);
  MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 45.*ns);
  MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  // properties that depend on energy
  MPT->AddProperty("RINDEX", energy, rindex,nEntries)->SetSpline(true);
  MPT->AddProperty("ABSLENGTH", energy, absorption,nEntries)->SetSpline(true);
  LYSO->SetMaterialPropertiesTable(MPT);
G4cout << "LYSO G4MaterialPropertiesTable" << G4endl;
  MPT->DumpTable();

    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FlashDetectorConstruction::ConstructPhantom(){


G4NistManager* nist = G4NistManager::Instance();
G4Material* phantomMaterial = nist->FindOrBuildMaterial("G4_WATER");

// ------------ Generate & Add Material Properties Table ------------
// 
/*
  G4double photonEnergy[] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };

  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

//
// Water
//
  G4double refractiveIndex1[] =
            { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
              1.346,  1.3465, 1.347,  1.3475, 1.348,
              1.3485, 1.3492, 1.35,   1.3505, 1.351,
              1.3518, 1.3522, 1.3530, 1.3535, 1.354,
              1.3545, 1.355,  1.3555, 1.356,  1.3568,
              1.3572, 1.358,  1.3585, 1.359,  1.3595,
              1.36,   1.3608};


  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

  G4double absorption[] =
           {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
           15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
           45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
           52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
           30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
           17.500*m, 14.500*m };

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

  assert(sizeof(scintilSlow) == sizeof(photonEnergy));
  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
        ->SetSpline(true);
  myMPT1->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
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
     1.56962*eV, 1.58974*eV, 1.61039*eV, 1.63157*eV,
     1.65333*eV, 1.67567*eV, 1.69863*eV, 1.72222*eV,
     1.74647*eV, 1.77142*eV, 1.7971 *eV, 1.82352*eV,
     1.85074*eV, 1.87878*eV, 1.90769*eV, 1.93749*eV,
     1.96825*eV, 1.99999*eV, 2.03278*eV, 2.06666*eV,
     2.10169*eV, 2.13793*eV, 2.17543*eV, 2.21428*eV,
     2.25454*eV, 2.29629*eV, 2.33962*eV, 2.38461*eV,
     2.43137*eV, 2.47999*eV, 2.53061*eV, 2.58333*eV,
     2.63829*eV, 2.69565*eV, 2.75555*eV, 2.81817*eV,
     2.88371*eV, 2.95237*eV, 3.02438*eV, 3.09999*eV,
     3.17948*eV, 3.26315*eV, 3.35134*eV, 3.44444*eV,
     3.54285*eV, 3.64705*eV, 3.75757*eV, 3.87499*eV,
     3.99999*eV, 4.13332*eV, 4.27585*eV, 4.42856*eV,
     4.59258*eV, 4.76922*eV, 4.95999*eV, 5.16665*eV,
     5.39129*eV, 5.63635*eV, 5.90475*eV, 6.19998*eV
  };

  const G4int numentries_water = sizeof(energy_water)/sizeof(G4double);

  //assume 100 times larger than the rayleigh scattering for now.
  G4double mie_water[] = {
     167024.4*m, 158726.7*m, 150742  *m,
     143062.5*m, 135680.2*m, 128587.4*m,
     121776.3*m, 115239.5*m, 108969.5*m,
     102958.8*m, 97200.35*m, 91686.86*m,
     86411.33*m, 81366.79*m, 76546.42*m,
     71943.46*m, 67551.29*m, 63363.36*m,
     59373.25*m, 55574.61*m, 51961.24*m,
     48527.00*m, 45265.87*m, 42171.94*m,
     39239.39*m, 36462.50*m, 33835.68*m,
     31353.41*m, 29010.30*m, 26801.03*m,
     24720.42*m, 22763.36*m, 20924.88*m,
     19200.07*m, 17584.16*m, 16072.45*m,
     14660.38*m, 13343.46*m, 12117.33*m,
     10977.70*m, 9920.416*m, 8941.407*m,
     8036.711*m, 7202.470*m, 6434.927*m,
     5730.429*m, 5085.425*m, 4496.467*m,
     3960.210*m, 3473.413*m, 3032.937*m,
     2635.746*m, 2278.907*m, 1959.588*m,
     1675.064*m, 1422.710*m, 1200.004*m,
     1004.528*m, 833.9666*m, 686.1063*m
  };

  assert(sizeof(mie_water) == sizeof(energy_water));
  // gforward, gbackward, forward backward ratio
  G4double mie_water_const[3]={0.99,0.99,0.8};

  myMPT1->AddProperty("MIEHG",energy_water,mie_water,numentries_water)
        ->SetSpline(true);
  myMPT1->AddConstProperty("MIEHG_FORWARD",mie_water_const[0]);
  myMPT1->AddConstProperty("MIEHG_BACKWARD",mie_water_const[1]);
  myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO",mie_water_const[2]);

  G4cout << "Water G4MaterialPropertiesTable" << G4endl;
  myMPT1->DumpTable();

  phantomMaterial->SetMaterialPropertiesTable(myMPT1);
  phantomMaterial->GetIonisation()->SetBirksConstant(0.126*mm/MeV); */
  
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
    phant_phys =new G4PVPlacement(0,
	                                    phantomPosition,
					    "phantomPhys",
					    phantomLogicalVolume,
					    physicalTreatmentRoom,
					    false,
					    0);
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
  myMPT2->AddProperty("RINDEX", photonEnergy_air, refractiveIndex2, nEntries_air);

  G4cout << "Air G4MaterialPropertiesTable" << G4endl;
  myMPT2->DumpTable();

  airNist->SetMaterialPropertiesTable(myMPT2); */
  
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
  G4Material* PMMA = nist->FindOrBuildMaterial("G4_PLEXIGLASS", 
  isotopes);
  
  //--------------------------------------------------
  //  PMMA for fibers
  //--------------------------------------------------
G4double photonEnergy_fib[] =
  {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
   2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
   2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
   2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
   2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
   2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
   2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
   3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
   3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
   3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

  const G4int nEntries_fib = sizeof(photonEnergy_fib)/sizeof(G4double);
  
  G4double refractiveIndex_fiber[] =
  { 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
    1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
    1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
    1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
    1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60};

  assert(sizeof(refractiveIndex_fiber) == sizeof(photonEnergy_fib));

  G4double abs_fiber[] =
  {5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
   5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
   5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,1.10*m,
   1.10*m,1.10*m,1.10*m,1.10*m,1.10*m,1.10*m, 1.*mm, 1.*mm, 1.*mm, 1.*mm,
    1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm};

  assert(sizeof(abs_fiber) == sizeof(photonEnergy_fib));

 /* G4double emission_fib[] =
  {0.05, 0.10, 0.30, 0.50, 0.75, 1.00, 1.50, 1.85, 2.30, 2.75,
   3.25, 3.80, 4.50, 5.20, 6.00, 7.00, 8.50, 9.50, 11.1, 12.4,
   12.9, 13.0, 12.8, 12.3, 11.1, 11.0, 12.0, 11.0, 17.0, 16.9,
   15.0, 9.00, 2.50, 1.00, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00,
   0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

  assert(sizeof(emission_fib) == sizeof(photonEnergy_fib)); */

  // Add entries into properties table
  G4MaterialPropertiesTable* mpt_fiber = new G4MaterialPropertiesTable();
  mpt_fiber->
           AddProperty("RINDEX",photonEnergy_fib,refractiveIndex_fiber,nEntries_fib);

  mpt_fiber->AddProperty("ABSLENGTH",photonEnergy_fib,abs_fiber,nEntries_fib);
  //mptWLSfiber->AddProperty("WLSCOMPONENT",photonEnergy_fib,emissionFib,nEntries_fib);


  PMMA->SetMaterialPropertiesTable(mpt_fiber);
  G4cout << "PMMA G4MaterialPropertiesTable" << G4endl;
  mpt_fiber->DumpTable();

  
  G4Material* PE = nist->FindOrBuildMaterial("G4_POLYETHYLENE", 
  isotopes);  
  
  G4double refractiveIndex_Clad[] =
  { 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49};

  assert(sizeof(refractiveIndex_Clad) == sizeof(photonEnergy_fib));

  G4double abs_Clad[] =
  {20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m};

  assert(sizeof(abs_Clad) == sizeof(photonEnergy_fib));

  // Add entries into properties table
  G4MaterialPropertiesTable* mptClad = new G4MaterialPropertiesTable();
  mptClad->AddProperty("RINDEX",photonEnergy_fib,refractiveIndex_Clad,nEntries_fib);
  mptClad->AddProperty("ABSLENGTH",photonEnergy_fib,abs_Clad,nEntries_fib);

  PE->SetMaterialPropertiesTable(mptClad);
  
  G4Material* TEFLON = nist->FindOrBuildMaterial("G4_TEFLON",isotopes);  
                   
G4double photonEnergy_teflon[] =
            { 7.897*eV,7.208*eV, 6.702*eV,  4.999*eV};
  G4double refractiveIndex3[] =
            { 1.432, 1.308, 1.364, 1.329};
const G4int nEntries_teflon = sizeof(photonEnergy_teflon)/sizeof(G4double);
  G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();
  myMPT3->AddProperty("RINDEX", photonEnergy_teflon, refractiveIndex3, nEntries_teflon);

  G4cout << "Teflon G4MaterialPropertiesTable" << G4endl;
  myMPT3->DumpTable();

  TEFLON->SetMaterialPropertiesTable(myMPT3);
  
  
  
G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, dZ/2);
    //                 
   logicCryst = 
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name
G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg); 
    
    G4ThreeVector position = G4ThreeVector(-80.0*mm,  0.*mm,0.*mm);     
    G4Transform3D transform = G4Transform3D(rotm,position);
                                    
   G4VPhysicalVolume* phys_cryst = new G4PVPlacement(transform,logicCryst,            
                      "crystalphys",         
                      phantomLogicalVolume,
                      false,0);  
//OOOOOOOOOOOooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooo
//costruiamo il wrap in teflon

 G4double wrap=0.2*mm; 
G4Box* big_cover_top_bottom = new G4Box("big_cover", (dZ+wrap)/2+wrap/2, wrap, (dX+wrap)/2);
G4Box* big_cover_sides = new G4Box("sides_cover", wrap/2, dY/2, (dX+wrap)/2);
G4Box* little_cover = new G4Box("little_cover",dZ/2,dY/2,wrap);



                     
  G4LogicalVolume* logicwrapper_long = 
    new G4LogicalVolume(big_cover_top_bottom,          //its solid
                        TEFLON,           //its material
                        "wrapperLVbig");        //its name
G4RotationMatrix rotm_wrap_up  = G4RotationMatrix();
    rotm_wrap_up.rotateY(90*deg); 
    
    G4ThreeVector position_wrap_up = G4ThreeVector(-3*wrap/2,(dY)/2+wrap,0*mm);     
    G4Transform3D transform_warp_up = G4Transform3D(rotm_wrap_up,position_wrap_up);
                                    
   G4VPhysicalVolume* phys_wrap_up = new G4PVPlacement(transform_warp_up,logicwrapper_long,            
                      "wrapperphysup",         
                      logicCryst,
                      false,0); 

G4RotationMatrix rotm_wrap_down  = G4RotationMatrix();
    rotm_wrap_down.rotateY(90*deg); 
    
    G4ThreeVector position_wrap_down = G4ThreeVector(-3*wrap/2,-((dY)/2+wrap),0*mm);      
    G4Transform3D transform_warp_down = G4Transform3D(rotm_wrap_down,position_wrap_down);
                                    
   G4VPhysicalVolume* phys_wrap_down = new G4PVPlacement(transform_warp_down,logicwrapper_long,            
                      "wrapperphysdown",         
                      logicCryst,
                      false,0); 
                      
G4LogicalVolume* logicwrapper_side = 
    new G4LogicalVolume(big_cover_sides,          //its solid
                        TEFLON,           //its material
                        "wrapperLVside");        //its name       
G4RotationMatrix rotm_wrap_side_r  = G4RotationMatrix();
    rotm_wrap_side_r.rotateY(90*deg); 
    
    G4ThreeVector position_wrap_side_r = G4ThreeVector(-3*wrap/2,0.*mm,dZ/2+wrap/2);     
    G4Transform3D transform_warp_side_r = G4Transform3D(rotm_wrap_side_r,position_wrap_side_r);
                                    
   G4VPhysicalVolume* phys_wrap_r = new G4PVPlacement(transform_warp_side_r,logicwrapper_side,            
                      "wrapperphysright",         
                      logicCryst,
                      false,0); 

G4RotationMatrix rotm_wrap_side_l  = G4RotationMatrix();
    rotm_wrap_side_l.rotateY(90*deg); 
    
    
    G4ThreeVector position_wrap_side_l = G4ThreeVector(-3*wrap/2,0*mm,-(dZ/2+wrap/2));       
    G4Transform3D transform_warp_side_l = G4Transform3D(rotm_wrap_side_l,position_wrap_side_l);
                                    
   G4VPhysicalVolume* phys_wrap_l = new G4PVPlacement(transform_warp_side_l,logicwrapper_side,            
                      "wrapperphysleft",         
                      logicCryst,
                      false,0);        

G4LogicalVolume* logicwrapper_little = 
    new G4LogicalVolume(little_cover,          //its solid
                        TEFLON,           //its material
                        "wrapperLVlittle");        //its name
 
 G4RotationMatrix rotm_wrap_little  = G4RotationMatrix();
    rotm_wrap_little.rotateY(90*deg); 
    
    G4ThreeVector position_wrap_little = G4ThreeVector(-dX/2-wrap,0*mm,0*mm);     
    G4Transform3D transform_warp_little = G4Transform3D(rotm_wrap_little,position_wrap_little);
                                    
   G4VPhysicalVolume* phys_wrap_little = new G4PVPlacement(transform_warp_little,logicwrapper_little,            
                      "wrapperphyslittle",         
                      logicCryst,
                      false,0);           
    
    //OOOOOOOOOOOOOOOOOooooooooooooooOOOOOOOOOOOOOOOOOOOOOooooooooooooooooOOO 

  G4OpticalSurface* opteflonSurface_up = new G4OpticalSurface("teflonSurface_up");
  opteflonSurface_up->SetType(dielectric_LUTDAVIS);
  opteflonSurface_up->SetFinish(Rough_LUT);
  opteflonSurface_up->SetModel(DAVIS);

  G4LogicalBorderSurface* teflonSurface_up=
          new G4LogicalBorderSurface("teflonSurface_up",
                                 phys_wrap_up,phys_cryst,opteflonSurface_up);

  G4OpticalSurface* opticalSurface_2 = dynamic_cast <G4OpticalSurface*>
        (teflonSurface_up->GetSurface(phys_wrap_up,phys_cryst)->
                                                       GetSurfaceProperty());
  if (opticalSurface_2) opticalSurface_2->DumpInfo();  
  
  
  G4OpticalSurface* opteflonSurface_down = new G4OpticalSurface("teflonSurface_down");
  opteflonSurface_down->SetType(dielectric_LUTDAVIS);
  opteflonSurface_down->SetFinish(Rough_LUT);
  opteflonSurface_down->SetModel(DAVIS);

  G4LogicalBorderSurface* teflonSurface_down=
          new G4LogicalBorderSurface("teflonSurface_down",
                                 phys_wrap_down,phys_cryst,opteflonSurface_down);

  G4OpticalSurface* opticalSurface_3 = dynamic_cast <G4OpticalSurface*>
        (teflonSurface_down->GetSurface(phys_wrap_down,phys_cryst)->
                                                       GetSurfaceProperty());
  if (opticalSurface_3) opticalSurface_3->DumpInfo();  
  
  
  G4OpticalSurface* opteflonSurface_right = new G4OpticalSurface("teflonSurface_right");
  opteflonSurface_right->SetType(dielectric_LUTDAVIS);
  opteflonSurface_right->SetFinish(Rough_LUT);
  opteflonSurface_right->SetModel(DAVIS);

  G4LogicalBorderSurface* teflonSurface_right=
          new G4LogicalBorderSurface("teflonSurface_right",
                                 phys_wrap_r,phys_cryst,opteflonSurface_right);

  G4OpticalSurface* opticalSurface_4 = dynamic_cast <G4OpticalSurface*>
        (teflonSurface_right->GetSurface(phys_wrap_r,phys_cryst)->
                                                       GetSurfaceProperty());
  if (opticalSurface_4) opticalSurface_4->DumpInfo();  
  
  G4OpticalSurface* opteflonSurface_left = new G4OpticalSurface("teflonSurface_left");
  opteflonSurface_left->SetType(dielectric_LUTDAVIS);
  opteflonSurface_left->SetFinish(Rough_LUT);
  opteflonSurface_left->SetModel(DAVIS);

  G4LogicalBorderSurface* teflonSurface_left=
          new G4LogicalBorderSurface("teflonSurface_right",
                                 phys_wrap_l,phys_cryst,opteflonSurface_left);

  G4OpticalSurface* opticalSurface_5 = dynamic_cast <G4OpticalSurface*>
        (teflonSurface_left->GetSurface(phys_wrap_l,phys_cryst)->
                                                       GetSurfaceProperty());
  if (opticalSurface_5) opticalSurface_5->DumpInfo();  
  
 G4OpticalSurface* opteflonSurface_back = new G4OpticalSurface("teflonSurface_back");
  opteflonSurface_back->SetType(dielectric_LUTDAVIS);
  opteflonSurface_back->SetFinish(Rough_LUT);
  opteflonSurface_back->SetModel(DAVIS);

  G4LogicalBorderSurface* teflonSurface_back=
          new G4LogicalBorderSurface("teflonSurface_back",
                                 phys_wrap_little,phys_cryst,opteflonSurface_back);

  G4OpticalSurface* opticalSurface_6 = dynamic_cast <G4OpticalSurface*>
        (teflonSurface_back->GetSurface(phys_wrap_little,phys_cryst)->
                                                       GetSurfaceProperty());
  if (opticalSurface_6) opticalSurface_6->DumpInfo();  

//OOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOO  
                      
  G4Tubs* opticfiber_core =
    new G4Tubs("OF_core", 0.*cm, opticfiber_core_radius, opticfiber_core_dx/2, 0., CLHEP::twopi);
      
   opticfiber_core_log =                         
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
    
    

  //
   G4VPhysicalVolume* physclad = new G4PVPlacement(transform_clad,                     
                       //position
                      opticfiber_clad_log,             
                      "outerfiber",                
                      logicCryst,false,0);                 
                            
   G4VPhysicalVolume* physcore = new G4PVPlacement(transform_opt,                     
                           //position
                          opticfiber_core_log,            

                          "OF_clad_phys",                

                          opticfiber_clad_log,false,0); 
//OOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooo
// Definiamo ora le ultime tre superfici ottiche interessanti


G4OpticalSurface* opcore_scint = new G4OpticalSurface("OpticFiberandScintillator");
  opcore_scint->SetType(dielectric_LUTDAVIS);
  opcore_scint->SetFinish(Rough_LUT);
  opcore_scint->SetModel(DAVIS);

  G4LogicalBorderSurface* core_scint=
          new G4LogicalBorderSurface("opfibscint",
                                 physcore,phys_cryst,opcore_scint);

  G4OpticalSurface* opticalSurface_7 = dynamic_cast <G4OpticalSurface*>
        (core_scint->GetSurface(physcore,phys_cryst)->
                                                       GetSurfaceProperty());
  if (opticalSurface_7) opticalSurface_7->DumpInfo();   
  
  
  G4OpticalSurface* opcore_clad = new G4OpticalSurface("OpticFiberandClad");
  opcore_clad->SetType(dielectric_LUTDAVIS);
  opcore_clad->SetFinish(Rough_LUT);
  opcore_clad->SetModel(DAVIS);

  G4LogicalBorderSurface* core_clad=
          new G4LogicalBorderSurface("opfibclad",
                                 physcore,physclad,opcore_clad);

  G4OpticalSurface* opticalSurface_8 = dynamic_cast <G4OpticalSurface*>
        (core_clad->GetSurface(physcore,physclad)->
                                                       GetSurfaceProperty());
  if (opticalSurface_8) opticalSurface_8->DumpInfo();   
  

 /* 

G4OpticalSurface* opphantom_clad = new G4OpticalSurface("PhantomandClad");
  opphantom_clad->SetType(dielectric_LUTDAVIS);
  opphantom_clad->SetFinish(Rough_LUT);
  opphantom_clad->SetModel(DAVIS);

  G4LogicalBorderSurface* phantom_clad=
          new G4LogicalBorderSurface("opfibclad",
                                 phant_phys,physclad,opphantom_clad);

  G4OpticalSurface* opticalSurface_9 = dynamic_cast <G4OpticalSurface*>
        (phantom_clad->GetSurface(phant_phys,physclad)->
                                                       GetSurfaceProperty());
  if (opticalSurface_9) opticalSurface_9->DumpInfo();   

        */ 


   
    G4VisAttributes * skyBlue1 = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
     skyBlue1->SetVisibility(true);
   // logicCryst -> SetVisAttributes(red);
    logicwrapper_long->SetVisAttributes(skyBlue1);
        logicwrapper_side->SetVisAttributes(skyBlue1);
            logicwrapper_little->SetVisAttributes(skyBlue1);
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
  FSensitiveDetector* cryst = new FSensitiveDetector("crystalSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
  // G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  // cryst->RegisterPrimitive(primitiv1);
  SetSensitiveDetector("CrystalLV",cryst);
  
  
 /* 
  G4String SDname_cr_kin;
    G4String SDname_cr_opt;


  FSensitiveDetector* sd_cr_kin = new FSensitiveDetector(SDname_cr_kin = "Kinetic_crystal",true);
  G4SDManager* SDman_cr_k = G4SDManager::GetSDMpointer();

  SDman_cr_k->AddNewDetector( sd_cr_kin );

  logicCryst->SetSensitiveDetector(sd_cr_kin);
  
  FSensitiveDetector* sd_cr_opt = new FSensitiveDetector(SDname_cr_opt = "Optic_crystal",false);
  G4SDManager* SDman_cr_opt = G4SDManager::GetSDMpointer();

  SDman_cr_opt->AddNewDetector( sd_cr_opt );

  logicCryst->SetSensitiveDetector(sd_cr_opt);
  */
  
  G4String SDname_of_opt="Optic_fiber";
   G4SDManager* SDman_of_opt = G4SDManager::GetSDMpointer();
   sd_of_opt = new FSensitiveDetector(SDname_of_opt);


  SDman_of_opt->AddNewDetector( sd_of_opt );
  SetSensitiveDetector("OF_core_LV",sd_of_opt,true);
  //opticfiber_core_log->SetSensitiveDetector(sd_of_opt);
  
  
  //G4MultiFunctionalDetector* optfiber = new G4MultiFunctionalDetector("fiberSD");
 // G4SDManager::GetSDMpointer()->AddNewDetector(optfiber); 
  //SetSensitiveDetector("OF_core_LV",optfiber);
  
  // declare patient as a MultiFunctionalDetector scorer
  //  
  //G4MultiFunctionalDetector* patient = new G4MultiFunctionalDetector("Phantom");
  //G4SDManager::GetSDMpointer()->AddNewDetector(phantom);
  //G4VPrimitiveScorer* primitiv2 = new G4PSDoseDeposit("dose");
  //patient->RegisterPrimitive(primitiv2);
  //SetSensitiveDetector("phantomLog",phantom);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
