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
/// \file FlashSteppingAction.cc
/// \brief Implementation of the FlashSteppingAction class

#include "FlashSteppingAction.hh"
#include "FlashEventAction.hh"
#include "FlashDetectorConstruction.hh"
#include "FSensitiveDetector.hh"
#include "OpticFiberSD.hh"
#include "PhotoDiodeSD.hh"

#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Threading.hh"
#include <string>
#include <sstream>
//#include "FlashAnalysis.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashSteppingAction::FlashSteppingAction(FlashEventAction* eventAction)
: G4UserSteppingAction(),fEventAction(eventAction){
  /*
  //G4int ThreadNumber=G4Threading::G4GetThreadId();
std::ostringstream oss;
oss << "Kinetic_E_crystal_" << ThreadNumber << ".csv" ;
std::string filename = oss.str();
  //G4String filename=("Kinetic_E_crystal_%d.csv",ThreadNumber);
  KinEnFile.open(filename,std::ios_base::app);*/
  
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashSteppingAction::~FlashSteppingAction()
{//KinEnFile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FlashSteppingAction::UserSteppingAction(const G4Step* theStep)
{
//first of all let's kill all optical photons in the phantom
/*
G4Track* track = Step->GetTrack();
G4ParticleDefinition* particleType = track->GetDefinition();

if (particleType == G4OpticalPhoton::OpticalPhotonDefinition() && preStepPoint->GetTouchableHandle()->GetVolume()->GetName()=="phantomPhys") track->SetTrackStatus(fStopAndKill);

*/
  G4Track* theTrack = theStep->GetTrack();
  
 G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
  G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();

  G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  static G4ThreadLocal G4OpBoundaryProcess* boundary = nullptr;

 
 G4ParticleDefinition* particleType = theTrack->GetDefinition();
  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){
    //Optical photon only

    if(thePrePV->GetName()!="OF_core_phys"&& thePostPV->GetName()=="photodiode"){
     
      theTrack->SetTrackStatus(fStopAndKill);}

    //Was the photon absorbed by the absorption process
    if(thePostPoint->GetProcessDefinedStep()->GetProcessName()
       =="OpAbsorption"){
      fEventAction->IncAbsorption();

    }

    boundaryStatus=boundary->GetStatus();

    //Check to see if the partcile was actually at a boundary
    //Otherwise the boundary status may not be valid

    if(thePostPoint->GetStepStatus()==fGeomBoundary){
      if(fExpectedNextStatus==StepTooSmall){
        if(boundaryStatus!=StepTooSmall){
          G4ExceptionDescription ed;
          ed << "FlashSteppingAction::UserSteppingAction(): "
                << "No reallocation step after reflection!"
                << G4endl;
          G4Exception("FlashSteppingAction::UserSteppingAction()", "Flash",
          FatalException,ed,
          "Something is wrong with the surface normal or geometry");
        }
      }
      fExpectedNextStatus=Undefined;
      switch(boundaryStatus){
      case Absorption:
        fEventAction->IncBoundaryAbsorption();
        break;
      case Detection: //Note, this assumes that the volume causing detection
                      //is the photocathode because it is the only one with
                      //non-zero efficiency
        {
        //Triger sensitive detector manually since photon is
        //absorbed but status was Detection
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        G4String sdName="PhotoDiodeSD";
        PhotoDiodeSD* pdSD = (PhotoDiodeSD*)SDman->FindSensitiveDetector(sdName);
        if(pdSD)pdSD->ProcessHits_constStep(theStep, nullptr);
        //trackInformation->AddTrackStatusFlag(hitPMT);
        break;
        }
      case FresnelReflection:
      case TotalInternalReflection:
      case LambertianReflection:
      case LobeReflection:
      case SpikeReflection:
      case BackScattering:
        //trackInformation->IncReflections();
        fExpectedNextStatus=StepTooSmall;
        break;
      default:
        break;
      }
      
    }
  }
 
 /*
 if(KinEnFile.is_open()){
 
 KinEnFile<< eventid<< "\t" << kineticEnergy << "\t" <<trackID<<G4endl;
 }*/
 	}
 
 
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

