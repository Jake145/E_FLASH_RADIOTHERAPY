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

#include "FLASHGeometryController.hh"
#include "FLASHDetectorConstruction.hh"
//#include "Applicator1point6BeamLine.hh" //For the first version only the max SSD is available
//#include "Applicator7BeamLine.hh"
//#include "Applicator30BeamLine.hh"

#include "Applicator80BeamLine.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"

#endif
#include "FLASHGeometryMessenger.hh"

/////////////////////////////////////////////////////////////////////////////
FLASHGeometryController::FLASHGeometryController()
{}

/////////////////////////////////////////////////////////////////////////////
FLASHGeometryController::~FLASHGeometryController()
{}

/////////////////////////////////////////////////////////////////////////////
void FLASHGeometryController::SetGeometry(G4String name)
{

    
     if (name == "app80")
    {
	registerGeometry(new Applicator80BeamLine());
	G4cout <<"Applicator 80 cm geometry activated" << G4endl;
      
    } 
  /*   
    else if (name == "app30")
    {
	registerGeometry(new Applicator30BeamLine());
	G4cout <<"Applicator 30 cm geometry activated" << G4endl;
      
    } 
    else if (name == "app7")
    {
	registerGeometry(new Applicator7BeamLine());
	G4cout <<"Applicator 7 cm geometry activated" << G4endl;
      
    } 
    else if (name == "app1point6")
    {
	registerGeometry(new Applicator1point6BeamLine());
	G4cout <<"Applicator 1.6 cm geometry activated" << G4endl;
      
    } */
    else if(name == "default") 
        
    {
	registerGeometry(new Applicator80BeamLine());
    } 
    else
    {
	G4cout <<"Unknown geometry: " << name << ". Geometry not changed." << G4endl;
    }
}
	
/////////////////////////////////////////////////////////////////////////////
void FLASHGeometryController::registerGeometry(G4VUserDetectorConstruction *detector)
{
	G4RunManager *runManager = G4RunManager::GetRunManager();
	runManager->SetUserInitialization(detector);
	runManager->GeometryHasBeenModified();
}
/////////////////////////////////////////////////////////////////////////////

void FLASHGeometryController::UpdateGeometry()
{
        G4RunManager::GetRunManager() -> GeometryHasBeenModified();

}




