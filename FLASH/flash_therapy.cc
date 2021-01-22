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
//
// ----------------------------------------------------------------------------
//                 GEANT 4 - FLASH V1
// ----------------------------------------------------------------------------
// Main Authors:Jake Harold Pensavalle
// Based on IORT by G.Russo, C.Casarino, G.C. Candiano, G.A.P. Cirrone, F.Romano
// 
//  ----------------------------------------------------------------------------

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4UImanager.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "FLASHEventAction.hh"
#include "FLASHPhysicsList.hh"
#include "FLASHPrimaryGeneratorAction.hh"
#include "FLASHRunAction.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4UImessenger.hh"
#include "globals.hh"
#include "FLASHSteppingAction.hh"
#include "FLASHGeometryController.hh"
#include "FLASHGeometryMessenger.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include <ctime>
#include "G4ScoringManager.hh"

//////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc ,char ** argv)
{
  // Set the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());

//  G4int seed = time(NULL); 
 // CLHEP::HepRandom::setTheSeed(seed);
  
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  G4ScoringManager::GetScoringManager();
  // Scoring mesh

  // Geometry controller is responsible for instantiating the
  // geometries. All geometry specific setup tasks are now in class
  // FLASHGeometryController.
  FLASHGeometryController *geometryController = new FLASHGeometryController();
	
  // Connect the geometry controller to the G4 user interface
  FLASHGeometryMessenger *geometryMessenger = new FLASHGeometryMessenger(geometryController);
		
  // Initialize the default FLASH geometry
  geometryController->SetGeometry("default");  

  runManager->SetUserInitialization(new FLASHPhysicsList());

  // Initialize the primary particles
  FLASHPrimaryGeneratorAction *pPrimaryGenerator = new FLASHPrimaryGeneratorAction();
  runManager -> SetUserAction(pPrimaryGenerator);
	
  // Optional UserActions: run, event, stepping
  FLASHRunAction* pRunAction = new FLASHRunAction();
  runManager -> SetUserAction(pRunAction);
	
  FLASHEventAction* pEventAction = new FLASHEventAction();
  runManager -> SetUserAction(pEventAction);
	
  FLASHSteppingAction* steppingAction = new FLASHSteppingAction(pRunAction); 
  runManager -> SetUserAction(steppingAction);    

  // Visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager -> Initialize();
	
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);    
    }
  else
    {  // interactive mode : define UI session
       
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
  
      UImanager->ApplyCommand("/control/execute defaultMacro.mac"); 

      if (ui->IsGUI())
      ui->SessionStart();
      delete ui;
    }


  delete visManager;           

  delete geometryMessenger;
  delete geometryController;
  delete runManager;
  return 0;
  
}

