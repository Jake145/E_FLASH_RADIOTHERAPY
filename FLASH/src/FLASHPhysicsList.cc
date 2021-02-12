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
//////////////////////////////////////////////////////////////////////////////////////////////
//

#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh" 
#include "G4Region.hh"     
#include "G4RegionStore.hh"   
#include "FLASHPhysicsList.hh"
#include "FLASHPhysicsListMessenger.hh"
#include "FLASHStepMax.hh"
#include "G4VPhysicsConstructor.hh"

// Physic lists (contained inside the Geant4 source code, in the 'physicslists folder')
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"  
#include "G4EmPenelopePhysics.hh"   
#include "G4EmExtraPhysics.hh"   
#include "G4EmProcessOptions.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProcessManager.hh"
#include "globals.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4UnitsTable.hh"
#include "G4DecayPhysics.hh"
#include "G4AutoLock.hh"

namespace { G4Mutex myFLASHPhysicsList = G4MUTEX_INITIALIZER; }

FLASHPhysicsList::FLASHPhysicsList() : G4VModularPhysicsList()
{
  defaultCutValue = 0.1 *mm; 
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  stepMaxProcess  = 0;

  pMessenger = new FLASHPhysicsListMessenger(this);

  SetVerboseLevel(1);

  // EM physics
  //emPhysicsList = new G4EmStandardPhysics_option4(1);
  emPhysicsList = new G4EmPenelopePhysics(1);
  //emName = G4String("emstandard_opt4");
  emName = G4String("empenelope");
  decPhysicsList = new G4DecayPhysics();
}

FLASHPhysicsList::~FLASHPhysicsList()
{
  delete pMessenger;
  delete emPhysicsList;
  delete decPhysicsList;
}

void FLASHPhysicsList::ConstructParticle()
{
  decPhysicsList->ConstructParticle();
}

void FLASHPhysicsList::ConstructProcess()
{
  // transportation
  AddTransportation();

  // electromagnetic physics list
  emPhysicsList->ConstructProcess();
  em_config.AddModels();

  // step limitation (as a full process)
  //
  AddStepMax();
}

void FLASHPhysicsList::AddPhysicsList(const G4String& name)
{

    if (verboseLevel>1) {
	G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
    }
    if (name == emName) return;

  if (name == "standard_opt3") {
	emName = name;
	delete emPhysicsList;
	emPhysicsList = new G4EmStandardPhysics_option3();
	G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option3" << G4endl;


    } else if (name == "livermore") {   
	emName = name;
	delete emPhysicsList;
	emPhysicsList = new G4EmLivermorePhysics();
	G4RunManager::GetRunManager()-> PhysicsHasBeenModified();
	G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmLivermorePhysics" << G4endl;

    } else if (name == "standard_opt4") {    
	emName = name;
	delete emPhysicsList;
	emPhysicsList = new G4EmStandardPhysics_option4();
	G4RunManager::GetRunManager()-> PhysicsHasBeenModified();
	G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option4" << G4endl;}

 else if (name == "standard_opt4") {
	emName = name;
	delete emPhysicsList;
	emPhysicsList = new G4EmStandardPhysics_option4();
	G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option4" << G4endl;
        }

    else {

	G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
	    << " is not defined"
	    << G4endl;
    }
}

void FLASHPhysicsList::AddStepMax()
{
  G4AutoLock lock(&myFLASHPhysicsList);
  
  // Step limitation seen as a process
  stepMaxProcess = new FLASHStepMax();

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (stepMaxProcess->IsApplicable(*particle) && pmanager)
      {
	pmanager ->AddDiscreteProcess(stepMaxProcess);
      }
  }
}

void FLASHPhysicsList::SetCuts()
{

  if (verboseLevel >0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

  // Set cuts for detector
  SetDetectorCut(defaultCutValue); 
  if (verboseLevel>0) DumpCutValuesTable();
}

void FLASHPhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

void FLASHPhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

void FLASHPhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

void FLASHPhysicsList::SetDetectorCut(G4double cut)
{
  G4String regionName = "DetectorLog";
  G4Region* region = G4RegionStore::GetInstance()->GetRegion(regionName);

  G4ProductionCuts* cuts = new G4ProductionCuts ;
  cuts -> SetProductionCut(cut,G4ProductionCuts::GetIndex("gamma"));
  cuts -> SetProductionCut(cut,G4ProductionCuts::GetIndex("e-"));
  cuts -> SetProductionCut(cut,G4ProductionCuts::GetIndex("e+"));
  region -> SetProductionCuts(cuts);
}

