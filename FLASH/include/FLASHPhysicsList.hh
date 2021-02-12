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
// * work  make  any repjhgresentation or  warranty, express or implied, *
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
// Main Authors:Jake Harold Pensavalle
//////////////////////////////////////////////////////////////////////////////////////////////
#ifndef FLASHPhysicsList_h
#define FLASHPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4EmConfigurator.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class FLASHStepMax;
class FLASHPhysicsListMessenger;

class FLASHPhysicsList: public G4VModularPhysicsList
{
public:

  FLASHPhysicsList();
  virtual ~FLASHPhysicsList();

  void ConstructParticle();

  void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);
  void SetDetectorCut(G4double cut);
  void AddPhysicsList(const G4String& name);
  void ConstructProcess();

  void AddStepMax();
  FLASHStepMax* GetStepMaxProcess() {return stepMaxProcess;};
  void AddPackage(const G4String& name);

private:

  G4EmConfigurator em_config;

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;

  G4String emName;
  G4VPhysicsConstructor* emPhysicsList;
  G4VPhysicsConstructor* decPhysicsList;
  FLASHStepMax* stepMaxProcess;

  FLASHPhysicsListMessenger* pMessenger;
};

#endif
