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
//

#include "FlashPrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Threading.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"
#include "globals.hh"
#include <assert.h>
FlashPrimaryGeneratorAction::FlashPrimaryGeneratorAction() {
  
  
  
  i=G4Threading::G4GetThreadId(); // helper index
  
  n=G4Threading::G4GetNumberOfCores(); // number of pulses of the beam, which is equal to the number of threads
  
  dt= 0.2 ; // duration of pulse
  
  T= 10; // total irradiation time
  
  dr = T/n - dt; // time between pulses
  
  assert (dr > 0);

  particleGun = new G4GeneralParticleSource();

  SetDefaultPrimaryParticle();
}

FlashPrimaryGeneratorAction::~FlashPrimaryGeneratorAction() {
  delete particleGun;
}

void FlashPrimaryGeneratorAction::SetDefaultPrimaryParticle() {
  
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *particle = particleTable->FindParticle("e-");
  particleGun->SetParticleDefinition(particle);

  G4double defaultX0 = -1000.0 * CLHEP::mm;
  X0 = defaultX0;

  G4double defaultY0 = 0.0 * CLHEP::mm;
  Y0 = defaultY0;

  G4double defaultZ0 = 0.0 * CLHEP::mm;
  Z0 = defaultZ0;

  G4double defaultsigmaY = 1.115 * CLHEP::mm;
  sigmaY = 2 * CLHEP::mm;

  G4double defaultsigmaZ = 1.115 * CLHEP::mm;
  sigmaX = 2 * CLHEP::mm;

  G4double defaultTheta = 3 * CLHEP::deg;
  Theta = defaultTheta;
}

void FlashPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
G4bool Source = false;
Pencil = true;
G4bool Rate = true;
if (Source){

	particleGun->GeneratePrimaryVertex(anEvent);
}

else{
  G4double x = X0;
  G4double y = Y0;
  G4double z = Z0;

  if (sigmaY > 0.0) {
    y += G4RandGauss::shoot(Y0, sigmaY);
  }

  if (sigmaX > 0.0) {
    x += G4RandGauss::shoot(X0, sigmaX);
  }

  particleGun->SetParticlePosition(G4ThreeVector(x, y, z));
if (Pencil == false){
  G4double Mx;
  G4double My;
  G4double Mz;
  G4double condition;

  while (true) {

    Mz = CLHEP::RandFlat::shoot(0.7, 1);
    My = CLHEP::RandFlat::shoot(-0.3, 0.3); // ranges good for 0<Theta<20
    Mx = CLHEP::RandFlat::shoot(-0.3, 0.3);

    condition = std::sqrt(Mx * Mx + My * My + Mz * Mz);

    if (condition < 1) {
      Mx = Mx / condition;
      My = My / condition;
      Mz = Mz / condition;

      if (Mz > std::cos(Theta)) {
        break;
      }
    }
  }

  particleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(
      G4ThreeVector(Mx, My, Mz));}
    else if(Pencil==true){
     particleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(
      G4ThreeVector(0., 0., 1.));}
      
      if (Rate){

      
       assert (i*(dt+dr) + dr > 0);
       assert ((i+1)*(dt+dr) >0);
       
      G4double Shot_time = CLHEP::RandFlat::shoot(i*(dt+dr) + dr,(i+1)*(dt+dr));
      
      
      particleGun->SetParticleTime(Shot_time*ns);
      
      
      }
      

  particleGun->GeneratePrimaryVertex(anEvent);
  
}}
