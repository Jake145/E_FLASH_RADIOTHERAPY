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
// This is the second version of Flash, a Geant4-based application
//
// 
//////////////////////////////////////////////////////////////////////////////////////////////

#include "G4SystemOfUnits.hh"
#include "FlashPrimaryGeneratorAction.hh"


#include "globals.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
  
FlashPrimaryGeneratorAction::FlashPrimaryGeneratorAction()
{
  
  particleGun  = new G4GeneralParticleSource();

  SetDefaultPrimaryParticle();  
}  

FlashPrimaryGeneratorAction::~FlashPrimaryGeneratorAction()
{
  delete particleGun;


}
  
void FlashPrimaryGeneratorAction::SetDefaultPrimaryParticle()
{    
  // ****************************
  // Default primary particle
  // ****************************
  
  // Define primary particles: electrons 
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable -> FindParticle("e-");   
  particleGun -> SetParticleDefinition(particle); 

  //qui definiamo la posizione della sorgente. Il fascio ha forma gaussiana e divergenza 6 gradi. Lo spettro energetico viene dato dalla macro a gps
  
  G4double defaultX0 = -1000.0 *CLHEP::mm;                 
  X0 = defaultX0;

  G4double defaultY0 = 0.0 *CLHEP::mm;  
  Y0 = defaultY0;

  G4double defaultZ0 = 0.0 *CLHEP::mm;  
  Z0 = defaultZ0;

  G4double defaultsigmaY = 1.115  *CLHEP::mm;  
  sigmaY = defaultsigmaY;

  G4double defaultsigmaZ = 1.115  *CLHEP::mm;  
  sigmaZ = defaultsigmaZ;

  
  G4double defaultTheta = 5.0 *CLHEP::deg;  
  Theta = defaultTheta;

}

void FlashPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // ****************************************
  // Set the beam angular apread 
  // and spot size
  // beam spot size
  // ****************************************

  // Set the position of the primary particles
  G4double x = X0;
  G4double y = Y0;
  G4double z = Z0;

  if ( sigmaY > 0.0 )
    {
      y += G4RandGauss::shoot( Y0, sigmaY );
    }
 
  if ( sigmaZ > 0.0 )
    {
      z += G4RandGauss::shoot( Z0, sigmaZ );
    }
  
  particleGun -> SetParticlePosition(G4ThreeVector( x , y , z ) );
 
  // ********************************************
  // Set the beam energy and energy spread
  // ********************************************

               
  G4double Mx;    
  G4double My;
  G4double Mz;
  G4double condition;
  
while (true)  {

  

  Mx =  CLHEP::RandFlat::shoot(0.7,1);
  My =  CLHEP::RandFlat::shoot(-0.3,0.3); // ranges good for 0<Theta<20
  Mz =  CLHEP::RandFlat::shoot(-0.3,0.3);

  condition = std::sqrt(Mx*Mx + My*My + Mz*Mz);

 
  if (condition < 1)  {
    Mx = Mx/condition;
    My = My/condition;
    Mz = Mz/condition;


    if (Mx > std::cos(Theta)) { 
      break;
        }
    }
}
  
 
  particleGun ->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection( G4ThreeVector(Mx,My,Mz) );
  

  // Generate a primary particle
  particleGun -> GeneratePrimaryVertex( anEvent ); 
} 




