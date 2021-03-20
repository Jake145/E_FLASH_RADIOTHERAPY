#include "FSensitiveDetector.hh"

#include "G4Step.hh"

#include "G4Track.hh"

#include "G4VProcess.hh"

#include "G4HCofThisEvent.hh"

#include "G4TouchableHistory.hh"

 FSensitiveDetector::FSensitiveDetector(G4String name):G4VSensitiveDetector(name),collectionID(-1)
{
collectionName.insert("FlashHitsCollection");

}

FSensitiveDetector::~FSensitiveDetector() {};

void FSensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{

 if (collectionID < 0) collectionID = GetCollectionID(0);
  // Argument : order of collection// as stored in the collectionName 
	HitsCollection = new FlashHitsCollection(SensitiveDetectorName,collectionName[0]);
  HCE -> AddHitsCollection(collectionID, HitsCollection);
  }
  
  
G4bool FSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*ROhist) {



  
  G4Track* aTrack= step->GetTrack();
  if(aTrack->GetTrackID!=1)return false;
  
  FlashHit* newHit = new FlashHit();
  newHit->SetParticle( aTrack->GetDefinition() );

  newHit->SetPosition( step->GetPreStepPoint()->GetPosition() );

  newHit->SetMomentum( step->GetPreStepPoint()->GetMomentum() );

  newHit->SetEnergy( step->GetPreStepPoint()->GetTotalEnergy() );

  newHit->SetParticle( aTrack->GetDefinition() );
  
  newHit->SetEdep(step->GetTotalEnergyDeposit());
  
         
    if ( preStep->GetStepStatus() == fGeomBoundary ){newHit->SetStatus("incident particle");}
    else {newHit->SetStatus("not incident particle");}

  
   


    


  
   

  HitsCollection->insert( newHit );

 

  return true;
}
