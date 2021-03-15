#include "FlashOpticsSensitiveDetector.hh"
#include "FlashHit.hh"
FlashSensitiveDetector::FlashSensitiveDetector(G4String name):G4VSensitiveDetector(name),collectionID(-1)
{
collectionName.insert("Optic_KIN_Collection");
}

void MySensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{

 if (collectionID < 0) collectionID = GetCollectionID(0);
  // Argument : order of collection// as stored in the collectionName 
hitsCollection = new FlashHitsCollection(name,collectionName[0]);
  HCE -> AddHitsCollection(collectionID, hitsCollection);
  }
  
  
 void MySensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*ROhist) {


 G4StepPoint* preStep = aStep->GetPreStepPoint();

  G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());

 

  FlashHit* newHit = new FlashHit();

  newHit->SetStripNo(  touchable->GetReplicaNumber(0) );

  newHit->SetPosition( aStep->GetPreStepPoint()->GetPosition() );

  newHit->SetMomentum( aStep->GetPreStepPoint()->GetMomentum() );

  newHit->SetEnergy( aStep->GetPreStepPoint()->GetTotalEnergy() );

  newHit->SetParticle( aStep->GetTrack()->GetDefinition() );

  hitsCollection->insert( newHit );

 

  return true;}
 
 void MySensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE) {
 
 }
