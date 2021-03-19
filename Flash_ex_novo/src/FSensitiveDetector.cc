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

G4int scintillation = 0;
G4int cherenkov = 0;
G4StepPoint* preStep = step->GetPreStepPoint();

G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());

 

  FlashHit* newHit = new FlashHit();
  G4Track* aTrack= step->GetTrack();
  
  
	
  //newHit->SetStripNo(  touchable->GetReplicaNumber(0) );
  newHit->SetParticle( aTrack->GetDefinition() );

  newHit->SetPosition( step->GetPreStepPoint()->GetPosition() );

  newHit->SetMomentum( step->GetPreStepPoint()->GetMomentum() );

  newHit->SetEnergy( step->GetPreStepPoint()->GetTotalEnergy() );

  newHit->SetParticle( aTrack->GetDefinition() );
  
  newHit->SetEdep(step->GetTotalEnergyDeposit());
  
  if (aTrack->GetTrackID()==1)
  {
    newHit->SetProcess("Primary");
  }
  else if (aTrack->GetTrackID()!=1 && aTrack->GetParentID()==1)
  {
    if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation")
    {
    scintillation+=1;
    newHit->SetProcess(aTrack->GetCreatorProcess()->GetProcessName() );
    }

    else if (aTrack->GetCreatorProcess()->GetProcessName() == "Cherenkov")
    {
    cherenkov+=1;
    newHit->SetProcess(aTrack->GetCreatorProcess()->GetProcessName());
    }
    else
    {
      newHit->SetProcess("Irrelevant secondary");
    }
  }
 
    
    newHit->SetScintilCount(scintillation);
     newHit->SetCherenkovCount(cherenkov);
         
    if ( preStep->GetStepStatus() == fGeomBoundary ){newHit->SetStatus("incident particle");}
    else {newHit->SetStatus("not incident particle");}

  
   


    


  
   

  HitsCollection->insert( newHit );

 

  return true;
}
