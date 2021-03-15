#include "FlashOpticsSensitiveDetector.hh"
#include "FlashHit.hh"
FlashSensitiveDetector::FlashSensitiveDetector(G4String name,G4bool Kinetic_or_Optic):G4VSensitiveDetector(name,Kinetic_or_Optic),collectionID(-1), Kin_Opt(Kinetic_or_Optic)
{
if (Kin_Opt==0) collectionName.insert("Optic_Collection");
else if (Kin_Opt==1) collectionName.insert("Kinetic_Collection");
else G4cout <<"Error in defining Collection Name"<<endl;

}

void FlashSensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{

 if (collectionID < 0) collectionID = GetCollectionID(0);
  // Argument : order of collection// as stored in the collectionName 
hitsCollection = new FlashHitsCollection(name,collectionName[0]);
  HCE -> AddHitsCollection(collectionID, hitsCollection);
  }
  
  
 void FlashSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*ROhist) {

G4int scintillation = 0;
G4 cherenkov = 0;
 G4StepPoint* preStep = aStep->GetPreStepPoint();

  G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());

 

  FlashHit* newHit = new FlashHit();
  G4Track aTrack= aStep->GetTrack()
  
  if (Kin_Opt==0){

  newHit->SetStripNo(  touchable->GetReplicaNumber(0) );

  newHit->SetPosition( aStep->GetPreStepPoint()->GetPosition() );

  newHit->SetMomentum( aStep->GetPreStepPoint()->GetMomentum() );

  newHit->SetEnergy( aStep->GetPreStepPoint()->GetTotalEnergy() );

  newHit->SetParticle( aTrack->GetDefinition() );
  
  newHit->SetEdep(aStep->GetTotalEnergyDeposit());
  
  if( aTrack->GetTrackID()==1){ newHit->SetProcess("Primary");};
  else if(aTrack->GetTrackID()!=1 && aTrack->GetParentID()==1 ){
  if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation"){
  scintillation+=1;

  newHit->SetProcess(aTrack->GetCreatorProcess()->GetProcessName());}
  else if (aTrack->GetCreatorProcess()->GetProcessName() == "Cherenkov"){
  cherenkov+=1;

    newHit->SetProcess(aTrack->GetCreatorProcess()->GetProcessName());}
  else newHit->SetProcess("Irrelevant secondary");
    
    };
    
    
    newHit->SetScintilCount(scintillation);
    newHit->SetCherenkovCount(cherenkov);};
    
    if (Kin_Opt==1){
    
  newHit->SetEdep(aStep->GetTotalEnergyDeposit());
  
  newHit->SetStripNo(  touchable->GetReplicaNumber(0) );

  newHit->SetPosition( aStep->GetPreStepPoint()->GetPosition() );

  newHit->SetMomentum( aStep->GetPreStepPoint()->GetMomentum() );

  newHit->SetEnergy( aStep->GetPreStepPoint()->GetTotalEnergy() );

  newHit->SetParticle( aTrack->GetDefinition() );
  

  hitsCollection->insert( newHit );

 

  return true;}}
 
 void FlashSensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE) {
 
 }
