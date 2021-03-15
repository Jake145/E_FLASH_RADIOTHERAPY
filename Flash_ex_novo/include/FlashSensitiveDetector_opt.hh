#ifndef FlashSensitiveDetector_h
#define FlashsensitiveDetector_h 1
#include "G4VSensitiveDetector.hh"
#include "FlashHit.hh"
class FlashSensitiveDetector_opt : public G4VSensitiveDetector{
public:
	FlashSensitiveDetector_opt(G4String name);
	virtual ~FlashSensitiveDetector_opt();
	virtual void Initialize(G4HCofThisEvent*HCE);
	virtual G4bool ProcessHits(G4Step* step,G4TouchableHistory* ROhist);
	virtual void EndOfEvent(G4HCofThisEvent*HCE);
	private:FlashHitsCollection * hitsCollection;
	int collectionID;};
	
#endif
