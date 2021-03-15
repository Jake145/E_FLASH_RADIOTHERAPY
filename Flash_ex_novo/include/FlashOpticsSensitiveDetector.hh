#include "G4VSensitiveDetector.hh"
class FlashSensitiveDetector : public G4VSensitiveDetector{
public:
	FlashSensitiveDetector(G4String name,G4bool Kinetic_or_Optic);
	virtual ~FlashSensitiveDetector();
	G4bool Kin_Opt;
	virtual void Initialize(G4HCofThisEvent*HCE);
	virtual G4bool ProcessHits(G4Step* step,G4TouchableHistory* ROhist);
	virtual void EndOfEvent(G4HCofThisEvent*HCE);
	private:FlashHitsCollection * hitsCollection;
	int collectionID;}
