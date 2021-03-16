#ifndef FSensitiveDetector_h
#define FSensitiveDetector_h 1
#include "G4VSensitiveDetector.hh"
#include "FlashHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

/// Flash sensitive detector

class FSensitiveDetector : public G4VSensitiveDetector
{
public:
    FSensitiveDetector(G4String name);
    virtual ~FSensitiveDetector();
    
    virtual void Initialize(G4HCofThisEvent*HCE);
    virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    
private:
    FlashHitsCollection* HitsCollection;
    G4int collectionID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
