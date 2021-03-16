#ifndef FlashHit_h
#define FlashHit_h 1

 

#include "G4VHit.hh"

#include "G4THitsCollection.hh"

#include "G4Allocator.hh"

#include "G4ThreeVector.hh"

#include "G4ParticleDefinition.hh"


 

class FlashHit : public G4VHit

{

  public:

 

      FlashHit();

      ~FlashHit();

      G4int operator==(const FlashHit &right) const;

 

      inline void *operator new(size_t);

      inline void operator delete(void *aHit);

      void Draw();

      void Print();

 

  private:
	G4String    Stat;
      G4String      Pname;
      
      G4int         ChK;
      
      G4int 	    Scint;
	
      G4int         stripNo;

      G4ThreeVector position;

      G4ThreeVector momentum;

      G4double      energy;
      
      G4double 	    Edep;

      G4ParticleDefinition* particle;

 

  public:
      inline void SetStatus(G4String Status)
      {Stat=Status;}
      
      inline G4String GetStatus_Kin()
      {return Stat;}
      inline void SetCherenkovCount(G4int CherenkovCount)
     
      { ChK=CherenkovCount; }
  	
      inline  G4int GetCherenkovCount()
     
      { return ChK; }
  	
      inline void SetScintilCount(G4int ScintilCount)
     
      { Scint=ScintilCount; }
  	
      inline  G4int GetScintilCount()
     
      { return Scint; }
  	  	
      inline void SetProcess(G4String ProcessName)
  	
      { Pname=ProcessName; }
  	
      inline G4String GetProcessName()
	
      { return Pname; }
	
      inline void SetStripNo(G4int strip)

      { stripNo=strip; }

      inline G4int GetStripNo()

      { return stripNo; }

      inline void SetPosition(G4ThreeVector pos)

      { position=pos; }

      inline G4ThreeVector GetPosition()

      { return position; }

      inline void SetMomentum(G4ThreeVector mom)

      { momentum = mom; }

      inline G4ThreeVector GetMomentum()

      { return momentum; }

      inline void SetEnergy(G4double ene)

      { energy = ene; }

      inline G4double GetEnergy()

      { return energy; }

      inline void SetParticle(G4ParticleDefinition* pdef)

      { particle = pdef; }

      inline G4ParticleDefinition* GetParticle()

      { return particle; }
      
      inline void SetEdep(G4double edep)

      { Edep = edep; }

      inline G4double GetEdep()

      { return Edep; }

};

typedef G4THitsCollection<FlashHit> FlashHitsCollection;

extern G4ThreadLocal G4Allocator<FlashHit>* FlashHitAllocator;

inline void* FlashHit::operator new(size_t)
{
    if (!FlashHitAllocator)
        FlashHitAllocator = new G4Allocator<FlashHit>;
    return (void*)FlashHitAllocator->MallocSingle();
}

inline void FlashHit::operator delete(void*aHit)
{
    FlashHitAllocator->FreeSingle((FlashHit*) aHit);
}

#endif
