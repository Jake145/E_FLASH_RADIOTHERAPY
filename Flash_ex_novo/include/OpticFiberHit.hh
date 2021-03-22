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
/// \file optical/LXe/include/OpticFiberHit.hh
/// \brief Definition of the OpticFiberHit class
//
//
#ifndef OpticFiberHit_h
#define OpticFiberHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"

#include "tls.hh"

class G4VTouchable;

class OpticFiberHit : public G4VHit
{
  public:
 
    OpticFiberHit();
    virtual ~OpticFiberHit();
    OpticFiberHit(const OpticFiberHit &right);

    const OpticFiberHit& operator=(const OpticFiberHit &right);
    G4bool operator==(const OpticFiberHit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
 
    virtual void Draw();
    virtual void Print();

    

    inline void IncPhotonCount(){fPhotons++;}
    inline G4int GetPhotonCount(){return fPhotons;}
    
    inline void IncCerenkovCount(){fCerenkov++;}
    inline G4int GetCerenkovCount(){return fCerenkov;}
    
    inline void IncScintillationCount(){fScintillation++;}
    inline G4int GetScintillationCount(){return fScintillation;}
    


  

  private:


    G4int fPhotons;
        G4int fCerenkov;
            G4int fScintillation;


};

typedef G4THitsCollection<OpticFiberHit> OpticFiberHitsCollection;

extern G4ThreadLocal G4Allocator<OpticFiberHit>* OpticFiberHitAllocator;

inline void* OpticFiberHit::operator new(size_t){
  if(!OpticFiberHitAllocator)
      OpticFiberHitAllocator = new G4Allocator<OpticFiberHit>;
  return (void *) OpticFiberHitAllocator->MallocSingle();
}

inline void OpticFiberHit::operator delete(void *aHit){
  OpticFiberHitAllocator->FreeSingle((OpticFiberHit*) aHit);
}

#endif
