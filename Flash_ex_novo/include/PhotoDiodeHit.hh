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
/// \file optical/LXe/include/PhotoDiodeHit.hh
/// \brief Definition of the PhotoDiodeHit class
//
//
#ifndef PhotoDiodeHit_h
#define PhotoDiodeHit_h 1

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

class PhotoDiodeHit : public G4VHit
{
  public:
 
    PhotoDiodeHit();
    virtual ~PhotoDiodeHit();
    PhotoDiodeHit(const PhotoDiodeHit &right);

    const PhotoDiodeHit& operator=(const PhotoDiodeHit &right);
    G4bool operator==(const PhotoDiodeHit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
 
    virtual void Draw();
    virtual void Print();

    

    inline void IncPhotonCount_photo(){fPhotons_photo++;}
    inline G4int GetPhotonCount_photo(){return fPhotons_photo;}
    
    inline void IncCerenkovCount_photo(){fCerenkov_photo++;}
    inline G4int GetCerenkovCount_photo(){return fCerenkov_photo;}
    
    inline void IncScintillationCount_photo(){fScintillation_photo++;}
    inline G4int GetScintillationCount_photo(){return fScintillation_photo;}
    


  

  private:


    G4int fPhotons_photo;
        G4int fCerenkov_photo;
            G4int fScintillation_photo;


};

typedef G4THitsCollection<PhotoDiodeHit> PhotoDiodeHitsCollection;

extern G4ThreadLocal G4Allocator<PhotoDiodeHit>* PhotoDiodeHitAllocator;

inline void* PhotoDiodeHit::operator new(size_t){
  if(!PhotoDiodeHitAllocator)
      PhotoDiodeHitAllocator = new G4Allocator<PhotoDiodeHit>;
  return (void *) PhotoDiodeHitAllocator->MallocSingle();
}

inline void PhotoDiodeHit::operator delete(void *aHit){
  PhotoDiodeHitAllocator->FreeSingle((PhotoDiodeHit*) aHit);
}

#endif
