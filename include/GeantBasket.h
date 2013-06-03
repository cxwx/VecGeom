#ifndef GEANT_BASKET
#define GEANT_BASKET

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "TMutex.h"
#include "GeantVolumeBasket.h"
 
//==============================================================================
// Basket of tracks in the same volume which are transported by a single thread
//==============================================================================

class TGeoVolume;

//______________________________________________________________________________
class GeantBasket : public TObject {
protected:
   Int_t             fNtracks;               // Number of tracks
   Int_t             fMaxTracks;             // Max number of tracks
   Int_t            *fIndex;                 //[fNtracks] Track indices in the global stack
   // have to change it to hold local tracks

public:
   GeantBasket();
   GeantBasket(Int_t maxtracks);
   virtual ~GeantBasket();
   
   void              AddTrack(Int_t itrack);
   void              AddTracks(const Int_t *array, Int_t ntracks);
   virtual void      Clear(Option_t *option="");
   Bool_t            Contains(Int_t event) const;
   Int_t             GetNtracks() const {return fNtracks;}
   Int_t            *GetTracks() const {return fIndex;}
   virtual void      Print(Option_t *option="") const;
   void              Resize(Int_t newSize);
   
   ClassDef(GeantBasket,1)  // A basket containing tracks in the same geomety volume
};   

//==============================================================================
// GeantMainScheduler - dispatcher running in a single thread. Collects tracks
// from all threads via an input queue and fills baskets corresponding to each
// volume, which are then injected in the main work queue.
//==============================================================================

class concurrent_queue;

//______________________________________________________________________________
class GeantMainScheduler : public TObject {
protected:
   Int_t                fNvolumes;            // Number of active volumes in the geometry
   Int_t                fNpriority;           // Number of priority baskets held
   GeantBasket        **fBaskets;             // Array of baskets to be filled
   GeantBasket        **fPriorityBaskets;     // Array of priority baskets
   Int_t                fPriorityRange[2];    // Prioritized events
   concurrent_queue    *feeder_queue;         // Main feeder
   concurrent_queue    *empty_queue;          // Queue with empty baskets
   concurrent_queue    *collector_queue;      // Queue for collecting tracks
   
public:
   GeantMainScheduler();
   GeantMainScheduler(Int_t nvolumes);
   virtual ~GeantMainScheduler();
   Int_t                AddTrack(Int_t itrack, Int_t ibasket);
   Int_t                GetNpriority() const {return fNpriority;}
   void                 SetPriorityRange(Int_t min, Int_t max) {fPriorityRange[0]=min; fPriorityRange[1]=max;}
   Int_t                FlushBaskets(Int_t threshold=0);   
   Int_t                FlushPriorityBaskets();
   
   ClassDef(GeantMainScheduler, 1)      // Main basket scheduler
};   

//==============================================================================
// GeantTrackCollector - A class collecting all tracks propagated by a single
// thread across the boundaries of the current volume.
//==============================================================================

class GeantVolumeBasket;

//______________________________________________________________________________
class GeantTrackCollection : public TObject {
protected:
   Int_t               fNtracks;              // Number of crossing tracks
   Int_t               fSize;                 // Size of the arrays > fNtracks
   Int_t              *fTracks;               //! Indices of crossing tracks
   GeantVolumeBasket **fBaskets;              //! Volume basket to which the track goes

public:
   GeantTrackCollection();
   GeantTrackCollection(Int_t size);
   virtual ~GeantTrackCollection();
   
   Int_t               AddTrack(Int_t itrack, GeantVolumeBasket *basket);
   Int_t               GetNtracks() const {return fNtracks;}
   Int_t               FlushTracks(GeantMainScheduler *main);
   virtual void        Print(Option_t *option="") const;
   
   ClassDef(GeantTrackCollection, 1)     // Track collection per thread
};

#endif
