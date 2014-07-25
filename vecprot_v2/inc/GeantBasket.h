#ifndef GEANT_BASKET
#define GEANT_BASKET

#include "TObject.h"
#include "TGeoExtension.h"
#include "TGeoVolume.h"
#include "GeantTrack.h" 
#include "sync_objects.h" // includes <atomic.h>

//==============================================================================
// Basket of tracks in the same volume which are transported by a single thread
//==============================================================================

class TGeoVolume;
class GeantBasketMgr;

//______________________________________________________________________________
class GeantBasket : public TObject {
public:
  enum EBasketFlags {
    kMixed =  BIT(14)       // Basked mixing tracks from different volumes
  };
protected:
   GeantBasketMgr   *fManager;             // Manager for the basket
   GeantTrack_v      fTracksIn;            // Vector of input tracks
   GeantTrack_v      fTracksOut;           // Vector of output tracks
//   GeantHit_v        fHits;              // Vector of produced hits
#if __cplusplus >= 201103L
   std::atomic_int   fAddingOp;            // Number of track adding ops
#endif

private:
   GeantBasket(const GeantBasket&);     // Not implemented
   GeantBasket& operator=(const GeantBasket&);  // Not implemented
public:
   GeantBasket();
   GeantBasket(Int_t size, GeantBasketMgr *mgr);
   virtual ~GeantBasket();
   
   // Add track from generator or physics process
   void              AddTrack(GeantTrack &track);
   // Add track from a track_v array (copied)
   void              AddTrack(GeantTrack_v &tracks, Int_t itr);
   // Add multiple tracks from a track_v array
   void              AddTracks(GeantTrack_v &tracks, Int_t istart, Int_t iend);
   virtual void      Clear(Option_t *option="");
   Bool_t            Contains(Int_t evstart, Int_t nevents=1) const;
   Int_t             GetNinput() const {return fTracksIn.GetNtracks();}
   Int_t             GetNoutput() const {return fTracksOut.GetNtracks();}
   GeantTrack_v     &GetInputTracks() {return fTracksIn;}
   GeantTrack_v     &GetOutputTracks() {return fTracksOut;}
   GeantBasketMgr   *GetBasketMgr() const {return fManager;}
   TGeoVolume       *GetVolume() const;
   Bool_t            IsMixed() const {return TObject::TestBit(kMixed);}
   inline Bool_t     IsAddingOp() const           {return (fAddingOp.load());}
   inline Int_t      LockAddingOp()               {return ++fAddingOp;}
   inline Int_t      UnLockAddingOp()             {return --fAddingOp;}
   virtual void      Print(Option_t *option="") const;
   void              PrintTrack(Int_t itr, Bool_t input=kTRUE) const;
   void              Recycle();
   void              SetMixed(Bool_t flag) {TObject::SetBit(kMixed, flag);}
   
   ClassDef(GeantBasket,1)  // A basket containing tracks in the same geomety volume
};

//______________________________________________________________________________
// Basket manager for a given volume. Holds a list of free baskets stored in a
// concurrent queue
//______________________________________________________________________________

class GeantScheduler;

//______________________________________________________________________________
class GeantBasketMgr : public TGeoExtension {
protected:
   GeantScheduler   *fScheduler;             // Scheduler for this basket
   TGeoVolume       *fVolume;                // Volume for which applies
   Int_t             fNumber;                // Number assigned
   Int_t             fBcap;                  // max capacity of baskets held
#if __cplusplus >= 201103L
   std::atomic_int   fThreshold;             // Adjustable transportability threshold
   std::atomic_int   fNbaskets;              // Number of baskets for this volume
   std::atomic_int   fNused;                 // Number of baskets in use
   typedef std::atomic<GeantBasket*> atomic_basket;
   atomic_basket     fCBasket;               // Current basket being filled
   atomic_basket     fPBasket;               // Current priority basket being filled
   std::atomic_flag  fLock;                  // atomic lock
#endif
   dcqueue<GeantBasket> fBaskets;            // queue of available baskets
   dcqueue<GeantBasket> *fFeeder;            // feeder queue to which baskets get injected
   TMutex            fMutex;                 // Mutex for this basket manager
private:
   GeantBasketMgr(const GeantBasketMgr&);    // Not implemented
   GeantBasketMgr& operator=(const GeantBasketMgr&); // Not implemented
#if __cplusplus >= 201103L
   GeantBasket      *StealAndReplace(atomic_basket &current);
   GeantBasket      *StealAndPin(atomic_basket &current);
   Bool_t            StealMatching(atomic_basket &global, GeantBasket *content);
#endif
public:

   GeantBasket      *GetNextBasket();

public:
   GeantBasketMgr() : fScheduler(0), fVolume(0), fNumber(0), fBcap(0), fThreshold(0), fNbaskets(0), 
                         fNused(0), fCBasket(0), fPBasket(0), fBaskets(), fLock(), fFeeder(0), fMutex() {}
   GeantBasketMgr(GeantScheduler *sch, TGeoVolume *vol, Int_t number);
   virtual ~GeantBasketMgr();
   
   virtual TGeoExtension *Grab() {return this;}
   virtual void           Release() const {delete this;}
   Int_t             AddTrack(GeantTrack &track, Bool_t priority=kFALSE);
   Int_t             AddTrack(GeantTrack_v &trackv, Int_t itr, Bool_t priority=kFALSE);
   Int_t             CollectPrioritizedTracks(Int_t evmin, Int_t evmax);
   Int_t             FlushPriorityBasket();
   Int_t             GarbageCollect();
   void              SetBcap(Int_t capacity)      {fBcap = capacity;}
   Int_t             GetBcap() const              {return fBcap;}
#if __cplusplus >= 201103L
   Int_t             GetNbaskets() const          {return fNbaskets.load();}
   Int_t             GetNused() const             {return fNused.load();}
   Int_t             GetThreshold() const         {return fThreshold.load();}
   void              SetThreshold(Int_t thr)      {fThreshold.store(thr);}
   GeantBasket      *GetCBasket() const           {return fCBasket.load();}
   GeantBasket      *GetPBasket() const           {return fPBasket.load();}
   void             *SetCBasket(GeantBasket* basket) {fCBasket.store(basket);}
   void             *SetPBasket(GeantBasket* basket) {fPBasket.store(basket);}
#endif   
   GeantScheduler   *GetScheduler() const         {return fScheduler;}
   const char       *GetName() const              {return (fVolume)?fVolume->GetName():ClassName();}
   Int_t             GetNumber() const            {return fNumber;}
   TGeoVolume       *GetVolume() const            {return fVolume;}
   virtual void      Print(Option_t *option="") const;
   void              RecycleBasket(GeantBasket *b);
   void              SetFeederQueue(dcqueue<GeantBasket> *queue) {fFeeder = queue;}
   dcqueue<GeantBasket> *GetFeederQueue() const {return fFeeder;}
   
   ClassDef(GeantBasketMgr,0)  // A path in geometry represented by the array of indices
};   
#endif
