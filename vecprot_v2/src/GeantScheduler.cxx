#include "GeantScheduler.h"

#include "TMath.h"
#include "GeantBasket.h"
#include "globals.h"
#include "GeantTrack.h"
#include "GeantTaskData.h"
#include "GeantEvent.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/SimpleNavigator.h"
#include "base/Vector3D.h"
#include "management/GeoManager.h"
#endif
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoManager.h"

ClassImp(GeantScheduler)

//______________________________________________________________________________
GeantScheduler::GeantScheduler()
  : TObject(), fNvolumes(0), fNpriority(0), fBasketMgr(0), fGarbageCollector(0),
    fNstvol(0), fIstvol(0), fNvect(0), fNsteps(0), fCrtMgr(0), fCollecting(false), fLearning(ATOMIC_FLAG_INIT), 
    fGBCLock(ATOMIC_FLAG_INIT) {
  // Default constructor
  fPriorityRange[0] = fPriorityRange[1] = -1;
  SetLearning(false);
  fNvect = new Int_t[257];
  memset(fNvect, 0, 257*sizeof(int));
#ifdef __STAT_DEBUG
  fPStat.InitArrays(GeantPropagator::Instance()->fNevents);
  fQStat.InitArrays(GeantPropagator::Instance()->fNevents);
  fTStat.InitArrays(GeantPropagator::Instance()->fNevents);
#endif
}

//______________________________________________________________________________
GeantScheduler::~GeantScheduler() {
  // dtor.
  delete fGarbageCollector;
  if (fBasketMgr) {
    for (Int_t ib = 0; ib < fNvolumes; ib++) {
      fBasketMgr[ib]->GetVolume()->SetFWExtension(0);
      delete fBasketMgr[ib];
    }
  }
  delete[] fBasketMgr;
  delete[] fNstvol;
  delete[] fIstvol;
  delete[] fNvect;
}

//______________________________________________________________________________
void GeantScheduler::ActivateBasketManagers()
{
// Activate basket managers based on the distribution of steps in corresponding
// volumes.
  TMath::Sort(fNvolumes, fNstvol, fIstvol);
  int ntot = 0;
  int nsum = 0;
  // Percent of steps to activate basketized volumes
  double threshold = 0.9;
  int nactive = 0;
  TObjArray *vlist = gGeoManager->GetListOfVolumes();
  TGeoVolume *vol;
  for (auto i=0; i<fNvolumes; i++) {
    ntot += fNstvol[i];
    vol = (TGeoVolume*)vlist->At(i);
    GeantBasketMgr *mgr = (GeantBasketMgr*)vol->GetFWExtension();
    if (mgr->IsActive()) {
      nsum += fNstvol[i];
      nactive++;
    }
  }  
  int nthreshold = ntot*threshold;
  for (auto i=0; i<fNvolumes; ++i) {
    vol = (TGeoVolume*)vlist->At(fIstvol[i]);
    GeantBasketMgr *mgr = (GeantBasketMgr*)vol->GetFWExtension();
    if (!mgr->IsActive()) {
      mgr->Activate();
      nsum += fNstvol[fIstvol[i]];
      nactive++;
    }  
    if (nsum > nthreshold) break;
  }
  threshold = double(nsum)/ntot;
  Printf("Activated %d volumes accounting for %4.1f%% of track steps", nactive, 100*threshold);
  int nprint = 10;
  if (nprint > fNvolumes) nprint = fNvolumes;
  for (auto i=0; i<nprint; ++i) {
    vol = (TGeoVolume*)vlist->At(fIstvol[i]);
    Printf("  * %s: %d steps", vol->GetName(), fNstvol[fIstvol[i]]);
  }
}

//______________________________________________________________________________
void GeantScheduler::AdjustBasketSize() {
  // Adjust the basket size to converge to Ntracks/2*Nthreads
  return;     // !!!!!!!!!!!!!!! NOT needed anymore !!!!!!!!!!!!!!!
/*
  const Int_t min_size = 4;
  const Int_t max_size = gPropagator->fNperBasket;
  Int_t nthreads = gPropagator->fNthreads;
  Int_t nproposed;
  for (Int_t ib = 0; ib < fNvolumes; ib++) {
    if (!GetNtracks(ib))
      continue;
    nproposed = GetNtracks(ib) / (2 * nthreads);
    nproposed -= nproposed % 4;
    if (!nproposed)
      continue;
    if (nproposed < min_size)
      nproposed = min_size;
    else if (nproposed > max_size)
      nproposed = max_size;
    //      Printf("basket %s (%d tracks): crt=%d  proposed=%d",
    //       fBasketMgr[ib]->GetName(), fNtracks[ib], fBasketMgr[ib]->GetThreshold(), nproposed);
    fBasketMgr[ib]->SetThreshold(nproposed);
  }
*/
}

//______________________________________________________________________________
void GeantScheduler::CreateBaskets() {
  // Create the array of baskets
  if (fBasketMgr)
    return;
  fNvolumes = gGeoManager->GetListOfVolumes()->GetEntries();
  fBasketMgr = new GeantBasketMgr *[fNvolumes];
  fNstvol = new Int_t[fNvolumes];
  fIstvol = new Int_t[fNvolumes];
  memset(fNstvol, 0, fNvolumes*sizeof(Int_t));
  memset(fIstvol, 0, fNvolumes*sizeof(Int_t));
  Geant::priority_queue<GeantBasket *> *feeder = WorkloadManager::Instance()->FeederQueue();
  TIter next(gGeoManager->GetListOfVolumes());
  TGeoVolume *vol;
  GeantBasketMgr *basket_mgr;
  Int_t icrt = 0;
  Int_t nperbasket = gPropagator->fNperBasket;
  while ((vol = (TGeoVolume *)next())) {
    basket_mgr = new GeantBasketMgr(this, vol, icrt);
    basket_mgr->SetThreshold(nperbasket);
    vol->SetFWExtension(basket_mgr);
    basket_mgr->SetFeederQueue(feeder);
    //      Printf("basket %s: %p feeder=%p", basket_mgr->GetName(), basket_mgr,
    //      basket_mgr->GetFeederQueue());
    fBasketMgr[icrt++] = basket_mgr;
  }
  const size_t MB = 1048576;
  size_t size = Sizeof() / MB;
  Printf("Size of scheduler including initial baskets: %ld MB", size);
  //   PrintSize();
}

//______________________________________________________________________________
Int_t GeantScheduler::AddTrack(GeantTrack &track, GeantTaskData *td) {
  // Main method to inject generated tracks. Track status is kNew here.
  TGeoVolume *vol = track.fPath->GetCurrentNode()->GetVolume();
  GeantBasketMgr *basket_mgr = static_cast<GeantBasketMgr *>(vol->GetFWExtension());
  Int_t ivol = basket_mgr->GetNumber();
  fNstvol[ivol]++;
  fNsteps++;
  // If no learning phase requested, activate all basket managers
  GeantPropagator *propagator = GeantPropagator::Instance();
  if ((propagator->fLearnSteps==0) &&
        !fLearning.test_and_set(std::memory_order_acquire)) {
    for (ivol=0; ivol < fNvolumes; ++ivol) fBasketMgr[ivol]->Activate();
  } else {
    // Activate the basket manager
    basket_mgr->Activate(); 
  }    	
  return basket_mgr->AddTrack(track, false, td);
}

//______________________________________________________________________________
Int_t GeantScheduler::AddTracks(GeantBasket *output, Int_t &ntot, Int_t &nnew, Int_t &nkilled, GeantTaskData *td) {
  // Main re-basketizing method. Add all tracks and inject baskets if above threshold.
  // Returns the number of injected baskets.
  Int_t ninjected = 0;
  Bool_t priority = kFALSE;
  GeantPropagator *propagator = GeantPropagator::Instance();
  GeantTrack_v &tracks = output->GetOutputTracks();
  Int_t ntracks = tracks.GetNtracks();
  ntot += ntracks;
  GeantBasketMgr *basket_mgr = 0;
  TGeoVolume *vol = 0;
  for (Int_t itr = 0; itr < ntracks; itr++) {
    // We have to collect the killed tracks
    if (tracks.fStatusV[itr] == kKilled || tracks.fStatusV[itr] == kExitingSetup ||
        tracks.fPathV[itr]->IsOutside()) {
      nkilled++;
      gPropagator->StopTrack(tracks, itr);
      tracks.DeleteTrack(itr);
      continue;
    }
    if (tracks.fStatusV[itr] == kNew) nnew++;
    tracks.fStatusV[itr] = kAlive;
    
    priority = kFALSE;
//    if (fPriorityRange[0] >= 0 && tracks.fEventV[itr] >= fPriorityRange[0] &&
//        tracks.fEventV[itr] <= fPriorityRange[1])
//      priority = kTRUE;

    vol = tracks.fPathV[itr]->GetCurrentNode()->GetVolume();
    basket_mgr = static_cast<GeantBasketMgr *>(vol->GetFWExtension());
    Int_t ivol = basket_mgr->GetNumber();
    tracks.fVindexV[itr] = ivol;
    fNstvol[ivol]++;
    Int_t nsteps = ++fNsteps;
    // Detect if the event the track is coming from is prioritized
    if (propagator->fEvents[tracks.fEvslotV[itr]]->IsPrioritized()) {
      ninjected += td->fBmgr->AddTrackSingleThread(tracks, itr, true, td);
      continue;
    }
    if (propagator->fLearnSteps && (nsteps%propagator->fLearnSteps)==0 &&
        !fLearning.test_and_set(std::memory_order_acquire)) {
      Printf("=== Learning phase of %d steps completed ===", propagator->fLearnSteps);
      // Here comes the algorithm activating basket managers...
      ActivateBasketManagers();
      propagator->fLearnSteps *= 4;
      fLearning.clear();
    }  
    if (basket_mgr->IsActive()) 
      ninjected += basket_mgr->AddTrack(tracks, itr, priority, td);
    else 
      ninjected += td->fBmgr->AddTrackSingleThread(tracks, itr, false, td);
  }
  tracks.Clear();
  return ninjected;
}

//______________________________________________________________________________
Int_t GeantScheduler::CollectPrioritizedPerThread(GeantBasketMgr *collector, GeantTaskData *td) {
  // Collect all tracks from prioritized events. Called concurrently by worker
  // threads. The thread getting to process the last basket manager resets the
  // collection process. The method performs work steal.
  Int_t ncollected = 0;
  Int_t imgr;
  // The IsCollecting flag can only be set by the scheduler thread, but can be
  // reset by any worker
  while (IsCollecting()) {
    imgr = ++fCrtMgr;
    if (imgr >= fNvolumes) return ncollected;
    if (imgr == fNvolumes-1) {
      SetCollecting(false);
      fCrtMgr.store(0);
    }  
    // Process current basket manager
    ncollected += fBasketMgr[imgr]->CollectPrioritizedTracksNew(collector, td);
  }
  return ncollected;
}

//______________________________________________________________________________
Int_t GeantScheduler::CollectPrioritizedTracks(GeantTaskData *td) {
  // Send the signal to all basket managers to prioritize all pending tracks
  // if any within the priority range.
  //   PrintSize();
  Int_t ninjected = 0;
  for (Int_t ibasket = 0; ibasket < fNvolumes; ibasket++)
    ninjected +=
        fBasketMgr[ibasket]->CollectPrioritizedTracks(fPriorityRange[0], fPriorityRange[1], td);
  return ninjected;
}

//______________________________________________________________________________
Int_t GeantScheduler::FlushPriorityBaskets() {
  // Flush all non-empty priority baskets.
  Int_t ninjected = 0;
  for (Int_t ibasket = 0; ibasket < fNvolumes; ibasket++) {
    ninjected += fBasketMgr[ibasket]->FlushPriorityBasket();
  }
  return ninjected;
}

//______________________________________________________________________________
Int_t GeantScheduler::GarbageCollect(GeantTaskData *td, bool force) {
  // Flush all filled baskets in the work queue.
  //   PrintSize();
  // Only one thread at a time
  if (force) {
    while (fGBCLock.test_and_set(std::memory_order_acquire))
      ;
  } else {
    if (fGBCLock.test_and_set(std::memory_order_acquire)) return 0;
  }
//  Printf("=== Garbage collect");
  Int_t ninjected = 0;
  for (Int_t ibasket = 0; ibasket < fNvolumes; ibasket++) {
    if (fBasketMgr[ibasket]->IsActive())
      ninjected += fBasketMgr[ibasket]->GarbageCollect(td);
  }
  fGBCLock.clear(std::memory_order_release);
  return ninjected;
}

//______________________________________________________________________________
size_t GeantScheduler::Sizeof() const {
  // Returns size of the scheduler, including allocated baskets
  size_t size = sizeof(GeantScheduler);
  for (auto i = 0; i < fNvolumes; ++i)
    size += fBasketMgr[i]->Sizeof();
  if (fGarbageCollector)
    size += fGarbageCollector->Sizeof();
  return size;
}

//______________________________________________________________________________
void GeantScheduler::PrintSize() const {
  // Prints detailed breakdown of size allocated
  size_t size = Sizeof();
  Printf("Size of scheduler: %ld bytes", size);
  for (auto i = 0; i < fNvolumes; ++i)
    fBasketMgr[i]->PrintSize();
  if (fGarbageCollector)
    fGarbageCollector->PrintSize();
}
