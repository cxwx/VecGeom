#include "GeantThreadData.h"
#include "globals.h"
#include "GeantBasket.h"
#include "GeantPropagator.h"

#include "TGeoMatrix.h"
#include "TArrayI.h"
#include "TGeoVolume.h"
#include "TRandom.h"
#include "TGeoHelix.h"

ClassImp(GeantThreadData)

//______________________________________________________________________________
GeantThreadData::GeantThreadData()
            :TObject(),
             fMaxPerBasket(0),
             fNprocesses(0),
             fMatrix(0),
             fVolume(0),
             fRndm(0),
             fDblArray(0),
             fProcStep(0),
             fPartInd(0),
             fPartNext(0),
             fPartTodo(0),
             fPartCross(0),
             fFieldPropagator(0),
             fRotation(0),
             fCollection(0)
{
// I/O ctor.
}

//______________________________________________________________________________
GeantThreadData::GeantThreadData(Int_t maxperbasket, Int_t maxprocesses)
            :TObject(),
             fMaxPerBasket(maxperbasket),
             fNprocesses(maxprocesses),
             fMatrix(new TGeoHMatrix()),
             fVolume(0),
             fRndm(new TRandom()),
             fDblArray(0),
             fProcStep(0),
             fPartInd(0),
             fPartNext(0),
             fPartTodo(0),
             fPartCross(0),
             fFieldPropagator(0),
             fRotation(new TGeoRotation()),
             fCollection(0)
{
// Constructor
   GeantPropagator *propagator = GeantPropagator::Instance();
   fDblArray  = new Double_t[5*fMaxPerBasket];
   fProcStep  = new Double_t[fNprocesses*fMaxPerBasket];
   fPartInd   = new TArrayI(fMaxPerBasket);
   fPartNext  = new TArrayI(fMaxPerBasket);
   fPartTodo  = new TArrayI(fMaxPerBasket);
   fPartCross = new TArrayI(fMaxPerBasket);
   fFieldPropagator = new TGeoHelix(1,1);
   fFieldPropagator->SetField(0,0,propagator->fBmag, kFALSE);
   fCollection = new GeantTrackCollection(100);
}

//______________________________________________________________________________
GeantThreadData::~GeantThreadData()
{
// Destructor
   delete fMatrix;
   delete fRndm;
   delete [] fDblArray;
   delete [] fProcStep;
   delete fPartInd;
   delete fPartNext;
   delete fPartTodo;
   delete fPartCross;
   delete fFieldPropagator;
   delete fRotation;
   delete fCollection;
}
