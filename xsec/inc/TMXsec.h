// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMXsec
#define ROOT_TMXsec


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMXSec                                                               //
//                                                                      //
// X-section for G5 per material                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Geant/Config.h"

#include <TEXsec.h>
#include <vector>

class TPDecay;
class GeantTrack_v;
class GeantTaskData;

class TMXsec {

public:
   TMXsec();
   TMXsec(const Char_t* name, const Char_t *title,
	  const Int_t z[], const Int_t a[], const Float_t w[], 
	  Int_t nel, Float_t dens, Bool_t weight=kFALSE, const TPDecay *decaytable=0);
   virtual ~TMXsec();
   const Char_t* GetName() const {return fName;}
   const Char_t* GetTitle() const {return fTitle;}
   Float_t Xlength(Int_t part, Float_t en, Double_t ptot);
//   Bool_t Xlength_v(Int_t npart, const Int_t part[], const Float_t en[], Double_t lam[]);
   Float_t DEdx(Int_t part, Float_t en);
//   Bool_t DEdx_v(Int_t npart, const Int_t part[], const Float_t en[], Float_t de[]);
   Float_t Range(Int_t part, Float_t en);
   Double_t InvRange(Int_t part, Float_t step);

   GEANT_CUDA_DEVICE_CODE
   void  Eloss(Int_t ntracks, GeantTrack_v &tracks);
   GEANT_CUDA_DEVICE_CODE
   void  ElossSingle(Int_t itrack, GeantTrack_v &tracks);
   void	ProposeStep(Int_t ntracks, GeantTrack_v &tracks, GeantTaskData *td);
   void	ProposeStepSingle(Int_t itr, GeantTrack_v &tracks, GeantTaskData *td);
   void	SampleInt(Int_t ntracks, GeantTrack_v &tracksin, GeantTaskData *td);
   void	SampleSingleInt(Int_t itr, GeantTrack_v &tracksin, GeantTaskData *td);
   GEANT_CUDA_DEVICE_CODE
   Float_t MS(Int_t ipart, Float_t energy);

   TEXsec *SampleInt(Int_t part, Double_t en, Int_t &reac, Double_t ptotal);
   Int_t SampleElement(GeantTaskData *td); // based on # atoms/vol. for the prototype
   Int_t SampleElement(); // based on # atoms/vol. for Geant4 with tab.phys. 

   Int_t SelectElement(Int_t pindex, Int_t rindex, Double_t energy);

//   static Bool_t Prune();
   void Print(Option_t * opt="") const;

private:
   TMXsec(const TMXsec&);      // Not implemented
   TMXsec& operator=(const TMXsec&);      // Not implemented

   char            fName[32];  // cross section name
   char            fTitle[128];// cross section title
   Int_t           fNEbins;    // number of energy bins
   Int_t           fNTotXL;    // dimension of fTotXL
   Int_t           fNCharge;   // dimension of tables for charged particles
   Int_t           fNRelXS;    // dimension of fRelXS
   Double_t        fEilDelta;  // logarithmic energy delta
   const Double_t *fEGrid;     //! Energy grid

   Int_t           fNElems;    // Number of elements
   TEXsec        **fElems;     // [fNElems] List of elements composing this material
   Float_t        *fTotXL;     // [fNTotXL] Total x-sec for this material
   Float_t        *fRelXS;     // [fNRelXS] Relative x-sec for this material
   Float_t        *fDEdx;      // [fNCharge] Ionisation energy loss for this material
   Float_t        *fMSangle;   // [fNCharge] table of MS average angle
   Float_t        *fMSansig;   // [fNCharge] table of MS sigma angle
   Float_t        *fMSlength;  // [fNCharge] table of MS average lenght correction
   Float_t        *fMSlensig;  // [fNCharge] table of MS sigma lenght correction
   Double_t       *fRatios;    // [fNElems]  relative #atoms/volume; normalized
   Float_t        *fRange;     // [fNCharge] ranges of the particle in this material
   std::vector< std::pair<Float_t,Double_t> > **fInvRangeTable; // [fNCharge] 
   const TPDecay  *fDecayTable;// pointer to the decay table

   ClassDef(TMXsec,1)  //Material X-secs

};

class TOMXsec : public TObject {
public:
 TOMXsec(TMXsec* mxsec): fMXsec(mxsec) {}
 ~TOMXsec() {delete fMXsec;}
 TMXsec* MXsec() const {return fMXsec;}

private:
   TOMXsec(const TOMXsec&); // Not implemented
   TOMXsec& operator=(const TOMXsec&); // Not implemented

   TMXsec *fMXsec;
};

#endif