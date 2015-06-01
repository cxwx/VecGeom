// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEFstate
#define ROOT_TEFstate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TEXSec                                                               //
//                                                                      //
// X-section for G5 per material                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TPartIndex.h"

class TFile;

class TPFstate;
class TFinState;

class TEFstate {
public:
  
  TEFstate();
  TEFstate(Int_t z, Int_t a, Float_t dens);
  ~TEFstate();
  Bool_t AddPart(Int_t kpart, Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[]);
  Bool_t AddPart(Int_t kpart, Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[], TFinState vecfs[]);
  
  Bool_t AddPartFS(Int_t kpart, Int_t ibin, Int_t reac, const Int_t npart[], const Float_t weight[],
                   const Float_t kerma[], const Float_t en[], const Char_t surv[], const Int_t pid[],
                   const Float_t mom[]);
  
  Int_t Ele() const {return fEle;}
  Double_t Dens() const {return fDens;}
  Double_t Emin() const {return fEmin;}
  Double_t Emax() const {return fEmax;}
  Int_t NEbins() const {return fNEbins;}
  Double_t EilDelta() const {return fEilDelta;}
  Int_t NEFstat() const {return fNEFstat;}
  
  Int_t NRpart() const {return fNRpart;}
  
  void SetRestCaptFstate(Int_t kpart, const TFinState &fstate);
  Bool_t HasRestCapture(Int_t partindex);
 
  Bool_t SampleReac(Int_t pindex, Int_t preac, Float_t en, Int_t& npart, Float_t& weight,
                    Float_t& kerma, Float_t &enr, const Int_t *&pid, const Float_t *&mom,
                    Int_t& ebinindx) const;
  Bool_t SampleReac(Int_t pindex, Int_t preac, Float_t en, Int_t& npart, Float_t& weight,
                    Float_t& kerma, Float_t &enr, const Int_t *&pid, const Float_t *&mom, 
                    Int_t& ebinindx, Double_t randn1, Double_t randn2) const;
  Bool_t SampleRestCaptFstate(Int_t kpart,Int_t& npart, Float_t& weight,
                    Float_t& kerma, Float_t &enr, const Int_t *&pid, const Float_t *&mom) const;
  Bool_t SampleRestCaptFstate(Int_t kpart,Int_t& npart, Float_t& weight,
                    Float_t& kerma, Float_t &enr, const Int_t *&pid, const Float_t *&mom,
                    Double_t randn) const;
  Bool_t GetReac(Int_t pindex, Int_t preac, Float_t en, Int_t ifs, Int_t& npart, Float_t& weight,
                 Float_t& kerma, Float_t &enr, const Int_t *&pid, const Float_t *&mom) const;
  
  static Bool_t FloatDiff(Double_t a, Double_t b, Double_t prec) {
    return TMath::Abs(a-b)>0.5*TMath::Abs(a+b)*prec;
  }
  
  void Draw(Option_t *option);
  Bool_t Resample();
  
  Bool_t Prune();
  
  static Int_t NLdElems() {return fNLdElems;}
  
  static TEFstate *Element(Int_t i) {
    if(i<0 || i>=fNLdElems) return 0; return fElements[i];}
  
  static TEFstate *GetElement(Int_t z, Int_t a=0, TFile *f=0);
  static TEFstate **GetElements() {return fElements;}
  
private:
  TEFstate(const TEFstate &); // Not implemented
  TEFstate& operator=(const TEFstate &); // Not implemented
  
  Int_t          fEle;     // Element code Z*10000+A*10+metastable level
  Float_t        fDens;    // Density in g/cm3
  Double_t       fAtcm3;   // Atoms per cubic cm unit density
  Double_t       fEmin;    // Minimum of the energy Grid
  Double_t       fEmax;    // Maximum of the energy Grid
  Int_t          fNEbins;  // Number of log steps in energy
  Double_t       fEilDelta; // Inverse log energy step
  const Double_t *fEGrid;  //! Common energy grid
  Int_t          fNEFstat; // Number of sampled final states
  Int_t          fNRpart;  // Number of particles with reaction
  TPFstate      *fPFstate;   // [fNRpart] Final state table per particle
  
  static Int_t   fNLdElems; //! number of loaded elements
  static TEFstate *fElements[NELEM]; //! databases of elements
  
  ClassDefNV(TEFstate,1)  // Element X-secs
  
};


#endif