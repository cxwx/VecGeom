#include "TGeoManager.h"
#include "TGeoBBox_v.h"
#include "TRandom.h"
#include "TMath.h"

#include <iostream>
#include "tbb/tick_count.h"" // timing from Intel TBB 
#include <cassert>

struct TStopWatch 
{
  tbb::tick_count t1;
  tbb::tick_count t2;
  void Start(){ t1=tbb::tick_count::now(); }
  void Stop(){ t2=tbb::tick_count::now(); }
  void Reset(){ /* */ ;}
  void Print(){  std::cerr << (t2 - t1).seconds() << std::endl; }
  double getDeltaSecs() { return (t2-t1).seconds(); }
};

#define NREP 10

main(int argc, char *argv[])
{
  int npoints=100;
  if(argc>1) sscanf(argv[1],"%d",&npoints);
  printf("npoints = %d\n",npoints);

  const Double_t dx=10; // these are half-distances
  const Double_t dy=20;
  const Double_t dz=30;

  TGeoManager *testvec = new TGeoManager("Test","This is a naive test");
  TGeoMaterial *vacmat = new TGeoMaterial("vacuum",0,0,0);
  TGeoMedium *vacmed = new TGeoMedium("vacuum",0,vacmat);

  TGeoVolume *world = testvec->MakeBox("world",vacmed,100,100,100);
  testvec->SetTopVolume(world);

  TGeoVolume *tbox = testvec->MakeBox("tbox",vacmed,dx,dy,dz);
  tbox->SetLineColor(kRed);
  tbox->SetFillColor(kRed);
  tbox->SetVisibility(1);
  world->AddNode(tbox,1,0);
  
  testvec->CloseGeometry();

  Double_t origin[3]={0,0,0};
  TGeoBBox_v *box = new TGeoBBox_v(dx, dy, dz,origin);
  const Double_t r3two = TMath::Power(2,1./3.);

  npoints=10;
  for(int i = 0 ;i < 14; i++) 
    {
      Double_t *points = new Double_t[3*npoints];
      TStopWatch tt;

      for(int i=0; i<npoints; ++i) {
	points[3*i  ]=r3two*(1-2.*gRandom->Rndm())*dx;
	points[3*i+1]=r3two*(1-2.*gRandom->Rndm())*dy;
	points[3*i+2]=r3two*(1-2.*gRandom->Rndm())*dz;
      }

      double *safety_v = new double[npoints];
      Bool_t *in_v = new Bool_t[npoints];
      double DeltaT=0., DeltaT_v=0.;
      for ( unsigned int repetitions = 0; repetitions < NREP; repetitions ++ ) 
	{
	  // assert correctness of result (simple checksum check)
	  {
	    double checksum=0., checksum_v=0.;
	    for(int i=0; i<npoints; ++i) {
	      in_v[i] = box->Contains(&points[3*i]); 
	      safety_v[i] = box->Safety(&points[3*i], in_v[i]);
	      checksum+=safety_v[i];
	    }
	    
	    
	    box->Contains_v(points,in_v,npoints);
	    box->Safety_v(points,safety_v,in_v,npoints);
	    for(int i=0; i<npoints; ++i) {
	      checksum_v+=safety_v[i];
	    }
	    //	    std::cerr << checksum_v << " vs " << checksum << std::endl;
	    //	    assert(checksum_v == checksum);
	  }


	  // measure timings here separately
	  tt.Start();
	  for(int i=0; i<npoints; ++i) {
	    in_v[i] = box->Contains(&points[3*i]); 
	    safety_v[i] = box->Safety(&points[3*i], in_v[i]);
	  }
	  tt.Stop();
	  DeltaT+= tt.getDeltaSecs();
	  tt.Reset();
	  
	  tt.Start();
	  box->Contains_v(points,in_v,npoints);
	  box->Safety_v(points,safety_v,in_v,npoints);
	  tt.Stop();
	  DeltaT_v+= tt.getDeltaSecs(); //      tt.Print();
	  tt.Reset();
	}

      std::cerr << npoints << " " << DeltaT/NREP << " " << DeltaT_v/NREP << " " << DeltaT/DeltaT_v << std::endl;
      
      delete[] safety_v;
      delete[] points;
      npoints*=2;
    }
  return 0;
}