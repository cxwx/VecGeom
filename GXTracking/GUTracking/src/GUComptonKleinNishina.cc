#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "GUComptonKleinNishina.h"
#include <iostream>

#include "backend/Backend.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_HOST
GUComptonKleinNishina::GUComptonKleinNishina(Random_t* states, int threadId) 
  :
  fRandomState(states), fThreadId(threadId),
  fMinX(1.e-8),  fMaxX(1000.), // fDeltaX(0.1), 
  // fMinY(1.e-8),  fMaxY(1001.), fDeltaY(0.1),
  fMaxZelement(100),       // Elements up to Z=100
  fNrow(100), fNcol(100) 
{
  //replace hard coded numbers by default constants
  fAliasSampler = new GUAliasSampler(fRandomState, fThreadId,
                                     fMaxZelement,
                                     fMinX, fMaxX,
                                     fNrow, fNcol);

  for( int z= 1 ; z < fMaxZelement; ++z)
  {
    //eventually arguments of BuildTable should be replaced by members of *this
    //  and dropped from the function signature. Same for BuildPdfTable
    BuildOneTable(z, fMinX, fMaxX, fNrow, fNcol);
  }
}

VECPHYS_CUDA_HEADER_BOTH 
GUComptonKleinNishina::GUComptonKleinNishina(Random_t* states, int threadId,
                                             GUAliasSampler* sampler) 
  :
  fRandomState(states), fThreadId(threadId),
  fMinX(1.e-8),  fMaxX(1000.), // fDeltaX(0.1), 
  // fMinY(1.e-8),  fMaxY(1001.), fDeltaY(0.1),  
  fNrow(100), fNcol(100) 
{
  //replace hard coded numbers by default constants
  fAliasSampler = sampler;
}

//need another Ctor with setable parameters

VECPHYS_CUDA_HEADER_BOTH 
GUComptonKleinNishina::~GUComptonKleinNishina() 
{
  if(fAliasSampler) delete fAliasSampler;
}


VECPHYS_CUDA_HEADER_HOST void 
GUComptonKleinNishina::BuildOneTable( int Z, 
                                   const double xmin, 
                                   const double xmax,
                                   const int nrow,
                                   const int ncol)
{
  //for now, the model does not own pdf.  Otherwise, pdf should be the 
  //data member of *this and set its point to the fpdf of fAliasSampler 
  double *pdf = new double [(nrow+1)*ncol];

  BuildLogPdfTable(Z,xmin,xmax,nrow,ncol,pdf); 
  fAliasSampler->BuildAliasTable(Z,nrow,ncol,pdf);

  delete [] pdf;
}

VECPHYS_CUDA_HEADER_HOST void 
GUComptonKleinNishina::BuildPdfTable(int Z, 
                                     const double xmin, 
                                     const double xmax,
                                     const int nrow,
                                     const int ncol,
                                     double *p
                                     )
{
  // Build the probability density function (KleinNishina pdf) 
  // in the energy randge [xmin,xmax] with an equal bin size
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  //           xmin (miminum energy)
  //           xmax (maxinum energy)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  //build pdf  
  double dx = (xmax - xmin)/nrow;
  //  double xo =  xmin + 0.5*dx;

  for(int i = 0; i <= nrow ; ++i) {
    //for each input energy bin
    double x = dx*i;

    double ymin = x/(1+2.0*x*inv_electron_mass_c2);

    double dy = (x - ymin)/(ncol-1);
    double yo = ymin + 0.5*dy;
  
    double sum = 0.;

    for(int j = 0; j < ncol ; ++j) {
      //for each output energy bin
      double y = yo + dy*j;
      double xsec = CalculateDiffCrossSection(0,x,y);
      p[i*ncol+j] = xsec;
      sum += xsec;
    }

    //normalization
    sum = 1.0/sum;

    for(int j = 0; j < ncol ; ++j) {
      p[i*ncol+j] *= sum;
    }
  }
}

VECPHYS_CUDA_HEADER_HOST void 
GUComptonKleinNishina::BuildLogPdfTable(int Z, 
                                        const double xmin, 
                                        const double xmax,
                                        const int nrow,
                                        const int ncol,
                                        double *p)
{
  // Build the probability density function (KleinNishina pdf) 
  // in the energy randge [xmin,xmax] with an equal bin size
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  //           xmin (miminum energy)
  //           xmax (maxinum energy)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  //build pdf: logarithmic scale in the input energy bin  

  double logxmin = log(xmin);
  double dx = (log(xmax) - logxmin)/nrow;

  for(int i = 0; i <= nrow ; ++i) {
    //for each input energy bin
    double x = exp(logxmin + dx*i);

    double ymin = x/(1+2.0*x*inv_electron_mass_c2);
    double dy = (x - ymin)/(ncol-1);
    double yo = ymin + 0.5*dy;
  
    double sum = 0.;

    for(int j = 0; j < ncol ; ++j) {
      //for each output energy bin
      double y = yo + dy*j;
      double xsec = CalculateDiffCrossSection(Z,x,y);
      p[i*ncol+j] = xsec;
      sum += xsec;
    }

    //normalization
    sum = 1.0/sum;

    for(int j = 0; j < ncol ; ++j) {
      p[i*ncol+j] *= sum;
    }
  }
}


// function implementing the cross section for KleinNishina
// TODO: need to get electron properties from somewhere

VECPHYS_CUDA_HEADER_BOTH double 
GUComptonKleinNishina::CalculateDiffCrossSection( int Zelement, 
                                                  double energy0, 
                                                  double energy1 ) const
{
  // based on Geant4 : KleinNishinaCompton
  // input  : energy0 (incomming photon energy)
  //          energy1 (scattered photon energy)
  // output : dsigma  (differential cross section) 

  double E0_m = energy0*inv_electron_mass_c2 ;
  double epsilon = energy1/energy0;

  double onecost = (1.- epsilon)/(epsilon*E0_m);
  double sint2   = onecost*(2.-onecost);
  double greject = 1. - epsilon*sint2/(1.+ epsilon*epsilon);
  double dsigma = (epsilon + 1./epsilon)*greject;

  return dsigma;
}

} // end namespace impl
} // end namespace vecphys