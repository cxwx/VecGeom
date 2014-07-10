//===-- volumes/UnplacedHype.h - Instruction class definition -------*- C++ -*-===//
//
//                     GeantV - VecGeom
//
// This file is distributed under the LGPL
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
///
/// \file volumes/UnplacedHype.h
/// \author Marilena Bandieramonte (marilena.bandieramonte@cern.ch)
/// \brief This file contains the declaration of the UnplacedHype class
///
/// _____________________________________________________________________________
///Hyperboloid class is defined by 5 parameters
///A Hype is the solid bounded by the following surfaces:
/// - 2 planes parallel with XY cutting the Z axis at Z=-dz and Z=+dz
///- Inner and outer lateral surfaces. These represent the surfaces
///described by the revolution of 2 hyperbolas about the Z axis:
/// r^2 - (t*z)^2 = a^2 where:
/// r = distance between hyperbola and Z axis at coordinate z
/// t = tangent of the stereo angle (angle made by hyperbola asimptotic lines and Z axis). t=0 means cylindrical surface.
/// a = distance between hyperbola and Z axis at z=0
//===----------------------------------------------------------------------===//

#ifndef VECGEOM_VOLUMES_UNPLACEDHYPE_H_
#define VECGEOM_VOLUMES_UNPLACEDHYPE_H_

#include "base/Global.h"

#include "base/AlignedBase.h"
#include "volumes/UnplacedVolume.h"

namespace VECGEOM_NAMESPACE {

class UnplacedHype : public VUnplacedVolume, AlignedBase {

private:
    Precision fRmin;    //Inner radius
    Precision fStIn;    //Stereo angle for inner surface
    Precision fRmax;    //Outer radius
    Precision fStOut;   //Stereo angle for outer surface
    Precision fDz;      //z-coordinate of the cutting planes
    
    //Precomputed Values
    Precision fTIn;     //Tangent of the Inner stereo angle
    Precision fTOut;    //Tangent of the Outer stereo angle
    Precision fTIn2;    //Squared value of fTIn
    Precision fTOut2;   //Squared value of fTOut
    Precision fRmin2;   //Squared Inner radius
    Precision fRmax2;   //Squared Outer radius
    

    
public:
    
    //constructor
    VECGEOM_CUDA_HEADER_BOTH
    UnplacedHype(const Precision rmin, const Precision stIn, const Precision rmax, const Precision stOut, const Precision dz);
    
    //get
    VECGEOM_CUDA_HEADER_BOTH
    Precision GetRmin() const{ return fRmin;}

    VECGEOM_CUDA_HEADER_BOTH
    Precision GetRmax() const{ return fRmax;}
    
    VECGEOM_CUDA_HEADER_BOTH
    Precision GetRmin2() const{ return fRmin2;}
    
    VECGEOM_CUDA_HEADER_BOTH
    Precision GetRmax2() const{ return fRmax2;}

    VECGEOM_CUDA_HEADER_BOTH
    Precision GetStIn() const{ return fStIn;}

    VECGEOM_CUDA_HEADER_BOTH
    Precision GetStOut() const{ return fStOut;}

    VECGEOM_CUDA_HEADER_BOTH
    Precision GetTIn() const{ return fTIn;}

    VECGEOM_CUDA_HEADER_BOTH
    Precision GetTOut() const{ return fTOut;}

    VECGEOM_CUDA_HEADER_BOTH
    Precision GetTIn2() const{ return fTIn2;}
    
    VECGEOM_CUDA_HEADER_BOTH
    Precision GetTOut2() const{ return fTOut2;}
    
    VECGEOM_CUDA_HEADER_BOTH
    Precision GetDz() const{ return fDz;}
    
    //set
    VECGEOM_CUDA_HEADER_BOTH
    void SetParameters(const Precision rMin, const Precision stIn, const Precision rMax, const Precision stOut, const Precision dz);

    virtual int memory_size() const { return sizeof(*this); }

//__________________________________________________________________
    

    VECGEOM_CUDA_HEADER_BOTH
    virtual void Print() const;
//__________________________________________________________________
    

    virtual void Print(std::ostream &os) const;
    
//__________________________________________________________________
    


    template <TranslationCode transCodeT, RotationCode rotCodeT>
    VECGEOM_CUDA_HEADER_DEVICE
    static VPlacedVolume* Create(LogicalVolume const *const logical_volume,
                               Transformation3D const *const transformation,
#ifdef VECGEOM_NVCC
                               const int id,
#endif
                               VPlacedVolume *const placement = NULL);

#ifdef VECGEOM_CUDA_INTERFACE
  virtual VUnplacedVolume* CopyToGpu() const;
  virtual VUnplacedVolume* CopyToGpu(VUnplacedVolume *const gpu_ptr) const;
#endif

private:

  
    //Specialized Volume
  VECGEOM_CUDA_HEADER_DEVICE
  virtual VPlacedVolume* SpecializedVolume(
      LogicalVolume const *const volume,
      Transformation3D const *const transformation,
      const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
      const int id,
#endif
      VPlacedVolume *const placement = NULL) const;

};

} // End global namespace

#endif // VECGEOM_VOLUMES_UNPLACEDHYPE_H_