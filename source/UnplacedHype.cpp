/// \file UnplacedHype.cpp
/// \author Marilena Bandieramonte (marilena.bandieramonte@cern.ch)

#include "volumes/UnplacedHype.h"

#include "management/VolumeFactory.h"
#include "volumes/SpecializedHype.h"

#include <stdio.h>
#include "base/RNG.h"


namespace VECGEOM_NAMESPACE {

    
    VECGEOM_CUDA_HEADER_BOTH
    void UnplacedHype::SetParameters(const Precision rMin, const Precision stIn,
                                     const Precision rMax, const Precision stOut,
                                     const Precision dz){
        
        //TODO: add eventual check
        fRmin=rMin;
        fStIn=stIn;
        fRmin=rMax;
        fStIn=stOut;
        fDz=dz;
    }
    
    VECGEOM_CUDA_HEADER_BOTH
    UnplacedHype::UnplacedHype(const Precision rMin, const Precision stIn,
                               const Precision rMax, const Precision stOut,
                               const Precision dz){
        
        SetParameters(rMin, stIn, rMax, stOut, dz);
        
        fTIn=tan(stIn*kDegToRad);       //Tangent of the Inner stereo angle
        fTOut=tan(stOut*kDegToRad);     //Tangent of the Outer stereo angle
        fTIn2=fTIn*fTIn;                //squared value of fTIn
        fTOut2=fTOut*fTOut;             //squared value of fTOut
        
        fRmin2=fRmin*fRmin;
        fRmax2=fRmax*fRmax;
        
    }


//__________________________________________________________________
    
    void UnplacedHype::Print() const {
        
    }
//__________________________________________________________________
    
    void UnplacedHype::Print(std::ostream &os) const {
        
    }

template <TranslationCode transCodeT, RotationCode rotCodeT>
VECGEOM_CUDA_HEADER_DEVICE
VPlacedVolume* UnplacedHype::Create(
    LogicalVolume const *const logical_volume,
    Transformation3D const *const transformation,
#ifdef VECGEOM_NVCC
    const int id,
#endif
    VPlacedVolume *const placement) {
  if (placement) {
    return new(placement) SpecializedHype<transCodeT, rotCodeT>(
#ifdef VECGEOM_NVCC
        logical_volume, transformation, NULL, id); // TODO: add bounding box?
#else
        logical_volume, transformation);
#endif
  }
  return new SpecializedHype<transCodeT, rotCodeT>(
#ifdef VECGEOM_NVCC
      logical_volume, transformation, NULL, id); // TODO: add bounding box?
#else
      logical_volume, transformation);
#endif
}

VECGEOM_CUDA_HEADER_DEVICE
VPlacedVolume* UnplacedHype::SpecializedVolume(
    LogicalVolume const *const volume,
    Transformation3D const *const transformation,
    const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
    const int id,
#endif
    VPlacedVolume *const placement) const {
  return VolumeFactory::CreateByTransformation<
      UnplacedHype>(volume, transformation, trans_code, rot_code,
#ifdef VECGEOM_NVCC
                              id,
#endif
                              placement);
}

} // End global namespace

namespace vecgeom {

#ifdef VECGEOM_CUDA_INTERFACE

void UnplacedHype_CopyToGpu(VUnplacedVolume *const gpu_ptr);

VUnplacedVolume* UnplacedHype::CopyToGpu(
    VUnplacedVolume *const gpu_ptr) const {
  UnplacedHype_CopyToGpu(gpu_ptr);
  CudaAssertError();
  return gpu_ptr;
}

VUnplacedVolume* UnplacedHype::CopyToGpu() const {
  VUnplacedVolume *const gpu_ptr = AllocateOnGpu<UnplacedHype>();
  return this->CopyToGpu(gpu_ptr);
}

#endif

#ifdef VECGEOM_NVCC

class VUnplacedVolume;

__global__
void UnplacedHype_ConstructOnGpu(VUnplacedVolume *const gpu_ptr) {
  new(gpu_ptr) vecgeom_cuda::UnplacedHype();
}

void UnplacedHype_CopyToGpu(VUnplacedVolume *const gpu_ptr) {
  UnplacedHype_ConstructOnGpu<<<1, 1>>>(gpu_ptr);
}

#endif

} // End namespace vecgeom
