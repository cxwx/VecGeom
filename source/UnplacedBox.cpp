/// \file UnplacedBox.cpp
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#include "volumes/UnplacedBox.h"

#include "backend/Backend.h"
#include "management/VolumeFactory.h"
#include "volumes/SpecializedBox.h"
#include "base/RNG.h"
#include <stdio.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECGEOM_CUDA_HEADER_BOTH
void UnplacedBox::Print() const {
  printf("UnplacedBox {%.2f, %.2f, %.2f}", x(), y(), z());
}

void UnplacedBox::Print(std::ostream &os) const {
  os << "UnplacedBox {" << x() << ", " << y() << ", " << z() << "}";
}


//______________________________________________________________________________
void UnplacedBox::Extent(Vector3D<Precision> & aMin, Vector3D<Precision> & aMax) const
{
    // Returns the full 3D cartesian extent of the solid.
  aMin.x() = -dimensions_[0];
  aMax.x() = dimensions_[0];
  aMin.y() = -dimensions_[1];
  aMax.y() = dimensions_[1];
  aMin.z() = -dimensions_[2];
  aMax.z() = dimensions_[2];
}

Vector3D<Precision> UnplacedBox::GetPointOnSurface() const
{
   //copy of original UBox algorithm
   double px, py, pz, select, sumS;
   double fDx = dimensions_[0];
   double fDy = dimensions_[1];
   double fDz = dimensions_[2];
   double Sxy = fDx * fDy, Sxz = fDx * fDz, Syz = fDy * fDz;

   sumS   = Sxy + Sxz + Syz;
   select = sumS * RNG::Instance().uniform(0.,1.);

   if (select < Sxy) {
      px = -fDx + 2 * fDx * RNG::Instance().uniform(0.,1.);
      py = -fDy + 2 * fDy * RNG::Instance().uniform(0.,1.);

      if (RNG::Instance().uniform(0.,1.) > 0.5)
      {
        pz = fDz;
      }
      else
      {
        pz = -fDz;
      }
   }
   else if ((select - Sxy) < Sxz) {
      px = -fDx + 2 * fDx * RNG::Instance().uniform(0.,1.);
      pz = -fDz + 2 * fDz * RNG::Instance().uniform(0.,1.);

      if (RNG::Instance().uniform(0.,1.) > 0.5)
      {
          py = fDy;
      }
      else
      {
          py = -fDy;
      }
   }
   else {
      py = -fDy + 2 * fDy * RNG::Instance().uniform(0.,1.);
      pz = -fDz + 2 * fDz * RNG::Instance().uniform(0.,1.);

      if (RNG::Instance().uniform(0.,1.) > 0.5) {
          px = fDx;
      }
      else {
          px = -fDx;
      }
   }
   return Vector3D<Precision>(px, py, pz);
}


#ifndef VECGEOM_NVCC

template <TranslationCode trans_code, RotationCode rot_code>
VPlacedVolume* UnplacedBox::Create(
    LogicalVolume const *const logical_volume,
    Transformation3D const *const transformation,
    VPlacedVolume *const placement) {
  if (placement) {
    new(placement) SpecializedBox<trans_code, rot_code>(logical_volume,
                                                        transformation);
    return placement;
  }
  return new SpecializedBox<trans_code, rot_code>(logical_volume,
                                                  transformation);
}

VPlacedVolume* UnplacedBox::CreateSpecializedVolume(
    LogicalVolume const *const volume,
    Transformation3D const *const transformation,
    const TranslationCode trans_code, const RotationCode rot_code,
    VPlacedVolume *const placement) {
  return VolumeFactory::CreateByTransformation<UnplacedBox>(
           volume, transformation, trans_code, rot_code, placement
         );
}

#else

template <TranslationCode trans_code, RotationCode rot_code>
__device__
VPlacedVolume* UnplacedBox::Create(
    LogicalVolume const *const logical_volume,
    Transformation3D const *const transformation,
    const int id, VPlacedVolume *const placement) {
  if (placement) {
    new(placement) SpecializedBox<trans_code, rot_code>(logical_volume,
                                                        transformation, id);
    return placement;
  }
  return new SpecializedBox<trans_code, rot_code>(logical_volume,
                                                  transformation, id);
}

__device__
VPlacedVolume* UnplacedBox::CreateSpecializedVolume(
    LogicalVolume const *const volume,
    Transformation3D const *const transformation,
    const TranslationCode trans_code, const RotationCode rot_code,
    const int id, VPlacedVolume *const placement) {
  return VolumeFactory::CreateByTransformation<UnplacedBox>(
           volume, transformation, trans_code, rot_code, id, placement
         );
}

#endif

#ifdef VECGEOM_CUDA_INTERFACE

void UnplacedBox_CopyToGpu(const Precision x, const Precision y,
                           const Precision z, VUnplacedVolume *const gpu_ptr);

VUnplacedVolume* UnplacedBox::CopyToGpu(VUnplacedVolume *const gpu_ptr) const {
  UnplacedBox_CopyToGpu(this->x(), this->y(), this->z(), gpu_ptr);
  vecgeom::CudaAssertError();
  return gpu_ptr;
}

VUnplacedVolume* UnplacedBox::CopyToGpu() const {
  VUnplacedVolume *const gpu_ptr = vecgeom::AllocateOnGpu<UnplacedBox>();
  return this->CopyToGpu(gpu_ptr);
}

#endif

#ifdef VECGEOM_NVCC

class VUnplacedVolume;

__global__
void UnplacedBox_ConstructOnGpu(const Precision x, const Precision y,
                                const Precision z, VUnplacedVolume *const gpu_ptr) {
  new(gpu_ptr) vecgeom::cuda::UnplacedBox(x, y, z);
}

void UnplacedBox_CopyToGpu(const Precision x, const Precision y,
                           const Precision z, VUnplacedVolume *const gpu_ptr) {
  UnplacedBox_ConstructOnGpu<<<1, 1>>>(x, y, z, gpu_ptr);
}

#endif

} } // End global namespace
