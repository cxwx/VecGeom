/**
 * @file unplaced_box.cpp
 * @author Johannes de Fine Licht (johannes.definelicht@cern.ch)
 */

#include <stdio.h>
#include "backend/backend.h"
#include "management/volume_factory.h"
#include "volumes/specialized_box.h"
#include "volumes/unplaced_box.h"

namespace VECGEOM_NAMESPACE {

template <TranslationCode trans_code, RotationCode rot_code>
VECGEOM_CUDA_HEADER_BOTH
VPlacedVolume* UnplacedBox::Create(
    LogicalVolume const *const logical_volume,
    TransformationMatrix const *const matrix,
    VPlacedVolume *const placement) {
  if (placement) {
    new(placement) SpecializedBox<trans_code, rot_code>(logical_volume, matrix);
    return placement;
  }
  return new SpecializedBox<trans_code, rot_code>(logical_volume, matrix);
}

VECGEOM_CUDA_HEADER_BOTH
VPlacedVolume* UnplacedBox::CreateSpecializedVolume(
    LogicalVolume const *const volume,
    TransformationMatrix const *const matrix,
    const TranslationCode trans_code, const RotationCode rot_code,
    VPlacedVolume *const placement) {
  return VolumeFactory::CreateByTransformation<UnplacedBox>(
           volume, matrix, trans_code, rot_code, placement
         );
}

VECGEOM_CUDA_HEADER_BOTH
void UnplacedBox::Print() const {
  printf("Box {%f, %f, %f}", dimensions_[0], dimensions_[1], dimensions_[2]);
}

} // End global namespace

namespace vecgeom {

#ifdef VECGEOM_CUDA_INTERFACE

void GpuInterface(const Precision x, const Precision y, const Precision z,
                  VUnplacedVolume *const gpu_ptr);

VUnplacedVolume* UnplacedBox::CopyToGpu(VUnplacedVolume *const gpu_ptr) const {
  GpuInterface(this->x(), this->y(), this->z(), gpu_ptr);
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
void ConstructOnGpu(const Precision x, const Precision y, const Precision z,
                    VUnplacedVolume *const gpu_ptr) {
  new(gpu_ptr) vecgeom_cuda::UnplacedBox(x, y, z);
}

void GpuInterface(const Precision x, const Precision y, const Precision z,
                  VUnplacedVolume *const gpu_ptr) {
  ConstructOnGpu<<<1, 1>>>(x, y, z, gpu_ptr);
}

#endif

} // End global namespace