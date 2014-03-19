/**
 * @file specialized_box.h
 * @author Johannes de Fine Licht (johannes.definelicht@cern.ch)
 */

#ifndef VECGEOM_VOLUMES_SPECIALIZEDBOX_H_
#define VECGEOM_VOLUMES_SPECIALIZEDBOX_H_

#ifdef VECGEOM_CUDA
#include <stdio.h>
#endif

#include "base/global.h"
#include "backend.h"
#include "implementation.h"

#include "base/transformation_matrix.h"
#include "volumes/placed_box.h"

namespace vecgeom {

template <TranslationCode trans_code, RotationCode rot_code>
class SpecializedBox : public PlacedBox {

public:

  VECGEOM_CUDA_HEADER_BOTH
  SpecializedBox(LogicalVolume const *const logical_volume,
                 TransformationMatrix const *const matrix)
      : PlacedBox(logical_volume, matrix) {}

  VECGEOM_CUDA_HEADER_BOTH
  virtual bool Inside(Vector3D<Precision> const &point) const;

  VECGEOM_CUDA_HEADER_BOTH
  virtual bool Inside(Vector3D<Precision> const &point,
                      Vector3D<Precision> &local) const;

  virtual void Inside(SOA3D<Precision> const &points,
                      bool *const output) const;

  virtual void Inside(AOS3D<Precision> const &points,
                      bool *const output) const;

  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision DistanceToIn(Vector3D<Precision> const &position,
                                 Vector3D<Precision> const &direction,
                                 const Precision step_max) const;

  virtual void DistanceToIn(SOA3D<Precision> const &position,
                            SOA3D<Precision> const &direction,
                            Precision const *const step_max,
                            Precision *const output) const;

  virtual void DistanceToIn(AOS3D<Precision> const &position,
                            AOS3D<Precision> const &direction,
                            Precision const *const step_max,
                            Precision *const output) const;

  virtual int memory_size() const { return sizeof(*this); }

  #ifdef VECGEOM_CUDA
  virtual VPlacedVolume* CopyToGpu(LogicalVolume const *const logical_volume,
                                   TransformationMatrix const *const matrix,
                                   VPlacedVolume *const gpu_ptr) const;
  virtual VPlacedVolume* CopyToGpu(
      LogicalVolume const *const logical_volume,
      TransformationMatrix const *const matrix) const;
  #endif

};

template <TranslationCode trans_code, RotationCode rot_code>
VECGEOM_CUDA_HEADER_BOTH
bool SpecializedBox<trans_code, rot_code>::Inside(
    Vector3D<Precision> const &point) const {
  return InsideDispatch<trans_code, rot_code, kScalar>(point);
}

template <TranslationCode trans_code, RotationCode rot_code>
VECGEOM_CUDA_HEADER_BOTH
bool SpecializedBox<trans_code, rot_code>::Inside(
    Vector3D<Precision> const &point,
    Vector3D<Precision> &local) const {
  bool output;
  BoxInside<trans_code, rot_code, kScalar>(
    unplaced_box()->dimensions(),
    *this->matrix_,
    point,
    local,
    &output
  );
  return output;
}

template <TranslationCode trans_code, RotationCode rot_code>
void SpecializedBox<trans_code, rot_code>::Inside(
    SOA3D<Precision> const &points,
    bool *const output) const {
  Inside_Looper<trans_code, rot_code>(*this, points, output);
}

template <TranslationCode trans_code, RotationCode rot_code>
void SpecializedBox<trans_code, rot_code>::Inside(
    AOS3D<Precision> const &points,
    bool *const output) const {
  Inside_Looper<trans_code, rot_code>(*this, points, output);
}

template <TranslationCode trans_code, RotationCode rot_code>
VECGEOM_CUDA_HEADER_BOTH
Precision SpecializedBox<trans_code, rot_code>::DistanceToIn(
    Vector3D<Precision> const &position,
    Vector3D<Precision> const &direction,
    const Precision step_max) const {

  return DistanceToInDispatch<trans_code, rot_code, kScalar>(
           position, direction, step_max
         );
                                                  
}

template <TranslationCode trans_code, RotationCode rot_code>
void SpecializedBox<trans_code, rot_code>::DistanceToIn(
    SOA3D<Precision> const &positions,
    SOA3D<Precision> const &directions,
    Precision const *const step_max,
    Precision *const output) const {
  DistanceToIn_Looper<trans_code, rot_code>(*this, positions, directions,
                                            step_max, output);
}

template <TranslationCode trans_code, RotationCode rot_code>
void SpecializedBox<trans_code, rot_code>::DistanceToIn(
    AOS3D<Precision> const &positions,
    AOS3D<Precision> const &directions,
    Precision const *const step_max,
    Precision *const output) const {
  DistanceToIn_Looper<trans_code, rot_code>(*this, positions, directions,
                                            step_max, output);
}

#ifdef VECGEOM_CUDA

namespace {

template <TranslationCode trans_code, RotationCode rot_code>
__global__
void ConstructOnGpu(LogicalVolume const *const logical_volume,
                    TransformationMatrix const *const matrix,
                    VPlacedVolume *const gpu_ptr) {
  new(gpu_ptr) SpecializedBox<trans_code, rot_code>(logical_volume, matrix);
}

} // End anonymous namespace

template <TranslationCode trans_code, RotationCode rot_code>
VPlacedVolume* SpecializedBox<trans_code, rot_code>::CopyToGpu(
    LogicalVolume const *const logical_volume,
    TransformationMatrix const *const matrix,
    VPlacedVolume *const gpu_ptr) const {

  ConstructOnGpu<trans_code, rot_code><<<1, 1>>>(
    logical_volume, matrix, gpu_ptr
  );
  CudaAssertError();
  return gpu_ptr;

}

template <TranslationCode trans_code, RotationCode rot_code>
VPlacedVolume* SpecializedBox<trans_code, rot_code>::CopyToGpu(
    LogicalVolume const *const logical_volume,
    TransformationMatrix const *const matrix) const {

  VPlacedVolume *const gpu_ptr =
      AllocateOnGpu<SpecializedBox<trans_code, rot_code> >();
  return CopyToGpu(logical_volume, matrix, gpu_ptr);  

}

#endif // VECGEOM_CUDA

} // End namespace vecgeom

#endif // VECGEOM_VOLUMES_SPECIALIZEDBOX_H_
