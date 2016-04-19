/// \file UnplacedParallelepiped.cpp
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#include "volumes/UnplacedParallelepiped.h"

#include <stdio.h>
#include "management/VolumeFactory.h"
#include "volumes/SpecializedParallelepiped.h"
#include "volumes/utilities/GenerationUtilities.h"

#ifndef VECGEOM_NVCC
#include "base/RNG.h"
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECGEOM_CUDA_HEADER_BOTH
UnplacedParallelepiped::UnplacedParallelepiped(Vector3D<Precision> const &dimensions, const Precision alpha,
                                               const Precision theta, const Precision phi)
    : fDimensions(dimensions), fAlpha(0), fTheta(0), fPhi(0), fTanAlpha(0), fTanThetaSinPhi(0), fTanThetaCosPhi(0),
      fCtx(0), fCty(0) {
  SetAlpha(alpha);
  SetThetaAndPhi(theta, phi);
}

//______________________________________________________________________________
VECGEOM_CUDA_HEADER_BOTH
UnplacedParallelepiped::UnplacedParallelepiped(const Precision x, const Precision y, const Precision z,
                                               const Precision alpha, const Precision theta, const Precision phi)
    : fDimensions(x, y, z), fAlpha(0), fTheta(0), fPhi(0), fTanAlpha(0), fTanThetaSinPhi(0), fTanThetaCosPhi(0),
      fCtx(0), fCty(0) {
  SetAlpha(alpha);
  SetThetaAndPhi(theta, phi);
}

//______________________________________________________________________________
VECGEOM_CUDA_HEADER_BOTH
void UnplacedParallelepiped::SetAlpha(const Precision alpha) {
  fAlpha = alpha;
  fTanAlpha = tan(kDegToRad * alpha);
  ComputeNormals();
}

//______________________________________________________________________________
VECGEOM_CUDA_HEADER_BOTH
void UnplacedParallelepiped::SetTheta(const Precision theta) { SetThetaAndPhi(theta, fPhi); }

//______________________________________________________________________________
VECGEOM_CUDA_HEADER_BOTH
void UnplacedParallelepiped::SetPhi(const Precision phi) { SetThetaAndPhi(fTheta, phi); }

//______________________________________________________________________________
VECGEOM_CUDA_HEADER_BOTH
void UnplacedParallelepiped::SetThetaAndPhi(const Precision theta, const Precision phi) {
  fTheta = theta;
  fPhi = phi;
  fTanThetaCosPhi = tan(kDegToRad * fTheta) * cos(kDegToRad * fPhi);
  fTanThetaSinPhi = tan(kDegToRad * fTheta) * sin(kDegToRad * fPhi);
  ComputeNormals();
}

//______________________________________________________________________________
VECGEOM_CUDA_HEADER_BOTH
void UnplacedParallelepiped::ComputeNormals() {
  Vector3D<Precision> v(sin(kDegToRad * fTheta) * cos(kDegToRad * fPhi),
                        sin(kDegToRad * fTheta) * sin(kDegToRad * fPhi), cos(kDegToRad * fTheta));
  Vector3D<Precision> vx(1., 0., 0.);
  Vector3D<Precision> vy(-sin(kDegToRad * fAlpha), -cos(kDegToRad * fAlpha), 0.);
  fNormals[0] = v.Cross(vy);
  fNormals[0].Normalize();
  fNormals[1] = v.Cross(vx);
  fNormals[1].Normalize();
  fNormals[2].Set(0., 0., 1.);
  fCtx = sqrt(1. + fTanAlpha * fTanAlpha + fTanThetaCosPhi * fTanThetaCosPhi);
  fCty = sqrt(1. + fTanThetaSinPhi * fTanThetaSinPhi);
  std::cout << "norm[0] = " << fNormals[0] << " norm[1] = " << fNormals[1] << "\n";
}

//______________________________________________________________________________
void UnplacedParallelepiped::Print() const {
  printf("UnplacedParallelepiped {%.2f, %.2f, %.2f, %.2f, %.2f, %.2f}", GetX(), GetY(), GetZ(), GetTanAlpha(),
         GetTanThetaCosPhi(), GetTanThetaSinPhi());
}

//______________________________________________________________________________
void UnplacedParallelepiped::Print(std::ostream &os) const {
  os << "UnplacedParallelepiped {" << GetX() << ", " << GetY() << ", " << GetZ() << ", " << GetTanAlpha() << ", "
     << GetTanThetaCosPhi() << ", " << GetTanThetaSinPhi();
}

//______________________________________________________________________________
VECGEOM_CUDA_HEADER_BOTH
void UnplacedParallelepiped::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const {
  // Returns the full 3D cartesian extent of the solid.
  Precision dx = fDimensions[0] + fDimensions[1] * Abs(fTanAlpha) + fDimensions[2] * Abs(fTanThetaCosPhi);
  Precision dy = fDimensions[1] + fDimensions[2] * Abs(fTanThetaSinPhi);
  Precision dz = fDimensions[2];
  aMin.Set(-dx, -dy, -dz);
  aMax.Set(dx, dy, dz);
}

//______________________________________________________________________________
Vector3D<Precision> UnplacedParallelepiped::GetPointOnSurface() const {
  // Generate randomly a point on one of the surfaces
  // Select randomly a surface
  Vector3D<Precision> point;
  int isurf = int(RNG::Instance().uniform(0., 6.));
  Precision dx = 0, dy = 0, dz = 0;
  switch (isurf) {
  case 0:
  case 1:
    // top/bottom
    dx = RNG::Instance().uniform(-fDimensions[0], fDimensions[0]);
    dy = RNG::Instance().uniform(-fDimensions[1], fDimensions[1]);
    dz = (2 * isurf - 1) * fDimensions[2];
    break;
  case 2:
  case 3:
    // front/behind
    dx = RNG::Instance().uniform(-fDimensions[0], fDimensions[0]);
    dy = (2 * (isurf - 2) - 1) * fDimensions[1];
    dz = RNG::Instance().uniform(-fDimensions[2], fDimensions[2]);
    break;
  case 4:
  case 5:
    // left/right
    dx = (2 * (isurf - 4) - 1) * fDimensions[0];
    dy = RNG::Instance().uniform(-fDimensions[1], fDimensions[1]);
    dz = RNG::Instance().uniform(-fDimensions[2], fDimensions[2]);
  }
  point[0] = dx + dy * fTanAlpha + dz * fTanThetaCosPhi;
  point[1] = dy + dz * fTanThetaSinPhi;
  point[2] = dz;
  return (point);
}

//______________________________________________________________________________
VECGEOM_CUDA_HEADER_BOTH
bool UnplacedParallelepiped::IsConvex() const {
  // A parallelepiped is convex shape
  return true;
}

//______________________________________________________________________________
template <TranslationCode transCodeT, RotationCode rotCodeT>
VECGEOM_CUDA_HEADER_DEVICE VPlacedVolume *UnplacedParallelepiped::Create(LogicalVolume const *const logical_volume,
                                                                         Transformation3D const *const transformation,
#ifdef VECGEOM_NVCC
                                                                         const int id,
#endif
                                                                         VPlacedVolume *const placement) {

  return CreateSpecializedWithPlacement<SpecializedParallelepiped<transCodeT, rotCodeT>>(
#ifdef VECGEOM_NVCC
      logical_volume, transformation, id, placement); // TODO: add bounding box?
#else
      logical_volume, transformation, placement);
#endif
}

//______________________________________________________________________________
VECGEOM_CUDA_HEADER_DEVICE
VPlacedVolume *UnplacedParallelepiped::SpecializedVolume(LogicalVolume const *const volume,
                                                         Transformation3D const *const transformation,
                                                         const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
                                                         const int id,
#endif
                                                         VPlacedVolume *const placement) const {
  return VolumeFactory::CreateByTransformation<UnplacedParallelepiped>(volume, transformation, trans_code, rot_code,
#ifdef VECGEOM_NVCC
                                                                       id,
#endif
                                                                       placement);
}

#ifdef VECGEOM_CUDA_INTERFACE

//______________________________________________________________________________
DevicePtr<cuda::VUnplacedVolume>
UnplacedParallelepiped::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const {
  return CopyToGpuImpl<UnplacedParallelepiped>(in_gpu_ptr, GetX(), GetY(), GetZ(), fAlpha, fTheta, fPhi);
}

//______________________________________________________________________________
DevicePtr<cuda::VUnplacedVolume> UnplacedParallelepiped::CopyToGpu() const {
  return CopyToGpuImpl<UnplacedParallelepiped>();
}

#endif // VECGEOM_CUDA_INTERFACE

} // End impl namespace

#ifdef VECGEOM_NVCC

namespace cxx {

template size_t DevicePtr<cuda::UnplacedParallelepiped>::SizeOf();
template void DevicePtr<cuda::UnplacedParallelepiped>::Construct(const Precision x, const Precision y,
                                                                 const Precision z, const Precision alpha,
                                                                 const Precision theta, const Precision phi) const;

} // End cxx namespace

#endif

} // End global namespace
