/*
 * PlacedCone.cpp
 *
 *  Created on: Jun 13, 2014
 *      Author: swenzel
 */

#include "volumes/PlacedCone.h"
#include "volumes/Cone.h"
#include "volumes/SpecializedCone.h"

#if defined(VECGEOM_ROOT)
#include "TGeoCone.h"
#endif

#if defined(VECGEOM_USOLIDS)
#include "UCons.hh"
#endif

#if defined(VECGEOM_GEANT4)
#include "G4Cons.hh"
#endif


namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECGEOM_NVCC
VPlacedVolume const* PlacedCone::ConvertToUnspecialized() const {
    return new SimpleCone(GetLabel().c_str(), logical_volume(), transformation());
}

#ifdef VECGEOM_ROOT
TGeoShape const* PlacedCone::ConvertToRoot() const {
    if( GetDPhi() == 2.*M_PI )
    {
       return new TGeoCone("RootCone",GetDz(),GetRmin1(),GetRmax1(), GetRmin2(), GetRmax2());
    }
    else
    {
       return new TGeoConeSeg("RootCone", GetDz(),GetRmin1(),GetRmax1(), GetRmin2(), GetRmax2(),
               GetSPhi()*(180/M_PI), GetSPhi()+180*GetDPhi()/(M_PI) );
    }
}
#endif

#ifdef VECGEOM_USOLIDS
::VUSolid const* PlacedCone::ConvertToUSolids() const {
  return new UCons("USolidCone", GetRmin1(), GetRmax1(), GetRmin2(), GetRmax2(), GetDz(), GetSPhi(), GetDPhi());
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const * PlacedCone::ConvertToGeant4() const {
  return new G4Cons("Geant4Cone", GetRmin1(), GetRmax1(), GetRmin2(), GetRmax2(), GetDz(), GetSPhi(), GetDPhi());
}
#endif

#endif

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VPlacedVolume> PlacedCone::CopyToGpu(
   DevicePtr<cuda::LogicalVolume> const logical_volume,
   DevicePtr<cuda::Transformation3D> const transform,
   DevicePtr<cuda::VPlacedVolume> const gpu_ptr) const
{
   return CopyToGpuImpl<PlacedCone>(logical_volume, transform, gpu_ptr);
}

DevicePtr<cuda::VPlacedVolume> PlacedCone::CopyToGpu(
      DevicePtr<cuda::LogicalVolume> const logical_volume,
      DevicePtr<cuda::Transformation3D> const transform) const
{
   return CopyToGpuImpl<PlacedCone>(logical_volume, transform);
}

#endif // VECGEOM_CUDA_INTERFACE

} // End impl namespace

#ifdef VECGEOM_NVCC

namespace cxx {

template size_t DevicePtr<cuda::PlacedCone>::SizeOf();
#ifdef HAS_PLACED_IMPL
template void DevicePtr<cuda::PlacedCone>::Construct(
   DevicePtr<cuda::LogicalVolume> const logical_volume,
   DevicePtr<cuda::Transformation3D> const transform,
   const int id) const;
#endif

} // End cxx namespace

#endif // VECGEOM_NVCC

} // End namespace vecgeom
