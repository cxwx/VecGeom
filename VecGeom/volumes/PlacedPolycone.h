/*
 * PlacedPolycone.h
 *
 *  Created on: May 14, 2014
 *      Author: swenzel
 */
#ifndef VECGEOM_VOLUMES_PLACEDPOLYCONE_H_
#define VECGEOM_VOLUMES_PLACEDPOLYCONE_H_

#include "VecGeom/base/Global.h"
#ifndef VECCORE_CUDA
#include "VecGeom/base/RNG.h"
#include <cmath>
#endif
#include "VecGeom/volumes/PlacedVolume.h"
#include "VecGeom/volumes/UnplacedVolume.h"
#include "VecGeom/volumes/kernel/PolyconeImplementation.h"
#include "VecGeom/volumes/PlacedVolImplHelper.h"
#include "VecGeom/volumes/UnplacedPolycone.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedPolycone;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedPolycone);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedPolycone : public VPlacedVolume {

public:
#ifndef VECCORE_CUDA
  PlacedPolycone(char const *const label, LogicalVolume const *const logicalVolume,
                 Transformation3D const *const transformation, vecgeom::PlacedBox const *const boundingBox)
      : VPlacedVolume(label, logicalVolume, transformation, boundingBox)
  {
  }

  PlacedPolycone(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
                 vecgeom::PlacedBox const *const boundingBox)
      : PlacedPolycone("", logicalVolume, transformation, boundingBox)
  {
  }
#else
  VECCORE_ATT_DEVICE PlacedPolycone(LogicalVolume const *const logicalVolume,
                                    Transformation3D const *const transformation, PlacedBox const *const boundingBox,
                                    const int id, const int copy_no, const int child_id)
      : VPlacedVolume(logicalVolume, transformation, boundingBox, id, copy_no, child_id)
  {
  }
#endif

  VECCORE_ATT_HOST_DEVICE
  virtual ~PlacedPolycone() {}

  VECCORE_ATT_HOST_DEVICE
  UnplacedPolycone const *GetUnplacedVolume() const
  {
    return static_cast<UnplacedPolycone const *>(GetLogicalVolume()->GetUnplacedVolume());
  }

  bool IsOpen() const { return (GetUnplacedVolume()->GetDeltaPhi() < kTwoPi); }
  Precision GetStartPhi() const { return GetUnplacedVolume()->GetStartPhi(); }
  Precision GetEndPhi() const { return GetUnplacedVolume()->GetEndPhi(); }
  int GetNumRZCorner() const { return 2 * (int)(GetUnplacedVolume()->GetNz()); } // in USolids nCorners = 2*nPlanes

  PolyconeHistorical *GetOriginalParameters() const { return GetUnplacedVolume()->GetOriginalParameters(); }

  void Reset() { const_cast<UnplacedPolycone *>(GetUnplacedVolume())->Reset(); }

#ifndef VECCORE_CUDA
  virtual VPlacedVolume const *ConvertToUnspecialized() const override;
#ifdef VECGEOM_ROOT
  virtual TGeoShape const *ConvertToRoot() const override;
#endif
#ifdef VECGEOM_GEANT4
  virtual G4VSolid const *ConvertToGeant4() const override;
#endif
#endif // VECGEOM_BENCHMARK

#if !defined(VECCORE_CUDA)
  virtual Precision Capacity() override { return GetUnplacedVolume()->Capacity(); }

  virtual void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override
  {
    GetUnplacedVolume()->Extent(aMin, aMax);
  }

  VECCORE_ATT_HOST_DEVICE
  virtual bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const override
  {
    return GetUnplacedVolume()->Normal(point, normal);
  }

  virtual double SurfaceArea() const override { return GetUnplacedVolume()->SurfaceArea(); }
#endif

}; // end class

template <typename UnplacedPolycone_t>
class SPlacedPolycone : public PlacedVolumeImplHelper<UnplacedPolycone_t, PlacedPolycone> {
  using Base = PlacedVolumeImplHelper<UnplacedPolycone_t, PlacedPolycone>;

public:
  typedef UnplacedPolycone UnplacedShape_t;
  using Base::Base;
};

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_PLACEDPOLYCONE_H_
