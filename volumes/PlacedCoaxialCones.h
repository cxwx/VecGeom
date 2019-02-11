/// @file PlacedCoaxialCones.h
/// @author Raman Sehgal (raman.sehgal@cern.ch), Evgueni Tcherniaev (evgueni.tcherniaev@cern.ch)

#ifndef VECGEOM_VOLUMES_PLACEDCOAXIALCONES_H_
#define VECGEOM_VOLUMES_PLACEDCOAXIALCONES_H_

#include "base/Global.h"

#include "volumes/PlacedVolume.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/kernel/CoaxialConesImplementation.h"
#include "volumes/PlacedVolImplHelper.h"
#include "volumes/UnplacedCoaxialCones.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedCoaxialCones;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedCoaxialCones);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedCoaxialCones : public PlacedVolumeImplHelper<UnplacedCoaxialCones, VPlacedVolume> {
  using Base = PlacedVolumeImplHelper<UnplacedCoaxialCones, VPlacedVolume>;

public:
#ifndef VECCORE_CUDA
  // constructor inheritance;
  using Base::Base;
  PlacedCoaxialCones(char const *const label, LogicalVolume const *const logicalVolume,
                       Transformation3D const *const transformation,
                       vecgeom::PlacedBox const *const boundingBox)
      : Base(label, logicalVolume, transformation, boundingBox)
  {
  }

  PlacedCoaxialCones(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
                       vecgeom::PlacedBox const *const boundingBox)
      : PlacedCoaxialCones("", logicalVolume, transformation, boundingBox)
  {
  }
#else
  VECCORE_ATT_DEVICE PlacedCoaxialCones(LogicalVolume const *const logicalVolume,
                                          Transformation3D const *const transformation,
                                          PlacedBox const *const boundingBox, const int id)
      : Base(logicalVolume, transformation, boundingBox, id)
  {
  }
#endif
  VECCORE_ATT_HOST_DEVICE
  virtual ~PlacedCoaxialCones() {}

  /*
   * Put the Required Getters and Setters here
   */


  VECCORE_ATT_HOST_DEVICE
  virtual void PrintType() const override;
  virtual void PrintType(std::ostream &os) const override;

// Comparison specific
#ifndef VECCORE_CUDA
  virtual VPlacedVolume const *ConvertToUnspecialized() const override;
#ifdef VECGEOM_ROOT
  virtual TGeoShape const *ConvertToRoot() const override;
#endif
#ifdef VECGEOM_GEANT4
  virtual G4VSolid const *ConvertToGeant4() const override;
#endif
#endif // VECCORE_CUDA
};

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
