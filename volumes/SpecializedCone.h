/*
 * SpecializedCone.h
 *
 *  Created on: May 14, 2014
 *      Author: swenzel
 */

#ifndef VECGEOM_VOLUMES_SPECIALIZEDCONE_H_
#define VECGEOM_VOLUMES_SPECIALIZEDCONE_H_

#include "base/Global.h"
#include "volumes/kernel/ConeImplementation.h"
#include "volumes/PlacedCone.h"
#include "volumes/SpecializedPlacedVolImplHelper.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <TranslationCode transCodeT, RotationCode rotCodeT, typename coneTypeT>
using SpecializedCone = SIMDSpecializedVolImplHelper<ConeImplementation<coneTypeT>, transCodeT, rotCodeT>;

using SimpleCone = SpecializedCone<translation::kGeneric, rotation::kGeneric, ConeTypes::UniversalCone>;

template <typename Type>
template <TranslationCode transCodeT, RotationCode rotCodeT>
VECCORE_ATT_DEVICE
VPlacedVolume *SUnplacedCone<Type>::Create(LogicalVolume const *const logical_volume,
                                           Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                                           const int id,
#endif
                                           VPlacedVolume *const placement)
{
  (void)placement;
  return new SpecializedCone<transCodeT, rotCodeT, Type>(logical_volume, transformation
#ifdef VECCORE_CUDA
                                                         ,
                                                         id
#endif
                                                         );
}

} // End VECGEOM_IMPL_NAMESPACE
} // End global namespace

#endif // VECGEOM_VOLUMES_SPECIALIZEDCONE_H_
