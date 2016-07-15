#include "volumes/PlacedSExtru.h"
#include "base/SOA3D.h"
#include "volumes/SpecializedSExtru.h"

#ifdef VECGEOM_ROOT
#include "TGeoXtru.h"
#endif

#ifdef VECGEOM_GEANT4
#include "G4ExtrudedSolid.hh"
#include "G4TwoVector.hh"
#include <vector>
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECGEOM_CUDA_HEADER_BOTH
void PlacedSExtru::PrintType() const
{
  printf("PlacedSExtru");
}

void PlacedSExtru::PrintType(std::ostream &s) const
{
  s << "PlacedSExtru";
}

Vector3D<Precision> PlacedSExtru::GetPointOnSurface() const
{
  // use generic implementation
  return VPlacedVolume::GetPointOnSurface();
}

// Comparison specific
#ifndef VECGEOM_NVCC
VPlacedVolume const *PlacedSExtru::ConvertToUnspecialized() const
{
  return this;
}
#ifdef VECGEOM_ROOT
TGeoShape const *PlacedSExtru::ConvertToRoot() const
{
  TGeoXtru *s = new TGeoXtru(2);

  // get vertices and make array to construct the ROOT volume
  auto &vertices = GetUnplacedStruct()->GetPolygon().GetVertices();
  const auto N   = vertices.size();
  s->DefinePolygon(N, vertices.x(), vertices.y());
  s->DefineSection(0, GetUnplacedStruct()->GetLowerZ(), 0., 0., 1.);
  s->DefineSection(1, GetUnplacedStruct()->GetUpperZ(), 0., 0., 1.);
  return s;
};
#endif
#if defined(VECGEOM_USOLIDS) && !defined(VECGEOM_REPLACE_USOLIDS)
::VUSolid const *PlacedSExtru::ConvertToUSolids() const
{
  return nullptr;
}
#endif
#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedSExtru::ConvertToGeant4() const
{
  std::vector<G4TwoVector> G4vertices;
  using ZSection = G4ExtrudedSolid::ZSection;
  std::vector<ZSection> zsections;
  auto &vertices = GetUnplacedStruct()->GetPolygon().GetVertices();
  for (size_t i = 0; i < vertices.size(); ++i) {
    G4vertices.push_back(G4TwoVector(vertices.x(i), vertices.y(i)));
  }
  zsections.push_back(ZSection(GetUnplacedStruct()->GetLowerZ(), 0., 1.));
  zsections.push_back(ZSection(GetUnplacedStruct()->GetUpperZ(), 0., 1.));
  G4String s("g4extru");
  return new G4ExtrudedSolid(s, G4vertices, zsections);
};
#endif
#endif // VECGEOM_NVCC

} // End impl namespace

#ifdef VECGEOM_NVCC

VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC(SpecializedSExtru)

#endif // VECGEOM_NVCC

} // End namespace vecgeom