#include <VecCore/VecCore>

#include "volumes/TessellatedStruct.h"
#include "test/benchmark/ArgParser.h"
#include "volumes/utilities/VolumeUtilities.h"

#ifdef VECGEOM_ROOT
#include "utilities/Visualizer.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "volumes/Box.h"
#endif

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

/* Simple test for the scalability of creation of the tessellated structure.
   An orb is split into ngrid theta and phi regions; each cell is represented
   as a quadrilateral. The solid will contain 2*(ngrid-1)*ngrid triangle facets */

using namespace vecgeom;

double r    = 10.;
double *sth = nullptr;
double *cth = nullptr;
double *sph = nullptr;
double *cph = nullptr;

void RandomDirection(Vector3D<double> &direction)
{
  double phi    = RNG::Instance().uniform(0., 2. * kPi);
  double theta  = std::acos(1. - 2. * RNG::Instance().uniform(0, 1));
  direction.x() = std::sin(theta) * std::cos(phi);
  direction.y() = std::sin(theta) * std::sin(phi);
  direction.z() = std::cos(theta);
}

VECGEOM_FORCE_INLINE
Vector3D<double> Vtx(int ith, int iph)
{
  return Vector3D<double>(r * sth[ith] * cph[iph], r * sth[ith] * sph[iph], r * cth[ith]);
}

#ifdef VECGEOM_ROOT
void AddFacetToVisualizer(TriangleFacet<double> const &facet, Visualizer &visualizer)
{
  TPolyLine3D pl(3);
  pl.SetLineColor(kBlue);
  for (int i = 0; i < 3; i++)
    pl.SetNextPoint(facet.fVertices[i].x(), facet.fVertices[i].y(), facet.fVertices[i].z());
  visualizer.AddLine(pl);
}
#endif

int main(int argc, char *argv[])
{
  using namespace vecgeom;
  //  using Real_v = typename VectorBackend::Real_v;

  OPTION_INT(ngrid, 100);

  // Sin/Cos tables
  double dth = kPi / ngrid;
  double dph = kTwoPi / ngrid;
  sth        = new double[ngrid + 1];
  cth        = new double[ngrid + 1];
  sph        = new double[ngrid + 1];
  cph        = new double[ngrid + 1];

  for (int i = 0; i <= ngrid; ++i) {
    sth[i] = vecCore::math::Sin(i * dth);
    cth[i] = vecCore::math::Cos(i * dth);
    sph[i] = vecCore::math::Sin(i * dph);
    cph[i] = vecCore::math::Cos(i * dph);
  }

  TessellatedStruct<double> tsl;
  for (int ith = 0; ith < ngrid; ++ith) {
    for (int iph = 0; iph < ngrid; ++iph) {
      // First/last rows - > triangles
      if (ith == 0) {
        tsl.AddTriangularFacet(Vector3D<double>(0, 0, r), Vtx(ith, iph), Vtx(ith, iph + 1));
      } else if (ith == ngrid - 1) {
        tsl.AddTriangularFacet(Vtx(ith, iph), Vector3D<double>(0, 0, -r), Vtx(ith, iph + 1));
      } else {
        tsl.AddQuadrilateralFacet(Vtx(ith, iph), Vtx(ith + 1, iph), Vtx(ith + 1, iph + 1), Vtx(ith, iph + 1));
      }
    }
  }

  // Close the solid
  tsl.Close();

// Visualize the facets
#ifdef VECGEOM_ROOT
  Visualizer visualizer;
  // Visualize bounding box
  Vector3D<double> deltas = 0.5 * (tsl.fMaxExtent - tsl.fMinExtent);
  Vector3D<double> origin = 0.5 * (tsl.fMaxExtent + tsl.fMinExtent);
  SimpleBox box("bbox", deltas.x(), deltas.y(), deltas.z());
  visualizer.AddVolume(box, Transformation3D(origin.x(), origin.y(), origin.z()));

  // Visualize facets
  for (auto facet : tsl.fFacets) {
    AddFacetToVisualizer(facet, visualizer);
  }
  visualizer.Show();
#endif

  return 0;
}
