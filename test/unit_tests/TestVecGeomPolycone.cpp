/*
 * TestVecGeomPolycone.cpp
 *
 *  Created on: Dec 8, 2014
 *      Author: swenzel
 */

// a file to test compilation of features during development

//.. ensure asserts are compiled in
#undef NDEBUG

#include "volumes/Polycone.h"
#include "volumes/Tube.h"
#include "volumes/Cone.h"
#include "volumes/LogicalVolume.h"
#include "volumes/PlacedVolume.h"
#include "base/Vector3D.h"
#include "ApproxEqual.h"

#ifdef VECGEOM_USOLIDS
#include "UPolycone.hh"
#include "UGenericPolycone.hh"
#include "UVector3.hh"
#endif

using Real_v = vecgeom::VectorBackend::Real_v;

#include <iostream>

using namespace vecgeom;

bool testingvecgeom = false;

typedef Vector3D<Precision> Vec3D_t;

template <class Polycone_t, class Vec_t = vecgeom::Vector3D<vecgeom::Precision>>

bool TestPolycone()
{
  return 0;
}

int main()
{
  int Nz = 4;
  // a tube and two cones
  double rmin[] = {0.1, 0.0, 0.0, 0.4};
  double rmax[] = {1., 2., 2., 1.5};
  double z[]    = {-1, -0.5, 0.5, 2};

  UnplacedPolycone poly1(0.,      /* initial phi starting angle */
                         kTwoPi,  /* total phi angle */
                         Nz,      /* number corners in r,z space */
                         z, rmin, /* r coordinate of these corners */
                         rmax);

  // poly1.Print();

  // lets make external separate tubes and cones representing the sections
  UnplacedCone section0(rmin[0], rmax[0], rmin[1], rmax[1], (z[1] - z[0]) / 2., 0, kTwoPi);
  UnplacedCone section1(rmin[1], rmax[1], rmin[2], rmax[2], (z[2] - z[1]) / 2., 0, kTwoPi);
  UnplacedCone section2(rmin[2], rmax[2], rmin[3], rmax[3], (z[3] - z[2]) / 2., 0, kTwoPi);

  assert(poly1.GetNz() == 4);
  assert(poly1.GetNSections() == 3);
  assert(poly1.GetSectionIndex(-0.8) == 0);
  assert(poly1.GetSectionIndex(0.51) == 2);
  assert(poly1.GetSectionIndex(0.) == 1);
  assert(poly1.GetSectionIndex(-2.) == -1);
  assert(poly1.GetSectionIndex(3.) == -2);
  assert(poly1.GetStartPhi() == 0.);
  assert((std::fabs(poly1.GetDeltaPhi() - kTwoPi)) < 1e-10);

  assert(poly1.GetStruct().fZs[0] == z[0]);
  assert(poly1.GetStruct().fZs[poly1.GetNSections()] == z[Nz - 1]);
  assert(poly1.Capacity() > 0);
  assert(std::fabs(poly1.Capacity() - (section0.Capacity() + section1.Capacity() + section2.Capacity())) < 1e-6);

  // create a place version
  VPlacedVolume const *placedpoly1 = (new LogicalVolume("poly1", &poly1))->Place(new Transformation3D());

  // test contains/inside
  assert(placedpoly1->Contains(Vec3D_t(0., 0., 0.)) == true);
  assert(placedpoly1->Contains(Vec3D_t(0., 0., -2.)) == false);
  assert(placedpoly1->Contains(Vec3D_t(0., 0., -0.8)) == false);
  assert(placedpoly1->Contains(Vec3D_t(0., 0., -1.8)) == false);
  assert(placedpoly1->Contains(Vec3D_t(0., 0., 10)) == false);
  assert(placedpoly1->Contains(Vec3D_t(0., 0., 1.8)) == false);

  // test DistanceToIn
  assert(ApproxEqual(placedpoly1->DistanceToIn(Vec3D_t(0., 0., -3.), Vec3D_t(0., 0., 1.)), 2.5));
  assert(placedpoly1->DistanceToIn(Vec3D_t(0., 0., -2.), Vec3D_t(0., 0., -1.)) == vecgeom::kInfLength);
  assert(ApproxEqual(placedpoly1->DistanceToIn(Vec3D_t(0., 0., 3), Vec3D_t(0., 0., -1.)), 2.5));
  assert(placedpoly1->DistanceToIn(Vec3D_t(0., 0., 3), Vec3D_t(0., 0., 1.)) == vecgeom::kInfLength);
  assert(ApproxEqual(placedpoly1->DistanceToIn(Vec3D_t(3., 0., 0), Vec3D_t(-1., 0., 0.)), 1.));
  assert(ApproxEqual(placedpoly1->DistanceToIn(Vec3D_t(0., 0., 1.999999999), Vec3D_t(1., 0., 0.)), 0.4));

  {
    // test vector interface
    using Vec3D_v = vecgeom::Vector3D<Real_v>;
    Vec3D_v p(Real_v(0.), Real_v(0.), Real_v(-3)), d(Real_v(0.), Real_v(0.), Real_v(1.));
    const auto dist = placedpoly1->DistanceToInVec(p, d);
    for (size_t lane = 0; lane < vecCore::VectorSize<Real_v>(); ++lane) {
      assert(ApproxEqual(vecCore::LaneAt(dist, lane), 2.5));
    }
  }

  // test SafetyToIn
  // assert( placedpoly1-> SafetyToIn( Vec3D_t(0.,0.,-3.)) == 2. );
  // assert( placedpoly1-> SafetyToIn( Vec3D_t(0.5,0.,-1.)) == 0. );
  // assert( placedpoly1-> SafetyToIn( Vec3D_t(0.,0.,3) ) == 1 );
  // assert( placedpoly1-> SafetyToIn( Vec3D_t(2.,0.,0.1) ) == 0 );
  // assert( placedpoly1-> SafetyToIn( Vec3D_t(3.,0.,0) ) == 1. );
  // std::cout << "test SIN " << placedpoly1->SafetyToIn(Vec3D_t(8., 0., 0.)) << std::endl;
  // test SafetyToOut
  // std::cout << "test SOUT " << placedpoly1->SafetyToOut(Vec3D_t(0., 0., 0.)) << std::endl;
  // assert( placedpoly1-> SafetyToOut( Vec3D_t(0.,0.,0.)) == 0.5 );
  //  assert( placedpoly1-> SafetyToOut( Vec3D_t(0.,0.,0.5)) == 0. );
  //  assert( placedpoly1-> SafetyToOut( Vec3D_t(1.9,0.,0) ) == 0.1 );
  // assert( placedpoly1-> SafetyToOut( Vec3D_t(0.2,0.,-1) ) == 0. );
  //  assert( placedpoly1-> SafetyToOut( Vec3D_t(1.4,0.,2) ) == 0. );

  // test DistanceToOut
  // [see VECGEOM-414] assert(ApproxEqual(placedpoly1->DistanceToOut(Vec3D_t(0., 0., 0.), Vec3D_t(0., 0., 1.)), 0.5));
  assert(ApproxEqual(placedpoly1->DistanceToOut(Vec3D_t(0., 0., 0.), Vec3D_t(0., 0., -1.)), 0.5));
  assert(ApproxEqual(placedpoly1->DistanceToOut(Vec3D_t(2., 0., 0.), Vec3D_t(1., 0., 0.)), 0.));
  assert(ApproxEqual(placedpoly1->DistanceToOut(Vec3D_t(2., 0., 0.), Vec3D_t(-1., 0., 0.)), 4.));
  assert(ApproxEqual(placedpoly1->DistanceToOut(Vec3D_t(1., 0., 2), Vec3D_t(0., 0., 1.)), 0.));
  assert(ApproxEqual(placedpoly1->DistanceToOut(Vec3D_t(0.5, 0., -1), Vec3D_t(0., 0., -1.)), 0.));
  assert(ApproxEqual(placedpoly1->DistanceToOut(Vec3D_t(0.5, 0., -1), Vec3D_t(0., 0., 1.)), 3.));
  return 0;
}
