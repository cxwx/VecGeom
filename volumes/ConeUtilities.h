/*
 * ConeUtilities.h
 *
 *  Created on: June 01, 2017
 *      Author: Raman Sehgal
 */
#ifndef VECGEOM_CONEUTILITIES_H_
#define VECGEOM_CONEUTILITIES_H_

#include "base/Global.h"
#include "volumes/Wedge_Evolution.h"
#include "base/Vector3D.h"
#include "volumes/ConeStruct.h"
#include "volumes/kernel/GenericKernels.h"
#include "volumes/kernel/shapetypes/ConeTypes.h"
#include <cstdio>

#define kConeTolerance 1e-7
#define kHalfConeTolerance 0.5 * kConeTolerance

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

class UnplacedCone;
template <typename T>
struct ConeStruct;
using UnplacedStruct_t = ConeStruct<double>;

namespace ConeUtilities {
/**
* Returns whether a point is inside a cylindrical sector, as defined
* by the two vectors that go along the endpoints of the sector
*
* The same could be achieved using atan2 to calculate the angle formed
* by the point, the origin and the X-axes, but this is a lot faster,
* using only multiplications and comparisons
*
* (-x*starty + y*startx) >= 0: calculates whether going from the start vector to
*the point
* we are traveling in the CCW direction (taking the shortest direction, of
*course)
*
* (-endx*y + endy*x) >= 0: calculates whether going from the point to the end
*vector
* we are traveling in the CCW direction (taking the shortest direction, of
*course)
*
* For a sector smaller than pi, we need that BOTH of them hold true - if going
*from start, to the
* point, and then to the end we are travelling in CCW, it's obvious the point is
*inside the
* cylindrical sector.
*
* For a sector bigger than pi, only one of the conditions needs to be true. This
*is less obvious why.
* Since the sector angle is greater than pi, it can be that one of the two
*vectors might be
* farther than pi away from the point. In that case, the shortest direction will
*be CW, so even
* if the point is inside, only one of the two conditions need to hold.
*
* If going from start to point is CCW, then certainly the point is inside as the
*sector
* is larger than pi.
*
* If going from point to end is CCW, again, the point is certainly inside.
*
* This function is a frankensteinian creature that can determine which of the
*two cases (smaller vs larger than pi)
* to use either at compile time (if it has enough information, saving an ifVolumeType
*statement) or at runtime.
**/

#if (1)
template <typename Real_v, typename ShapeType, bool onSurfaceT, bool includeSurface = true>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
static void PointInCyclicalSector(UnplacedStruct_t const &volume, Real_v const &x, Real_v const &y,
                                  typename vecCore::Mask_v<Real_v> &ret)
{

  using namespace ::vecgeom::ConeTypes;
  // assert(SectorType<ShapeType>::value != kNoAngle && "ShapeType without a
  // sector passed to PointInCyclicalSector");

  // typedef Real_v Real_v;
  // using vecgeom::ConeTypes::SectorType;
  // using vecgeom::ConeTypes::EAngleType;

  Real_v startx = volume.fAlongPhi1x; // GetAlongPhi1X();
  Real_v starty = volume.fAlongPhi1y; // GetAlongPhi1Y();

  Real_v endx = volume.fAlongPhi2x; // GetAlongPhi2X();
  Real_v endy = volume.fAlongPhi2y; // GetAlongPhi2Y();

  bool smallerthanpi;

  if (SectorType<ShapeType>::value == kUnknownAngle)
    smallerthanpi = volume.fDPhi <= M_PI;
  else
    smallerthanpi = SectorType<ShapeType>::value == kPi || SectorType<ShapeType>::value == kSmallerThanPi;

  Real_v startCheck = (-x * starty) + (y * startx);
  Real_v endCheck   = (-endx * y) + (endy * x);

  if (onSurfaceT) {
    // in this case, includeSurface is irrelevant
    ret = (Abs(startCheck) <= kHalfConeTolerance) | (Abs(endCheck) <= kHalfConeTolerance);
  } else {
    if (smallerthanpi) {
      if (includeSurface)
        ret = (startCheck >= -kHalfConeTolerance) & (endCheck >= -kHalfConeTolerance);
      else
        ret = (startCheck >= kHalfConeTolerance) & (endCheck >= kHalfConeTolerance);
    } else {
      if (includeSurface)
        ret = (startCheck >= -kHalfConeTolerance) | (endCheck >= -kHalfConeTolerance);
      else
        ret = (startCheck >= kHalfConeTolerance) | (endCheck >= kHalfConeTolerance);
    }
  }
}

#endif

#if (1)
template <typename Real_v, bool ForInnerRadius>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
static Real_v GetRadiusOfConeAtPoint(UnplacedStruct_t const &cone, Real_v const pointZ)
{

  if (ForInnerRadius) {
    return cone.fInnerSlope * pointZ + cone.fInnerOffset;

  } else {
    return cone.fOuterSlope * pointZ + cone.fOuterOffset;
  }
}

#endif

/*
* Check intersection of the trajectory with a phi-plane
* All points of the along-vector of a phi plane lie on
* s * (alongX, alongY)
* All points of the trajectory of the particle lie on
* (x, y) + t * (vx, vy)
* Thefore, it must hold that s * (alongX, alongY) == (x, y) + t * (vx, vy)
* Solving by t we get t = (alongY*x - alongX*y) / (vy*alongX - vx*alongY)
* s = (x + t*vx) / alongX = (newx) / alongX
*
* If we have two non colinear phi-planes, need to make sure
* point falls on its positive direction <=> dot product between
* along vector and hit-point is positive <=> hitx*alongX + hity*alongY > 0
*/

template <typename Real_v, typename ConeType, bool PositiveDirectionOfPhiVector, bool insectorCheck>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
static void PhiPlaneTrajectoryIntersection(Precision alongX, Precision alongY, Precision normalX, Precision normalY,
                                           UnplacedStruct_t const &cone, Vector3D<Real_v> const &pos,
                                           Vector3D<Real_v> const &dir, Real_v &dist,
                                           typename vecCore::Mask_v<Real_v> &ok)
{

  // typedef Real_v Real_v;

  dist = kInfLength;

  // approaching phi plane from the right side?
  // this depends whether we use it for DistanceToIn or DistanceToOut
  // Note: wedge normals poing towards the wedge inside, by convention!
  if (insectorCheck)
    ok = ((dir.x() * normalX) + (dir.y() * normalY) > 0.); // DistToIn  -- require tracks entering volume
  else
    ok = ((dir.x() * normalX) + (dir.y() * normalY) < 0.); // DistToOut -- require tracks leaving volume

  // if( /*Backend::early_returns &&*/ vecCore::MaskEmpty(ok) ) return;

  Real_v dirDotXY = (dir.y() * alongX) - (dir.x() * alongY);
  vecCore__MaskedAssignFunc(dist, dirDotXY != 0, ((alongY * pos.x()) - (alongX * pos.y())) / NonZero(dirDotXY));
  ok &= dist > -kHalfConeTolerance;
  // if( /*Backend::early_returns &&*/ vecCore::MaskEmpty(ok) ) return;

  if (insectorCheck) {
    Real_v hitx          = pos.x() + dist * dir.x();
    Real_v hity          = pos.y() + dist * dir.y();
    Real_v hitz          = pos.z() + dist * dir.z();
    Real_v r2            = (hitx * hitx) + (hity * hity);
    Real_v innerRadIrTol = GetRadiusOfConeAtPoint<Real_v, true>(cone, hitz) + kTolerance;
    Real_v outerRadIrTol = GetRadiusOfConeAtPoint<Real_v, false>(cone, hitz) - kTolerance;

    ok &= Abs(hitz) <= cone.fTolIz && (r2 >= innerRadIrTol * innerRadIrTol) && (r2 <= outerRadIrTol * outerRadIrTol);

    // GL: tested with this if(PosDirPhiVec) around if(insector), so
    // if(insector){} requires PosDirPhiVec==true to run
    //  --> shapeTester still finishes OK (no mismatches) (some cycles saved...)
    if (PositiveDirectionOfPhiVector) {
      ok = ok && ((hitx * alongX) + (hity * alongY)) > 0.;
    }
  } else {
    if (PositiveDirectionOfPhiVector) {
      Real_v hitx = pos.x() + dist * dir.x();
      Real_v hity = pos.y() + dist * dir.y();
      ok          = ok && ((hitx * alongX) + (hity * alongY)) >= 0.;
    }
  }
}

template <typename Real_v, bool ForInnerSurface>
VECCORE_ATT_HOST_DEVICE
static Vector3D<Real_v> GetNormal(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point)
{

  // typedef Real_v Real_v;
  Real_v rho = point.Perp();
  Vector3D<Real_v> norm(0., 0., 0.);

  if (ForInnerSurface) {
    // Handling inner conical surface
    Precision rmin1 = cone.fRmin1;
    Precision rmin2 = cone.fRmin2;
    if ((rmin1 == rmin2) && (rmin1 != 0.)) {
      // cone act like tube
      norm.Set(-point.x(), -point.y(), 0.);
    } else {
      Precision secRMin = cone.fSecRMin;
      norm.Set(-point.x(), -point.y(), cone.fZNormInner * (rho * secRMin));
    }
  } else {
    Precision rmax1 = cone.fRmax1;
    Precision rmax2 = cone.fRmax2;
    if ((rmax1 == rmax2) && (rmax1 != 0.)) {
      // cone act like tube
      norm.Set(point.x(), point.y(), 0.);
    } else {
      Precision secRMax = cone.fSecRMax;
      norm.Set(point.x(), point.y(), cone.fZNormOuter * (rho * secRMax));
    }
  }
  //  std::cout<<"Normal : "<< norm << " : " << __LINE__ << std::endl;
  return norm;
}

template <typename Real_v, bool ForInnerSurface>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
static typename vecCore::Mask_v<Real_v> IsOnConicalSurface(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point)
{

  using namespace ConeUtilities;
  using namespace ConeTypes;
  // typedef Real_v Real_v;

  Real_v rho = point.Perp2();
  // std::cout<<"ForInnerSurface : " <<ForInnerSurface <<" : " << __LINE__ << std::endl;
  Real_v coneRad  = GetRadiusOfConeAtPoint<Real_v, ForInnerSurface>(cone, point.z());
  Real_v coneRad2 = coneRad * coneRad;
  // std::cout<<"rho2 : "<< rho <<" : ConeRad2 : "<< coneRad2 << std::endl;
  return (rho >= (coneRad2 - kConeTolerance * coneRad)) && (rho <= (coneRad2 + kConeTolerance * coneRad)) &&
         (Abs(point.z()) < (cone.fDz + kTolerance));
}

template <typename Real_v, bool ForInnerSurface>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
static typename vecCore::Mask_v<Real_v> IsMovingOutsideConicalSurface(UnplacedStruct_t const &cone,
                                                                      Vector3D<Real_v> const &point,
                                                                      Vector3D<Real_v> const &direction)
{

  // std::cout<<" ==== Entered IsOnConicalSurface ==== " << std::endl;
  // typedef typename vecCore::Mask_v<Real_v> Bool_t;
  // Bool_t isOnConicalSurface = IsOnConicalSurface<Backend, ForInnerSurface>(cone, point);
  // Vector3D<Real_v> norm = GetNormal<Backend, ForInnerSurface>(cone, point);
  return IsOnConicalSurface<Real_v, ForInnerSurface>(cone, point) &&
         (direction.Dot(GetNormal<Real_v, ForInnerSurface>(cone, point)) >= 0.);
}

template <typename Real_v, bool ForInnerSurface>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
static typename vecCore::Mask_v<Real_v> IsMovingInsideConicalSurface(UnplacedStruct_t const &cone,
                                                                     Vector3D<Real_v> const &point,
                                                                     Vector3D<Real_v> const &direction)
{

  //  std::cout << "IsOnConicalSurface : "<< IsOnConicalSurface<Backend, ForInnerSurface>(cone, point) << " : " <<
  //  __LINE__ << std::endl;
  return IsOnConicalSurface<Real_v, ForInnerSurface>(cone, point) &&
         (direction.Dot(GetNormal<Real_v, ForInnerSurface>(cone, point)) <= 0.);
}

template <typename Real_v>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
static typename vecCore::Mask_v<Real_v> IsOnStartPhi(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point)
{
  //  class evolution::Wedge;
  return cone.fPhiWedge.IsOnSurfaceGeneric(cone.fPhiWedge.GetAlong1(), cone.fPhiWedge.GetNormal1(), point);
}

template <typename Real_v>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
static typename vecCore::Mask_v<Real_v> IsOnEndPhi(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point)
{

  return cone.fPhiWedge.IsOnSurfaceGeneric(cone.fPhiWedge.GetAlong2(), cone.fPhiWedge.GetNormal2(), point);
}

template <typename Real_v, bool ForTopPlane>
VECCORE_ATT_HOST_DEVICE
static typename vecCore::Mask_v<Real_v> IsOnZPlaneAndMovingInside(UnplacedStruct_t const &cone,
                                                                  Vector3D<Real_v> const &point,
                                                                  Vector3D<Real_v> const &direction)
{

  Real_v rho    = point.Perp2();
  Precision fDz = cone.fDz;

  if (ForTopPlane) {
    return (rho > (cone.fSqRmin2 - kHalfConeTolerance)) && (rho < (cone.fSqRmax2 + kHalfConeTolerance)) &&
           (point.z() < (fDz + kHalfConeTolerance)) && (point.z() > (fDz - kHalfConeTolerance)) && (direction.z() < 0.);
  } else {
    return (rho > (cone.fSqRmin1 - kHalfConeTolerance)) && (rho < (cone.fSqRmax1 + kHalfConeTolerance)) &&
           (point.z() < (-fDz + kHalfConeTolerance)) && (point.z() > (-fDz - kHalfConeTolerance)) &&
           (direction.z() > 0.);
  }
}

template <typename Real_v, bool ForTopPlane>
VECCORE_ATT_HOST_DEVICE
static typename vecCore::Mask_v<Real_v> IsOnZPlaneAndMovingOutside(UnplacedStruct_t const &cone,
                                                                   Vector3D<Real_v> const &point,
                                                                   Vector3D<Real_v> const &direction)
{

  Real_v rho    = point.Perp2();
  Precision fDz = cone.fDz;

  if (ForTopPlane) {
    return (rho > (cone.fSqRmin2 - kHalfConeTolerance)) && (rho < (cone.fSqRmax2 + kHalfConeTolerance)) &&
           (point.z() < (fDz + kHalfConeTolerance)) && (point.z() > (fDz - kHalfConeTolerance)) && (direction.z() > 0.);
  } else {
    return (rho > (cone.fSqRmin1 - kHalfConeTolerance)) && (rho < (cone.fSqRmax1 + kHalfConeTolerance)) &&
           (point.z() < (-fDz + kHalfConeTolerance)) && (point.z() > (-fDz - kHalfConeTolerance)) &&
           (direction.z() < 0.);
  }
}

} // End of NS ConeUtilities

/* This class is introduced to allow Partial Specialization of selected functions,
** and will be very much useful when running Cone and Polycone in Scalar mode
*/
template <class Real_v, class coneTypeT>
class ConeHelpers {

public:
  ConeHelpers() {}
  ~ConeHelpers() {}
  template <bool ForDistToIn, bool ForInnerSurface>
  VECCORE_ATT_HOST_DEVICE
  static typename vecCore::Mask_v<Real_v> DetectIntersectionAndCalculateDistanceToConicalSurface(
      UnplacedStruct_t const &cone, Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction, Real_v &distance)
  {

    using namespace ConeUtilities;
    using namespace ConeTypes;
    // typedef Real_v Real_v;
    typedef typename vecCore::Mask_v<Real_v> Bool_t;

    Bool_t ok(false);
    Bool_t done(false);
    if (ForDistToIn) {
      Bool_t isOnSurfaceAndMovingInside =
          ConeUtilities::IsMovingInsideConicalSurface<Real_v, ForInnerSurface>(cone, point, direction);

      if (!checkPhiTreatment<coneTypeT>(cone)) {
        vecCore__MaskedAssignFunc(distance, isOnSurfaceAndMovingInside, Real_v(0.));
        done |= isOnSurfaceAndMovingInside;
        if (vecCore::MaskFull(done)) return done;
      } else {
        Bool_t insector(false);
        ConeUtilities::PointInCyclicalSector<Real_v, coneTypeT, false, true>(cone, point.x(), point.y(), insector);
        vecCore__MaskedAssignFunc(distance, insector && isOnSurfaceAndMovingInside, Real_v(0.));
        done |= (insector && isOnSurfaceAndMovingInside);
        if (vecCore::MaskFull(done)) return done;
      }

    } else {
      Bool_t isOnSurfaceAndMovingOutside =
          IsMovingOutsideConicalSurface<Real_v, ForInnerSurface>(cone, point, direction);

      if (!checkPhiTreatment<coneTypeT>(cone)) {
        vecCore__MaskedAssignFunc(distance, isOnSurfaceAndMovingOutside, Real_v(0.));
        done |= isOnSurfaceAndMovingOutside;
        if (vecCore::MaskFull(done)) return done;
      } else {
        Bool_t insector(false);
        ConeUtilities::PointInCyclicalSector<Real_v, coneTypeT, false, true>(cone, point.x(), point.y(), insector);
        vecCore__MaskedAssignFunc(distance, insector && isOnSurfaceAndMovingOutside, Real_v(0.));
        done |= (insector && isOnSurfaceAndMovingOutside);
        if (vecCore::MaskFull(done)) return done;
      }
    }

    Real_v pDotV2D = point.x() * direction.x() + point.y() * direction.y();

    Real_v a(0.), b(0.), c(0.);
    Precision fDz = cone.fDz;
    if (ForInnerSurface) {

      Precision rmin1 = cone.fRmin1;
      Precision rmin2 = cone.fRmin2;
      if (rmin1 == rmin2) {
        b = pDotV2D;
        a = direction.Perp2();
        c = point.Perp2() - rmin2 * rmin2;
      } else {

        Precision t = cone.fTanInnerApexAngle;
        Real_v newPz(0.);
        if (cone.fRmin2 > cone.fRmin1)
          newPz = (point.z() + fDz + cone.fInnerConeApex) * t;
        else
          newPz = (point.z() - fDz - cone.fInnerConeApex) * t;

        Real_v dirT = direction.z() * t;
        a           = (direction.x() * direction.x()) + (direction.y() * direction.y()) - dirT * dirT;

        b = pDotV2D - (newPz * dirT);
        c = point.Perp2() - (newPz * newPz);
      }

      Real_v b2 = b * b;
      Real_v ac = a * c;
      if (vecCore::MaskFull(b2 < ac)) return Bool_t(false);
      Real_v d2 = b2 - ac;

      Real_v delta = Sqrt(vecCore::math::Abs(d2));
      if (ForDistToIn) {
        vecCore__MaskedAssignFunc(distance, !done && d2 >= 0. && (b >= 0.), (c / NonZero(-b - delta)));
        vecCore__MaskedAssignFunc(distance, !done && d2 >= 0. && (b < 0.), (-b + delta) / NonZero(a));
      } else {
        vecCore__MaskedAssignFunc(distance, !done && d2 >= 0. && (b >= 0.), (-b - delta) / NonZero(a));
        vecCore__MaskedAssignFunc(distance, !done && d2 >= 0. && (b < 0.), (c / NonZero(-b + delta)));
      }

      if (vecCore::MaskFull(distance < 0.)) return Bool_t(false);
      Real_v newZ = point.z() + (direction.z() * distance);
      ok          = (Abs(newZ) < fDz);

    } else {

      Precision rmax1 = cone.fRmax1;
      Precision rmax2 = cone.fRmax2;
      if (rmax1 == rmax2) {
        b = pDotV2D;
        a = direction.Perp2();
        c = point.Perp2() - rmax2 * rmax2;
      } else {

        Precision t = cone.fTanOuterApexAngle;
        Real_v newPz(0.);
        if (cone.fRmax2 > cone.fRmax1)
          newPz = (point.z() + fDz + cone.fOuterConeApex) * t;
        else
          newPz     = (point.z() - fDz - cone.fOuterConeApex) * t;
        Real_v dirT = direction.z() * t;
        a           = direction.x() * direction.x() + direction.y() * direction.y() - dirT * dirT;
        b           = pDotV2D - (newPz * dirT);
        c           = point.Perp2() - (newPz * newPz);
      }
      Real_v b2 = b * b;
      Real_v ac = a * c;
      if (vecCore::MaskFull(b2 < ac)) return Bool_t(false);
      Real_v d2    = b2 - ac;
      Real_v delta = Sqrt(vecCore::math::Abs(d2));

      if (ForDistToIn) {
        vecCore__MaskedAssignFunc(distance, !done && d2 >= 0. && (b > 0.), (-b - delta) / NonZero(a));
        vecCore__MaskedAssignFunc(distance, !done && d2 >= 0. && (b < 0.), (c / NonZero(-b + delta)));
      } else {
        vecCore__MaskedAssignFunc(distance, !done && d2 >= 0. && (b < 0.), (-b + delta) / NonZero(a));
        vecCore__MaskedAssignFunc(distance, !done && d2 >= 0. && (b >= 0.), (c / NonZero(-b - delta)));
      }

      if (vecCore::MaskFull(distance < 0.)) return Bool_t(false);
      Real_v newZ = point.z() + (direction.z() * distance);
      ok          = (Abs(newZ) < fDz + kHalfConeTolerance);
    }
    vecCore__MaskedAssignFunc(distance, distance < 0., Real_v(kInfLength));

    if (checkPhiTreatment<coneTypeT>(cone)) {
      Real_v hitx(0), hity(0), hitz(0);
      Bool_t insector(false); // = Backend::kFalse;
      vecCore__MaskedAssignFunc(hitx, distance < kInfLength, point.x() + distance * direction.x());
      vecCore__MaskedAssignFunc(hity, distance < kInfLength, point.y() + distance * direction.y());
      vecCore__MaskedAssignFunc(hitz, distance < kInfLength, point.z() + distance * direction.z());

      ConeUtilities::PointInCyclicalSector<Real_v, coneTypeT, false, true>(cone, hitx, hity, insector);
      ok &= ((insector) && (distance < kInfLength));
    }
    return ok;
  }

  template <bool ForInside>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void GenericKernelForContainsAndInside(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point,
                                                typename vecCore::Mask_v<Real_v> &completelyinside,
                                                typename vecCore::Mask_v<Real_v> &completelyoutside)
  {

    typedef typename vecCore::Mask_v<Real_v> Bool_t;

    // very fast check on z-height
    Real_v absz       = Abs(point[2]);
    completelyoutside = absz > MakePlusTolerant<ForInside>(cone.fDz, kConeTolerance);
    if (ForInside) {
      completelyinside = absz < MakeMinusTolerant<ForInside>(cone.fDz, kConeTolerance);
    }
    if (vecCore::MaskFull(completelyoutside)) {
      return;
    }

    // check on RMAX
    Real_v r2 = point.x() * point.x() + point.y() * point.y();
    // calculate cone radius at the z-height of position
    Real_v rmax  = cone.fOuterSlope * point.z() + cone.fOuterOffset;
    Real_v rmax2 = rmax * rmax;

    completelyoutside |= r2 > MakePlusTolerantSquare<ForInside>(rmax, rmax2, kConeTolerance);
    if (ForInside) {
      completelyinside &= r2 < MakeMinusTolerantSquare<ForInside>(rmax, rmax2, kConeTolerance);
    }
    if (vecCore::MaskFull(completelyoutside)) {
      return;
    }

    // check on RMIN
    if (ConeTypes::checkRminTreatment<coneTypeT>(cone)) {
      Real_v rmin  = cone.fInnerSlope * point.z() + cone.fInnerOffset;
      Real_v rmin2 = rmin * rmin;

      completelyoutside |= r2 <= MakeMinusTolerantSquare<ForInside>(rmin, rmin2, kConeTolerance);
      if (ForInside) {
        completelyinside &= r2 > MakePlusTolerantSquare<ForInside>(rmin, rmin2, kConeTolerance);
      }
      if (vecCore::MaskFull(completelyoutside)) {
        return;
      }
    }

    if (ConeTypes::checkPhiTreatment<coneTypeT>(cone)) {
      Bool_t completelyoutsidephi;
      Bool_t completelyinsidephi;
      cone.fPhiWedge.GenericKernelForContainsAndInside<Real_v, ForInside>(point, completelyinsidephi,
                                                                          completelyoutsidephi);
      completelyoutside |= completelyoutsidephi;
      if (ForInside) completelyinside &= completelyinsidephi;
    }
  }

  template <typename Inside_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Inside(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point, Inside_v &inside)
  {

    using Bool_v       = vecCore::Mask_v<Real_v>;
    using InsideBool_v = vecCore::Mask_v<Inside_v>;
    Bool_v completelyinside(false), completelyoutside(false);
    GenericKernelForContainsAndInside<true>(cone, point, completelyinside, completelyoutside);
    inside = EInside::kSurface;
    vecCore::MaskedAssign(inside, (InsideBool_v)completelyoutside, Inside_v(EInside::kOutside));
    vecCore::MaskedAssign(inside, (InsideBool_v)completelyinside, Inside_v(EInside::kInside));
  }
};

template <class coneTypeT>
class ConeHelpers<Precision, coneTypeT> {

public:
  ConeHelpers() {}
  ~ConeHelpers() {}

  template <bool ForDistToIn, bool ForInnerSurface>
  VECCORE_ATT_HOST_DEVICE
  static bool DetectIntersectionAndCalculateDistanceToConicalSurface(UnplacedStruct_t const &cone,
                                                                     Vector3D<Precision> const &point,
                                                                     Vector3D<Precision> const &direction,
                                                                     Precision &distance)
  {

    using namespace ConeUtilities;
    using namespace ConeTypes;
    bool ok(false);

    if (ForDistToIn) {
      bool isOnSurfaceAndMovingInside =
          ConeUtilities::IsMovingInsideConicalSurface<Precision, ForInnerSurface>(cone, point, direction);

      if (!checkPhiTreatment<coneTypeT>(cone)) {
        if (isOnSurfaceAndMovingInside) {
          distance = 0.;
          return true;
        }
      } else {
        bool insector(false);
        ConeUtilities::PointInCyclicalSector<Precision, coneTypeT, false, true>(cone, point.x(), point.y(), insector);
        if (insector && isOnSurfaceAndMovingInside) {
          distance = 0.;
          return true;
        }
      }
    } else {
      bool isOnSurfaceAndMovingOutside =
          IsMovingOutsideConicalSurface<Precision, ForInnerSurface>(cone, point, direction);

      if (!checkPhiTreatment<coneTypeT>(cone)) {
        if (isOnSurfaceAndMovingOutside) {
          distance = 0.;
          return true;
        }
      } else {
        bool insector(false);
        ConeUtilities::PointInCyclicalSector<Precision, coneTypeT, false, true>(cone, point.x(), point.y(), insector);

        if (insector && isOnSurfaceAndMovingOutside) {
          distance = 0.;
          return true;
        }
      }
    }

    Precision pDotV2D = point.x() * direction.x() + point.y() * direction.y();

    Precision a(0.), b(0.), c(0.);
    if (ForInnerSurface) {

      if (cone.fRmin1 == cone.fRmin2) {
        b = pDotV2D;
        a = direction.Perp2();
        c = point.Perp2() - cone.fRmin2 * cone.fRmin2;
      } else {

        Precision newPz(0.);
        if (cone.fRmin2 > cone.fRmin1)
          newPz = (point.z() + cone.fDz + cone.fInnerConeApex) * cone.fTanInnerApexAngle;
        else
          newPz = (point.z() - cone.fDz - cone.fInnerConeApex) * cone.fTanInnerApexAngle;

        Precision dirT = direction.z() * cone.fTanInnerApexAngle;
        a              = (direction.x() * direction.x()) + (direction.y() * direction.y()) - dirT * dirT;

        b = pDotV2D - (newPz * dirT);
        c = point.Perp2() - (newPz * newPz);
      }

      Precision b2 = b * b;
      Precision ac = a * c;
      if (b2 < ac) return false;

      Precision d2 = b2 - ac;

      Precision delta = Sqrt(d2);
      if (ForDistToIn) {
        if (b >= 0.) {
          distance = (c / NonZero(-b - delta));
        } else {
          distance = (-b + delta) / NonZero(a);
        }
      } else {
        if (b >= 0.) {
          distance = (-b - delta) / NonZero(a);
        } else {
          distance = (c / NonZero(-b + delta));
        }
      }

      if (distance < 0.) return false;
      Precision newZ = point.z() + (direction.z() * distance);
      ok             = (Abs(newZ) < cone.fDz);

    } else {

      if (cone.fRmax1 == cone.fRmax2) {
        b = pDotV2D;
        a = direction.Perp2();
        c = point.Perp2() - cone.fRmax2 * cone.fRmax2;
      } else {

        // Precision t = cone.fTanOuterApexAngle;
        Precision newPz(0.);
        if (cone.fRmax2 > cone.fRmax1)
          newPz = (point.z() + cone.fDz + cone.fOuterConeApex) * cone.fTanOuterApexAngle;
        else
          newPz        = (point.z() - cone.fDz - cone.fOuterConeApex) * cone.fTanOuterApexAngle;
        Precision dirT = direction.z() * cone.fTanOuterApexAngle;
        a              = direction.x() * direction.x() + direction.y() * direction.y() - dirT * dirT;
        b              = pDotV2D - (newPz * dirT);
        c              = point.Perp2() - (newPz * newPz);
      }
      Precision b2 = b * b;
      Precision ac = a * c;
      if (b2 < ac) return false;
      Precision d2    = b2 - ac;
      Precision delta = Sqrt(d2);

      if (ForDistToIn) {
        if (b > 0.) {
          distance = (-b - delta) / NonZero(a); // BE ATTENTIVE, not covers the condition for b==0.
        } else {
          distance = (c / NonZero(-b + delta));
        }

      } else {
        if (b < 0.) {
          distance = (-b + delta) / NonZero(a);
        } else {
          distance = (c / NonZero(-b - delta));
        }
      }

      if (distance < 0.) return false;
      Precision newZ = point.z() + (direction.z() * distance);
      ok             = (Abs(newZ) < cone.fDz + kHalfConeTolerance);
    }
    if (distance < 0.) {
      distance = kInfLength;
    }

    if (checkPhiTreatment<coneTypeT>(cone)) {
      Precision hitx(0), hity(0); //, hitz(0);
      bool insector(false);       // = Backend::kFalse;
      if (distance < kInfLength) {
        hitx = point.x() + distance * direction.x();
        hity = point.y() + distance * direction.y();
        // hitz = point.z() + distance * direction.z();
      }

      ConeUtilities::PointInCyclicalSector<Precision, coneTypeT, false, true>(cone, hitx, hity, insector);
      ok &= ((insector) && (distance < kInfLength));
    }
    return ok;
  }

  template <bool ForInside>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void GenericKernelForContainsAndInside(UnplacedStruct_t const &cone, Vector3D<Precision> const &point,
                                                bool &completelyinside, bool &completelyoutside)
  {

    // very fast check on z-height
    Precision absz    = Abs(point[2]);
    completelyoutside = absz > MakePlusTolerant<ForInside>(cone.fDz, kConeTolerance);
    if (ForInside) {
      completelyinside = absz < MakeMinusTolerant<ForInside>(cone.fDz, kConeTolerance);
    }
    if (completelyoutside) return;

    // check on RMAX
    Precision r2 = point.x() * point.x() + point.y() * point.y();
    // calculate cone radius at the z-height of position
    Precision rmax  = cone.fOuterSlope * point.z() + cone.fOuterOffset;
    Precision rmax2 = rmax * rmax;

    completelyoutside |= r2 > MakePlusTolerantSquare<ForInside>(rmax, rmax2, kConeTolerance);
    if (ForInside) {
      completelyinside &= r2 < MakeMinusTolerantSquare<ForInside>(rmax, rmax2, kConeTolerance);
    }
    if (completelyoutside) return;

    // check on RMIN
    if (ConeTypes::checkRminTreatment<coneTypeT>(cone)) {
      Precision rmin  = cone.fInnerSlope * point.z() + cone.fInnerOffset;
      Precision rmin2 = rmin * rmin;

      completelyoutside |= r2 <= MakeMinusTolerantSquare<ForInside>(rmin, rmin2, kConeTolerance);
      if (ForInside) {
        completelyinside &= r2 > MakePlusTolerantSquare<ForInside>(rmin, rmin2, kConeTolerance);
      }
      if (completelyoutside) return;
    }

    if (ConeTypes::checkPhiTreatment<coneTypeT>(cone)) {
      bool completelyoutsidephi;
      bool completelyinsidephi;
      cone.fPhiWedge.GenericKernelForContainsAndInside<Precision, ForInside>(point, completelyinsidephi,
                                                                             completelyoutsidephi);
      completelyoutside |= completelyoutsidephi;
      if (ForInside) completelyinside &= completelyinsidephi;
    }
  }

  template <typename Inside_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Inside(UnplacedStruct_t const &cone, Vector3D<Precision> const &point, Inside_v &inside)
  {
    bool completelyinside(false), completelyoutside(false);
    GenericKernelForContainsAndInside<true>(cone, point, completelyinside, completelyoutside);
    inside                        = EInside::kSurface;
    if (completelyoutside) inside = EInside::kOutside;
    if (completelyinside) inside  = EInside::kInside;
  }
};
}
}

#endif