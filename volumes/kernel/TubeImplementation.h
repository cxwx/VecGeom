/// @file TubeImplementation.h
/// @author Georgios Bitzes (georgios.bitzes@cern.ch)

#ifndef VECGEOM_VOLUMES_KERNEL_TUBEIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_TUBEIMPLEMENTATION_H_

#include "base/Global.h"
#include "base/Transformation3D.h"
#include "base/Vector3D.h"
#include "volumes/kernel/GenericKernels.h"
#include "volumes/kernel/shapetypes/TubeTypes.h"
#include "volumes/UnplacedTube.h"
#include "volumes/Wedge.h"
#include <cstdio>

#define TUBE_SAFETY_OLD  // use old (and faster) definitions of SafetyToIn() and SafetyToOut()

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v_1t(TubeImplementation, TranslationCode, translation::kGeneric, RotationCode, rotation::kGeneric, typename)

inline namespace VECGEOM_IMPL_NAMESPACE {

namespace TubeUtilities {

/**
 * Returns whether a point is inside a cylindrical sector, as defined
 * by the two vectors that go along the endpoints of the sector
 *
 * The same could be achieved using atan2 to calculate the angle formed
 * by the point, the origin and the X-axes, but this is a lot faster,
 * using only multiplications and comparisons
 *
 * (-x*starty + y*startx) >= 0: calculates whether going from the start vector to the point
 * we are traveling in the CCW direction (taking the shortest direction, of course)
 *
 * (-endx*y + endy*x) >= 0: calculates whether going from the point to the end vector
 * we are traveling in the CCW direction (taking the shortest direction, of course)
 *
 * For a sector smaller than pi, we need that BOTH of them hold true - if going from start, to the
 * point, and then to the end we are travelling in CCW, it's obvious the point is inside the
 * cylindrical sector. 
 * 
 * For a sector bigger than pi, only one of the conditions needs to be true. This is less obvious why.
 * Since the sector angle is greater than pi, it can be that one of the two vectors might be
 * farther than pi away from the point. In that case, the shortest direction will be CW, so even
 * if the point is inside, only one of the two conditions need to hold.
 * 
 * If going from start to point is CCW, then certainly the point is inside as the sector
 * is larger than pi.
 *
 * If going from point to end is CCW, again, the point is certainly inside.
 *
 * This function is a frankensteinian creature that can determine which of the two cases (smaller vs larger than pi)
 * to use either at compile time (if it has enough information, saving an if statement) or at runtime.
**/

  template<typename Backend, typename ShapeType, typename UnplacedVolumeType, bool onSurfaceT, bool includeSurface=true>
 VECGEOM_INLINE
 VECGEOM_CUDA_HEADER_BOTH
 void PointInCyclicalSector(UnplacedVolumeType const& volume, typename Backend::precision_v const &x, 
   typename Backend::precision_v const &y, typename Backend::bool_v &ret) {
  using namespace TubeTypes;
  // assert(SectorType<ShapeType>::value != kNoAngle && "ShapeType without a sector passed to PointInCyclicalSector");

  typedef typename Backend::precision_v Float_t;

  Float_t startx = volume.alongPhi1x();
  Float_t starty = volume.alongPhi1y();

  Float_t endx = volume.alongPhi2x();
  Float_t endy = volume.alongPhi2y();

  bool smallerthanpi;

  if(SectorType<ShapeType>::value == kUnknownAngle)
    smallerthanpi = volume.dphi() <= M_PI;
  else
    smallerthanpi = SectorType<ShapeType>::value == kOnePi || SectorType<ShapeType>::value == kSmallerThanPi;

  Float_t startCheck = (-x*starty + y*startx);
  Float_t endCheck = (-endx*y + endy*x);

  if(onSurfaceT) {
    //in this case, includeSurface is irrelevant
    ret = (Abs(startCheck) <= kHalfTolerance) | (Abs(endCheck) <= kHalfTolerance);
  }
  else {
    if(smallerthanpi) {
      if(includeSurface) ret = (startCheck >= -kHalfTolerance) & (endCheck >= -kHalfTolerance);
      else               ret = (startCheck >=  kHalfTolerance) & (endCheck >=  kHalfTolerance);
    }
    else {
      if(includeSurface) ret = (startCheck >= -kHalfTolerance) | (endCheck >= -kHalfTolerance);
      else               ret = (startCheck >=  kHalfTolerance) | (endCheck >=  kHalfTolerance);
    }
  }
}

template<typename Backend, typename TubeType, bool LargestSolution, bool insectorCheck>
VECGEOM_INLINE
VECGEOM_CUDA_HEADER_BOTH
void CircleTrajectoryIntersection(typename Backend::precision_v const &b,
                                  typename Backend::precision_v const &c,
                                  UnplacedTube const& tube,
                                  Vector3D<typename Backend::precision_v> const &pos,
                                  Vector3D<typename Backend::precision_v> const &dir,
                                  typename Backend::precision_v &dist,
                                  typename Backend::bool_v &ok) {
  using namespace TubeTypes;
  typedef typename Backend::precision_v Float_t;
  typedef typename Backend::bool_v Bool_t;

  Float_t delta = b*b - c;
  ok = delta > 0.0;
  if(LargestSolution)  ok |= delta == 0.;  // this takes care of scratching conventions

  MaskedAssign( !ok, 0., &delta);
  delta = Sqrt(delta);
  if(!LargestSolution) delta = -delta;

  dist = -b + delta;
  ok &= dist >= -kHalfTolerance;
  if(Backend::early_returns && IsEmpty(ok)) return;

  if(insectorCheck) {
    Float_t hitz = pos.z() + dist * dir.z();
    ok &= (Abs(hitz) <= tube.z());
    if(Backend::early_returns && IsEmpty(ok)) return;

    if(checkPhiTreatment<TubeType>(tube)) {
      Bool_t insector = Backend::kFalse;
      Float_t hitx = pos.x() + dist * dir.x();
      Float_t hity = pos.y() + dist * dir.y();
      PointInCyclicalSector<Backend, TubeType, UnplacedTube, false, true>(tube, hitx, hity, insector);
      //insector = tube.GetWedge().ContainsWithBoundary<Backend>( Vector3D<typename Backend::precision_v>(hitx, hity, hitz) );
      ok &= insector;
    }
  }
}

/*
 * Input: A point p and a unit vector v.
 * Returns the perpendicular distance between
 * (the infinite line defined by the vector) 
 * and (the point)
 * 
 * How does it work? Let phi be the angle formed
 * by v and the position vector of the point.
 * 
 * Let proj be the projection vector of point onto that
 * line.
 * 
 * We now have a right triangle formed by points
 * (0, 0), point and the projected point
 * 
 * For that triangle, it holds that
 * sin theta = perpendiculardistance / (magnitude of point vector)
 * 
 * So perpendiculardistance = sin theta * (magnitude of point vector)
 *
 * But.. the magnitude of the cross product between the point vector
 * and the v vector is:
 *
 * |p x v| = |p| * |v| * sin theta
 * 
 * Since |v| = 1, the magnitude of the cross product is exactly
 * what we're looking for, the formula for which is simply
 * p.x * v.y - p.y * v.x
 *
 */

template<typename Backend>
VECGEOM_INLINE
VECGEOM_CUDA_HEADER_BOTH
void PerpendicularDistance2D(typename Backend::precision_v const &px, typename Backend::precision_v const &py,
                            typename Backend::precision_v const &vx, typename Backend::precision_v const &vy,
                            typename Backend::precision_v &dist) {
  dist = px * vy - py * vx;
}

/*
 * Find safety distance from a point to the phi plane
 */

template<typename Backend, typename TubeType, bool inside>
VECGEOM_INLINE
VECGEOM_CUDA_HEADER_BOTH
void PhiPlaneSafety( UnplacedTube const& tube,
                     Vector3D<typename Backend::precision_v> const &pos,
                     typename Backend::precision_v &safety) {
  using namespace TubeTypes;
  typedef typename Backend::precision_v Float_t;

  safety = kInfinity;
  if(SectorType<TubeType>::value == kUnknownAngle) {
    // for phi-sections w/dphi>PI, the right result could be distance to origin
    MaskedAssign( tube.dphi() > M_PI, Sqrt(pos.x()*pos.x() + pos.y()*pos.y()), &safety);
  }
  else if( SectorType<TubeType>::value == kBiggerThanPi ) {
    safety = Sqrt( pos.x()*pos.x() + pos.y()*pos.y() );
  }

  Float_t phi1;
  PerpendicularDistance2D<Backend>(pos.x(), pos.y(), tube.alongPhi1x(), tube.alongPhi1y(), phi1);
  if(inside)
    phi1 *= -1;

  if(SectorType<TubeType>::value == kOnePi) {
    MaskedAssign( Abs(phi1)>kHalfTolerance, Abs(phi1), &safety );
    return;
  }

  // make sure point falls on positive part of projection
  MaskedAssign(phi1>-kHalfTolerance
               && pos.x()*tube.alongPhi1x()+pos.y()*tube.alongPhi1y() > 0
               && phi1<safety, phi1, &safety);

  Float_t phi2;
  PerpendicularDistance2D<Backend>(pos.x(), pos.y(), tube.alongPhi2x(), tube.alongPhi2y(), phi2);
  if(!inside)
    phi2 *= -1;

  // make sure point falls on positive part of projection
  MaskedAssign(phi2>-kHalfTolerance
               && pos.x()*tube.alongPhi2x()+pos.y()*tube.alongPhi2y() > 0
               && phi2<safety, phi2, &safety);

}

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

template<typename Backend, typename TubeType, bool PositiveDirectionOfPhiVector, bool insectorCheck>
VECGEOM_INLINE
VECGEOM_CUDA_HEADER_BOTH
void PhiPlaneTrajectoryIntersection(Precision alongX, Precision alongY,
                                    Precision normalX, Precision normalY,
                                    UnplacedTube const& tube,
                                    Vector3D<typename Backend::precision_v> const &pos,
                                    Vector3D<typename Backend::precision_v> const &dir,
                                    typename Backend::precision_v &dist,
                                    typename Backend::bool_v &ok) {

  typedef typename Backend::precision_v Float_t;
  dist = kInfinity;

  // approaching phi plane from the right side?
  // this depends whether we use it for DistanceToIn or DistanceToOut
  // Note: wedge normals poing towards the wedge inside, by convention!
  if(insectorCheck) ok = ( dir.x()*normalX + dir.y()*normalY > 0. );  // DistToIn  -- require tracks entering volume
  else              ok = ( dir.x()*normalX + dir.y()*normalY < 0. );  // DistToOut -- require tracks leaving volume

  // if( Backend::early_returns && IsEmpty(ok) ) return;

  Float_t dirDotXY = (dir.y()*alongX - dir.x()*alongY);
  MaskedAssign( dirDotXY!=0, (alongY*pos.x() - alongX*pos.y() ) / dirDotXY, &dist );
  ok &= dist > -kHalfTolerance;
  // if( Backend::early_returns && IsEmpty(ok) ) return;

  if(insectorCheck) {
    Float_t hitx = pos.x() + dist * dir.x();
    Float_t hity = pos.y() + dist * dir.y();
    Float_t hitz = pos.z() + dist * dir.z();
    Float_t r2 = hitx*hitx + hity*hity;
    ok &= Abs(hitz) <= tube.tolIz()
       && (r2 >= tube.tolIrmin2())
       && (r2 <= tube.tolIrmax2());

    // GL: tested with this if(PosDirPhiVec) around if(insector), so if(insector){} requires PosDirPhiVec==true to run
    //  --> shapeTester still finishes OK (no mismatches) (some cycles saved...)
    if(PositiveDirectionOfPhiVector) {
      ok = ok && (hitx*alongX + hity*alongY) > 0.;
    }
  }
  else {
    if(PositiveDirectionOfPhiVector) {
      Float_t hitx = pos.x() + dist * dir.x();
      Float_t hity = pos.y() + dist * dir.y();
      ok = ok && (hitx*alongX + hity*alongY) >= 0.;
    }
  }
}

} // End of NS TubeUtilities

class PlacedTube;

template <TranslationCode transCodeT, RotationCode rotCodeT, typename tubeTypeT>
struct TubeImplementation {

  static const int transC = transCodeT;
  static const int rotC   = rotCodeT;
  using PlacedShape_t = PlacedTube;
  using UnplacedShape_t = UnplacedTube;

  VECGEOM_CUDA_HEADER_BOTH
  static void PrintType() {
     printf("SpecializedTube<%i, %i, %s>", transCodeT, rotCodeT, tubeTypeT::toString());
  }

  template <typename Stream>
  static void PrintType( Stream & s ) {
       s << "SpecializedTube<" << transCodeT << "," << rotCodeT << ",TubeTypes::" << tubeTypeT::toString() << ">";
    }

  template <typename Stream>
  static void PrintImplementationType(Stream &s) {
    s << "TubeImplementation<" << transCodeT << "," << rotCodeT << ",TubeTypes::" << tubeTypeT::toString() << ">";
  }

  template <typename Stream>
  static void PrintUnplacedType(Stream &s) { s << "UnplacedTube"; }

  /////GenericKernel Contains/Inside implementation
    template <typename Backend, bool ForInside>
    VECGEOM_INLINE
    VECGEOM_CUDA_HEADER_BOTH
    static void GenericKernelForContainsAndInside(UnplacedTube const &tube,
              Vector3D<typename Backend::precision_v> const &point,
              typename Backend::bool_v &completelyinside,
              typename Backend::bool_v &completelyoutside) {

        typedef typename Backend::precision_v Float_t;
        typedef typename Backend::bool_v Bool_t;

        // very fast check on z-height
        Float_t absz = Abs(point[2]);
        completelyoutside = absz > MakePlusTolerant<ForInside>( tube.z() );
        if (ForInside)
        {
            completelyinside = absz < MakeMinusTolerant<ForInside>( tube.z() );
        }
        if (Backend::early_returns) {
            if ( IsFull(completelyoutside) ) {
              return;
            }
        }

        // check on RMAX
        Float_t r2 = point.x()*point.x()+point.y()*point.y();
        // calculate cone radius at the z-height of position

        completelyoutside |= r2 > MakePlusTolerantSquare<ForInside>( tube.rmax(), tube.rmax2() );
        if (ForInside)
        {
          completelyinside &= r2 < MakeMinusTolerantSquare<ForInside>( tube.rmax(), tube.rmax2() );
        }
        if (Backend::early_returns) {
                if ( IsFull(completelyoutside) ) {
                  return;
                }
        }

        // check on RMIN
        if (TubeTypes::checkRminTreatment<tubeTypeT>(tube)) {
          completelyoutside |= r2 <= MakeMinusTolerantSquare<ForInside>( tube.rmin(), tube.rmin2() );
          if (ForInside)
          {
           completelyinside &= r2 > MakePlusTolerantSquare<ForInside>( tube.rmin(), tube.rmin2() );
          }
            if (Backend::early_returns) {
               if ( IsFull(completelyoutside) ) {
                  return;
               }
            }
        }

        if(TubeTypes::checkPhiTreatment<tubeTypeT>(tube)) {
              Bool_t completelyoutsidephi;
              Bool_t completelyinsidephi;
              tube.GetWedge().GenericKernelForContainsAndInside<Backend,ForInside>( point,
                completelyinsidephi, completelyoutsidephi );

              completelyoutside |= completelyoutsidephi;
              if( ForInside )
                completelyinside &= completelyinsidephi;
           }
    }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void UnplacedContains(
      UnplacedTube const &tube,
      Vector3D<typename Backend::precision_v> const &point,
      typename Backend::bool_v &contains) {

      typedef typename Backend::bool_v Bool_t;
      Bool_t unused;
      Bool_t outside;
      GenericKernelForContainsAndInside<Backend, false>(tube, point, unused, outside);
      contains = !outside;
  }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void Contains(
      UnplacedTube const &tube,
      Transformation3D const &transformation,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> &localPoint,
      typename Backend::bool_v &inside) {

    localPoint = transformation.Transform<transCodeT, rotCodeT>(point);
    UnplacedContains<Backend>(tube, localPoint, inside);
  }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void Inside(UnplacedTube const &tube,
                     Transformation3D const &transformation,
                     Vector3D<typename Backend::precision_v> const &point,
                     typename Backend::inside_v &inside) {

      Vector3D<typename Backend::precision_v> localPoint
            = transformation.Transform<transCodeT, rotCodeT>(point);

      typedef typename Backend::bool_v Bool_t;
      Bool_t completelyinside, completelyoutside;
      GenericKernelForContainsAndInside<Backend,true>(tube, localPoint, completelyinside, completelyoutside);
      inside = EInside::kSurface;
      MaskedAssign(completelyoutside, EInside::kOutside, &inside);
      MaskedAssign(completelyinside,  EInside::kInside, &inside);
  }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void DistanceToIn(UnplacedTube const &tube,
                           Transformation3D const &transformation,
                           Vector3D<typename Backend::precision_v> const &masterPoint,
                           Vector3D<typename Backend::precision_v> const &masterDirection,
                           typename Backend::precision_v const &/*stepMax*/,
                           typename Backend::precision_v &distance)
  {
    using namespace TubeUtilities;
    using namespace TubeTypes;  
    typedef typename Backend::precision_v Float_t;
    typedef typename Backend::bool_v Bool_t;

    Vector3D<Float_t> point;
    Vector3D<Float_t> dir;
    transformation.Transform<transCodeT, rotCodeT>(masterPoint, point);
    transformation.TransformDirection<rotCodeT>(masterDirection, dir);

    Bool_t done = Backend::kFalse;

    //=== First, for points outside and moving away --> return infinity
    distance = kInfinity;

    // outside of Z range and going away?
    Float_t distz = Abs(point.z()) - tube.z();   // avoid a division for now
    done |= distz>kHalfTolerance && point.z()*dir.z() >= 0;

    // // outside of tube and going away?
    // done |= Abs(point.x()) > tube.rmax()+kHalfTolerance && point.x()*dir.x() >= 0;
    // done |= Abs(point.y()) > tube.rmax()+kHalfTolerance && point.y()*dir.y() >= 0;
    // if(Backend::early_returns && IsFull(done)) return;

    // outside of outer tube and going away?
    Float_t rsq = point.x()*point.x() + point.y()*point.y();
    Float_t rdotn = point.x()*dir.x() + point.y()*dir.y();
    done |= rsq > tube.tolIrmax2() && rdotn>=0;
    if( Backend::early_returns && IsFull(done) ) return;


    //=== Next, check all dimensions of the tube, whether points are inside --> return -1
    MaskedAssign( !done, -1.0, &distance );

    // For points inside z-range, return -1
    Bool_t inside = distz < -kHalfTolerance;

    inside &= rsq < tube.tolIrmax2();
    if(checkRminTreatment<tubeTypeT>(tube)) {
      inside &= rsq > tube.tolIrmin2();
    }
    if(checkPhiTreatment<tubeTypeT>(tube) && !IsEmpty(inside)) {
      Bool_t insector;
      PointInCyclicalSector<Backend, tubeTypeT, UnplacedTube, false, false>(tube, point.x(), point.y(), insector);
      inside &= insector;
      //inside &= tube.GetWedge().ContainsWithoutBoundary<Backend>( point );  // slower than PointInCyclicalSector()
    }
    done |= inside;
    if(Backend::early_returns && IsFull(done)) return;


    //=== Next step: check if z-plane is the right entry point (both r,phi should be valid at z-plane crossing)
    MaskedAssign( !done, kInfinity, &distance );

    distz /= Abs(dir.z());

    Float_t hitx = point.x() + distz*dir.x();
    Float_t hity = point.y() + distz*dir.y();
    Float_t r2 = hitx*hitx + hity*hity;  // radius of intersection with z-plane
    Bool_t okz = distz > -kHalfTolerance  && (point.z()*dir.z()<0);

    okz &= (r2 <= tube.rmax2());
    if(checkRminTreatment<tubeTypeT>(tube)) {
      okz &= (tube.rmin2() <= r2);
    }
    if(checkPhiTreatment<tubeTypeT>(tube) && !IsEmpty(okz)) {
      Bool_t insector;
      PointInCyclicalSector<Backend, tubeTypeT, UnplacedTube, false>(tube, hitx, hity, insector);
      okz &= insector;
      //okz &= tube.GetWedge().ContainsWithBoundary<Backend>(  Vector3D<Float_t>(hitx, hity, 0.0) );
    }
    MaskedAssign( !done && okz, distz, &distance);
    done |= okz;
    //if(Backend::early_returns && IsFull(done) ) return;


    //=== Next step: intersection of the trajectories with the two circles

    // Here for values used in both rmin and rmax calculations
    Float_t invnsq = 1.0 / ( 1.0 - dir.z()*dir.z() );
    Float_t b = invnsq * rdotn;

    /*
     * rmax
     * If the particle were to hit rmax, it would hit the closest point of the two
     * --> only consider the smallest solution of the quadratic equation
     */
    Float_t crmax = invnsq * (rsq - tube.rmax2());
    Float_t dist_rmax;
    Bool_t ok_rmax = Backend::kFalse;
    CircleTrajectoryIntersection<Backend, tubeTypeT, false, true>(b, crmax, tube, point, dir, dist_rmax, ok_rmax);
    ok_rmax &= dist_rmax<distance;
    MaskedAssign( !done && ok_rmax, dist_rmax, &distance);
    done |= ok_rmax;
    if(Backend::early_returns && IsFull(done)) return;

    /*
     * rmin
     * If the particle were to hit rmin, it would hit the farthest point of the two
     * --> only consider the largest solution to the quadratic equation
     */
    Float_t dist_rmin = -kInfinity;
    Bool_t ok_rmin = Backend::kFalse;
    if(checkRminTreatment<tubeTypeT>(tube)) {
      /*
       * What happens if both intersections are valid for the same particle?
       * This can only happen when particle is outside of the hollow space and will certainly hit rmax, not rmin
       * So rmax solution always takes priority over rmin, and will overwrite it in case both are valid
       */
      Float_t crmin = invnsq * (rsq - tube.rmin2());
      CircleTrajectoryIntersection<Backend, tubeTypeT, true, true>(b, crmin, tube, point, dir, dist_rmin, ok_rmin);
      ok_rmin &= dist_rmin<distance;
      MaskedAssign(!done && ok_rmin, dist_rmin, &distance);
      // done |= ok_rmin; // can't be done here, it's wrong in case phi-treatment is needed!
    }


    /*
     * Calculate intersection between trajectory and the two phi planes
     */
    if(checkPhiTreatment<tubeTypeT>(tube)) {

      Float_t dist_phi;
      Bool_t ok_phi;
      Wedge const& w = tube.GetWedge();
      PhiPlaneTrajectoryIntersection<Backend, tubeTypeT, SectorType<tubeTypeT>::value != kOnePi, true>(
              tube.alongPhi1x(), tube.alongPhi1y(),
              w.GetNormal1().x(), w.GetNormal1().y(),
              tube, point, dir, dist_phi, ok_phi);
      ok_phi &= dist_phi<distance;
      MaskedAssign(!done && ok_phi, dist_phi, &distance);
      done |= ok_phi;
//      if(Backend::early_returns && IsFull(done)) return;

      /*
       * If the tube is pi degrees, there's just one phi plane,
       * so no need to check again
       */

      if(SectorType<tubeTypeT>::value != kOnePi) {
        PhiPlaneTrajectoryIntersection<Backend, tubeTypeT, true, true>(
                tube.alongPhi2x(), tube.alongPhi2y(),
                w.GetNormal2().x(), w.GetNormal2().y(),
                tube, point, dir, dist_phi, ok_phi);
        MaskedAssign(ok_phi && dist_phi<distance, dist_phi, &distance);
      }
    }
  } // end of DistanceToIn()


  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void DistanceToOut(
      UnplacedTube const &tube,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> const &dir,
      typename Backend::precision_v const &/*stepMax*/,
      typename Backend::precision_v &distance) {

    using namespace TubeTypes;
    using namespace TubeUtilities;
    typedef typename Backend::precision_v Float_t;
    typedef typename Backend::bool_v Bool_t;

    distance = -1;
    Bool_t done = Backend::kFalse;

    //=== First we check all dimensions of the tube, whether points are outside --> return -1

    // For points outside z-range, return -1
    Float_t distz = tube.z() - Abs(point.z());   // avoid a division for now
    done |= distz<-kHalfTolerance;  // distance is already set to -1
//  if(Backend::early_returns && IsFull(done)) return;

    Float_t rsq = point.x()*point.x() + point.y()*point.y();
    Float_t rdotn = dir.x()*point.x() + dir.y()*point.y();
    Float_t crmax = rsq - tube.rmax2();  // avoid a division for now
    Float_t crmin = rsq;

    // if outside of Rmax, return -1
    done |= crmax > kTolerance*tube.rmax();
    if( Backend::early_returns && IsFull(done) ) return;

    if(checkRminTreatment<tubeTypeT>(tube)) {
      // if point is within inner-hole of a hollow tube, it is outside of the tube --> return -1
      crmin -= tube.rmin2();  // avoid a division for now
      done |=  crmin < -kTolerance*tube.rmin();
      if( Backend::early_returns && IsFull(done) ) return;
    }

    // TODO: add outside check for phi-sections here

    // OK, since we're here, then distance must be non-negative, and the smallest of possible intersections
    MaskedAssign( !done, kInfinity, &distance );

    Float_t invdirz = 1./dir.z();
    distz = (tube.z() - point.z()) * invdirz;
    MaskedAssign(dir.z() < 0, (-tube.z() - point.z()) * invdirz, &distz);
    MaskedAssign( !done && distz<distance, distz, &distance );


    /*
     * Find the intersection of the trajectories with the two circles.
     * Here I compute values used in both rmin and rmax calculations.
     */

    Float_t invnsq = 1 / ( dir.x()*dir.x() + dir.y()*dir.y() );
    Float_t b = invnsq * rdotn;

    /*
     * rmin
     */

    if(checkRminTreatment<tubeTypeT>(tube)) {
      Float_t dist_rmin = kInfinity;
      Bool_t ok_rmin = Backend::kFalse;
      crmin *= invnsq;
      CircleTrajectoryIntersection<Backend, tubeTypeT, false, false>(b, crmin, tube, point, dir, dist_rmin, ok_rmin);
      MaskedAssign(ok_rmin && dist_rmin < distance, dist_rmin, &distance);
    }

    /*
     * rmax
     */

    Float_t dist_rmax = kInfinity;
    Bool_t ok_rmax = Backend::kFalse;
    crmax *= invnsq;
    CircleTrajectoryIntersection<Backend, tubeTypeT, true, false>(b, crmax, tube, point, dir, dist_rmax, ok_rmax);
    MaskedAssign(ok_rmax && dist_rmax < distance, dist_rmax, &distance);


    /* Phi planes
     *
     * OK, this is getting weird - the only time I need to
     * check if hit-point falls on the positive direction
     * of the phi-vector is when angle is bigger than PI.
     *
     * Otherwise, any distance I get from there is guaranteed to
     * be larger - so final result would still be correct and no need to
     * check it
     */

    if(checkPhiTreatment<tubeTypeT>(tube)) {
      Float_t dist_phi = kInfinity;
      Bool_t ok_phi = Backend::kFalse;

      Wedge const& w = tube.GetWedge();
      if(SectorType<tubeTypeT>::value == kSmallerThanPi) {

        Precision normal1X = w.GetNormal1().x();
        Precision normal1Y = w.GetNormal1().y();
        PhiPlaneTrajectoryIntersection<Backend, tubeTypeT, false, false>(
          tube.alongPhi1x(), tube.alongPhi1y(), normal1X, normal1Y, tube, point, dir, dist_phi, ok_phi);
        MaskedAssign(ok_phi && dist_phi < distance, dist_phi, &distance);

        PhiPlaneTrajectoryIntersection<Backend, tubeTypeT, false, false>(
          tube.alongPhi2x(), tube.alongPhi2y(), w.GetNormal2().x(), w.GetNormal2().y(),
          tube, point, dir, dist_phi, ok_phi
          );
        MaskedAssign(ok_phi && dist_phi < distance, dist_phi, &distance);
      }
      else if(SectorType<tubeTypeT>::value == kOnePi) {
        PhiPlaneTrajectoryIntersection<Backend, tubeTypeT, false, false>(
          tube.alongPhi2x(), tube.alongPhi2y(),
          w.GetNormal2().x(), w.GetNormal2().x(),
          tube, point, dir, dist_phi, ok_phi);
        MaskedAssign(ok_phi && dist_phi < distance, dist_phi, &distance);
      }
      else {
        // angle bigger than pi or unknown
        // need to check that point falls on positive direction of phi-vectors
        PhiPlaneTrajectoryIntersection<Backend, tubeTypeT, true, false>(
            tube.alongPhi1x(), tube.alongPhi1y(), w.GetNormal1().x(), w.GetNormal1().y(), tube, point, dir, dist_phi, ok_phi);
        MaskedAssign(ok_phi && dist_phi < distance, dist_phi, &distance);

        PhiPlaneTrajectoryIntersection<Backend, tubeTypeT, true, false>(
          tube.alongPhi2x(), tube.alongPhi2y(), w.GetNormal2().x(), w.GetNormal2().y(), tube, point, dir, dist_phi, ok_phi);
        MaskedAssign(ok_phi && dist_phi < distance, dist_phi, &distance);
      }
    }
  }

  /// This function keeps track of both positive (outside) and negative (inside) distances separately
  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyAssign(typename Backend::precision_v safety,
                    typename Backend::precision_v& positiveSafety,
                    typename Backend::precision_v& negativeSafety) {
    MaskedAssign( safety>=0 && safety<positiveSafety, safety, &positiveSafety );
    MaskedAssign( safety<=0 && safety>negativeSafety, safety, &negativeSafety );
  }

/** SafetyKernel finds distances from point to each face of the tube, returning
    largest negative distance (w.r.t. faces which point is inside of ) and
    smallest positive distance (w.r.t. faces which point is outside of)
 */
  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyKernel(UnplacedTube const &tube,
                           Vector3D<typename Backend::precision_v> const &point,
                           typename Backend::precision_v &safePos,
                           typename Backend::precision_v &safeNeg) {

    // TODO: implement caching if input point is not changed
    using namespace TubeTypes;
    using namespace TubeUtilities;
    typedef typename Backend::precision_v Float_t;
    // typedef typename Backend::bool_v Bool_t;

    safePos = kInfinity;
    safeNeg = -safePos;   // reuse to avoid casting overhead

    Float_t safez = Abs(point.z()) - tube.z();
    SafetyAssign<Backend>( safez, safePos, safeNeg );

    Float_t r = Sqrt(point.x()*point.x() + point.y()*point.y());
    Float_t safermax = r - tube.rmax();
    SafetyAssign<Backend>( safermax, safePos, safeNeg );

    if(checkRminTreatment<tubeTypeT>(tube)) {
      Float_t safermin = tube.rmin() - r;
      SafetyAssign<Backend>(safermin, safePos, safeNeg);
    }

    if(checkPhiTreatment<tubeTypeT>(tube)) {
      Float_t safephi;
      PhiPlaneSafety<Backend, tubeTypeT, false>(tube, point, safephi);
      SafetyAssign<Backend>(safephi, safePos, safeNeg);
    }
  }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyToIn(UnplacedTube const &tube,
                         Transformation3D const &transformation,
                         Vector3D<typename Backend::precision_v> const &point,
                         typename Backend::precision_v &safety)
  {

#ifdef TUBE_SAFETY_OLD
    SafetyToInOld<Backend>(tube, transformation, point, safety);
#else
    typedef typename Backend::precision_v Float_t;

    Vector3D<Float_t> local_point = transformation.Transform<transCodeT,rotCodeT>(point);

    Float_t safetyInsidePoint, safetyOutsidePoint;
    SafetyKernel<Backend>(tube, local_point, safetyOutsidePoint, safetyInsidePoint);

    // Mostly called for points outside --> safetyOutside is finite --> return safetyOutside
    // If safetyOutside == infinity --> return safetyInside
    CondAssign( safetyOutsidePoint==kInfinity, safetyInsidePoint, safetyOutsidePoint, &safety );
#endif
  }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyToOut(UnplacedTube const &tube,
                          Vector3D<typename Backend::precision_v> const &point,
                          typename Backend::precision_v &safety)
  {

#ifdef TUBE_SAFETY_OLD
    SafetyToOutOld<Backend>(tube, point, safety);
#else
    typedef typename Backend::precision_v Float_t;

    Float_t safetyInsidePoint, safetyOutsidePoint;
    SafetyKernel<Backend>(tube, point, safetyOutsidePoint, safetyInsidePoint);

    // Mostly called for points inside --> safetyOutside==infinity, return |safetyInside| (flip sign)
    // If called for points outside -- return -safetyOutside
    CondAssign( safetyOutsidePoint==kInfinity, -safetyInsidePoint, -safetyOutsidePoint, &safety );
#endif
  }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyToInOld(UnplacedTube const &tube,
                            Transformation3D const &transformation,
                            Vector3D<typename Backend::precision_v> const &motherPoint,
                            typename Backend::precision_v &safety)
  {
    using namespace TubeTypes;
    using namespace TubeUtilities;
    typedef typename Backend::precision_v Float_t;
    typedef typename Backend::bool_v Bool_t;

    Vector3D<Float_t> point = transformation.Transform<transCodeT,rotCodeT>(motherPoint);

    safety = Abs(point.z()) - tube.z();

    Float_t r = Sqrt(point.x()*point.x() + point.y()*point.y());
    Float_t safermax = r - tube.rmax();
    MaskedAssign( safermax>safety, safermax, &safety );

    if(checkRminTreatment<tubeTypeT>(tube)) {
      Float_t safermin = tube.rmin() - r;
      MaskedAssign( safermin>safety, safermin, &safety);
    }

    if(checkPhiTreatment<tubeTypeT>(tube)) {
      Bool_t insector;
      PointInCyclicalSector<Backend, tubeTypeT, UnplacedTube, false, false>(tube, point.x(), point.y(), insector);
      if(Backend::early_returns && IsFull(insector)) return;

      Float_t safephi;
      PhiPlaneSafety<Backend, tubeTypeT, false>(tube, point, safephi);
      MaskedAssign(!insector && safephi<kInfinity && safephi>safety, safephi, &safety );
    }
  }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyToOutOld(UnplacedTube const &tube,
                             Vector3D<typename Backend::precision_v> const &point,
                             typename Backend::precision_v &safety)
  {
     using namespace TubeTypes;
     using namespace TubeUtilities;
     typedef typename Backend::precision_v Float_t;

     safety = tube.z() - Abs(point.z());
     Float_t r = Sqrt(point.x()*point.x() + point.y()*point.y());
     Float_t safermax = tube.rmax() - r;
     MaskedAssign( safermax<safety, safermax, &safety );

     if(checkRminTreatment<tubeTypeT>(tube)) {
       Float_t safermin = r - tube.rmin();
       MaskedAssign( safermin<safety, safermin, &safety );
     }

    if(checkPhiTreatment<tubeTypeT>(tube)) {
      Float_t safephi;
      PhiPlaneSafety<Backend, tubeTypeT, true>(tube, point, safephi);
      MaskedAssign( safephi<safety, safephi, &safety);
    }
  }

}; // End of struct TubeImplementation
} // end of inline NS
} // End global namespace

#endif // VECGEOM_VOLUMES_KERNEL_TUBEIMPLEMENTATION_H_
