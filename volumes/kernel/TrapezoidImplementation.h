/*
 * @file TrapezoidImplementation.h
 * @author Guilherme Lima (lima at fnal dot gov)
 *
 * 2014-05-14 G.Lima - Created using interface from Johannes's parallelepiped example and UTrap's implementation code
 */

#ifndef VECGEOM_VOLUMES_KERNEL_TRAPEZOIDIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_TRAPEZOIDIMPLEMENTATION_H_

#include "base/Global.h"

#include "base/Transformation3D.h"
#include "volumes/kernel/GenericKernels.h"
#include "volumes/UnplacedTrapezoid.h"

#include <cstdio>

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v(TrapezoidImplementation, TranslationCode, translation::kGeneric, RotationCode, rotation::kGeneric)

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedTrapezoid;
using TrapSidePlane = vecgeom::UnplacedTrapezoid::TrapSidePlane;

template <TranslationCode transCodeT, RotationCode rotCodeT>
struct TrapezoidImplementation {

  static const int transC = transCodeT;
  static const int rotC   = rotCodeT;

  using PlacedShape_t = PlacedTrapezoid;
  using UnplacedShape_t = UnplacedTrapezoid;

  VECGEOM_CUDA_HEADER_BOTH
  static void PrintType() {
     printf("SpecializedTrapezoid<%i, %i>", transCodeT, rotCodeT);
  }

  template <typename Stream>
  static void PrintType( Stream & s ) {
    s << "SpecializedTrapezoid<" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintImplementationType(Stream &s) {
    s << "TrapezoidImplemenation<" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintUnplacedType(Stream &s) { s << "UnplacedTrapezoid"; }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void UnplacedContains(
      UnplacedTrapezoid const &unplaced,
      Vector3D<typename Backend::precision_v> const &point,
      typename Backend::bool_v &inside);

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void Contains(
      UnplacedTrapezoid const &unplaced,
      Transformation3D const &transformation,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> &localPoint,
      typename Backend::bool_v &inside);

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void Inside(UnplacedTrapezoid const &unplaced,
      Transformation3D const &transformation,
      Vector3D<typename Backend::precision_v> const &masterPoint,
      typename Backend::inside_v &inside);

  template <typename Backend, bool ForInside>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void GenericKernelForContainsAndInside(
      UnplacedTrapezoid const &unplaced,
      Vector3D<typename Backend::precision_v> const &,
      typename Backend::bool_v &,
      typename Backend::bool_v &);

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void DistanceToIn(
      UnplacedTrapezoid const &unplaced,
      Transformation3D const &transformation,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> const &direction,
      typename Backend::precision_v const &stepMax,
      typename Backend::precision_v &distance);

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void DistanceToOut(
      UnplacedTrapezoid const &unplaced,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> const &direction,
      typename Backend::precision_v const &stepMax,
      typename Backend::precision_v &distance);

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyToIn(UnplacedTrapezoid const &unplaced,
                         Transformation3D const &transformation,
                         Vector3D<typename Backend::precision_v> const &point,
                         typename Backend::precision_v &safety);

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyToOut(UnplacedTrapezoid const &unplaced,
                          Vector3D<typename Backend::precision_v> const &point,
                          typename Backend::precision_v &safety);

};


template <TranslationCode transCodeT, RotationCode rotCodeT>
template <typename Backend, bool ForInside>
VECGEOM_CUDA_HEADER_BOTH
void TrapezoidImplementation<transCodeT, rotCodeT>::GenericKernelForContainsAndInside(
    UnplacedTrapezoid const &unplaced,
    Vector3D<typename Backend::precision_v> const &localPoint,
    typename Backend::bool_v &completelyInside,
    typename Backend::bool_v &completelyOutside) {

  // z-region
  completelyOutside = Abs(localPoint[2]) > MakePlusTolerant<ForInside>( unplaced.GetDz() );
  if ( Backend::early_returns && IsFull(completelyOutside) )  return;
  if (ForInside) {
    completelyInside = Abs(localPoint[2]) < MakeMinusTolerant<ForInside>( unplaced.GetDz() );
  }


#ifndef VECGEOM_PLANESHELL_DISABLE
  unplaced.GetPlanes() -> GenericKernelForContainsAndInside<Backend,ForInside>(
      localPoint, completelyInside, completelyOutside );
#else
  typedef typename Backend::precision_v Float_t;

  TrapSidePlane const* fPlanes = unplaced.GetPlanes();
  Float_t dist[4];
  for(unsigned int i=0; i<4; ++i) {
    dist[i] = fPlanes[i].fA*localPoint.x() + fPlanes[i].fB*localPoint.y()
            + fPlanes[i].fC*localPoint.z() + fPlanes[i].fD;
  }

  for(unsigned int i=0; i<4; ++i) {
    // is it outside of this side plane?
    completelyOutside |= dist[i] > MakePlusTolerant<ForInside>(0.);
    if ( Backend::early_returns && IsFull(completelyOutside) )  return;
    if ( ForInside ) {
      completelyInside &= dist[i] < MakeMinusTolerant<ForInside>(0.);
    }
  }
#endif

  return;
}

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void TrapezoidImplementation<transCodeT, rotCodeT>::UnplacedContains(
  UnplacedTrapezoid const &unplaced,
  Vector3D<typename Backend::precision_v> const &point,
  typename Backend::bool_v &inside) {

  typedef typename Backend::bool_v Bool_t;
  Bool_t unused, outside;
  GenericKernelForContainsAndInside<Backend, false>(unplaced, point, unused, outside);
  inside = !outside;
}

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void TrapezoidImplementation<transCodeT, rotCodeT>::Contains(
    UnplacedTrapezoid const &unplaced,
    Transformation3D const &transformation,
    Vector3D<typename Backend::precision_v> const &masterPoint,
    Vector3D<typename Backend::precision_v> &localPoint,
    typename Backend::bool_v &inside) {

  localPoint = transformation.Transform<transCodeT, rotCodeT>(masterPoint);

  typedef typename Backend::bool_v Bool_t;
  Bool_t unused, outside;
  GenericKernelForContainsAndInside<Backend,false>(unplaced, localPoint,
                                                   unused, outside);

  inside = !outside;
}

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void TrapezoidImplementation<transCodeT, rotCodeT>::Inside(
    UnplacedTrapezoid const &unplaced,
    Transformation3D const &transformation,
    Vector3D<typename Backend::precision_v> const &masterPoint,
    typename Backend::inside_v &inside) {

  typedef typename Backend::bool_v Bool_t;

  // convert from master to local coordinates
  Vector3D<typename Backend::precision_v> localPoint =
      transformation.Transform<transCodeT, rotCodeT>(masterPoint);

  Bool_t fullyinside, fullyoutside;
  GenericKernelForContainsAndInside<Backend,true>(
      unplaced, localPoint, fullyinside, fullyoutside);

  inside=EInside::kSurface;
  MaskedAssign(fullyinside,  EInside::kInside,  &inside);
  MaskedAssign(fullyoutside, EInside::kOutside, &inside);
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside - return kInfinity if no
// intersection.
//
// ALGORITHM: For each component (z-planes, side planes), calculate
// pair of minimum (smin) and maximum (smax) intersection values for
// which the particle is in the extent of the shape.  The point of
// entrance (exit) is found by the largest smin (smallest smax).
//
//  If largest smin > smallest smax, the trajectory does not reach
//  inside the shape.
//
template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void TrapezoidImplementation<transCodeT, rotCodeT>::DistanceToIn(
    UnplacedTrapezoid const &unplaced,
    Transformation3D const &transformation,
    Vector3D<typename Backend::precision_v> const &masterPoint,
    Vector3D<typename Backend::precision_v> const &masterDir,
    typename Backend::precision_v const &/*stepMax*/,
    typename Backend::precision_v &distance) {

  typedef typename Backend::precision_v Float_t;
  typedef typename Backend::bool_v Bool_t;

  Vector3D<Float_t> point = transformation.Transform<transCodeT, rotCodeT>(masterPoint);
  Vector3D<Float_t> dir   = transformation.TransformDirection<rotCodeT>(masterDir);

  distance = kInfinity;    // set to early returned value
  Bool_t done = Backend::kFalse;

  //
  // Step 1: find range of distances along dir between Z-planes (smin, smax)
  //

  // convenience variables for direction pointing to +z or -z.
  // Note that both posZdir and NegZdir may be false, if dir.z() is zero!
  Bool_t posZdir  = dir.z() > 0.0;
  Bool_t negZdir  = dir.z() < 0.0;
  Float_t zdirSign = Backend::kOne;  // z-direction
  MaskedAssign( negZdir, -Backend::kOne, &zdirSign);

  Float_t max = zdirSign*unplaced.GetDz() - point.z();     // z-dist to farthest z-plane

  // step 1.a) input particle is moving away --> return infinity

  // check if moving away towards +z
  done |= ( posZdir && max < MakePlusTolerant<true>(0.));

  // check if moving away towards -z
  done |= ( negZdir && max > MakeMinusTolerant<true>(0.));

  // if all particles moving away, we're done
  if (Backend::early_returns && IsFull(done)) return;

  // Step 1.b) General case:
  //   smax,smin are range of distances within z-range, taking direction into account.
  //   smin<smax - smax is positive, but smin may be either positive or negative
  Float_t dirFactor = Backend::kOne / dir.z();     // convert distances from z to dir
  Float_t smax = max * dirFactor;
  Float_t smin = (-zdirSign*unplaced.GetDz() - point.z())*dirFactor;

  // Step 1.c) special case: if dir is perpendicular to z-axis...
  Bool_t test = (!posZdir) && (!negZdir);

  // ... and out of z-range, then trajectory will not intercept volume
  Bool_t zrange = Abs(point.z()) < MakeMinusTolerant<true>(unplaced.GetDz());
  done |= ( test && !zrange );
  if (Backend::early_returns && IsFull(done)) return;

  // ... or within z-range, then smin=-1, smax=infinity for now
  MaskedAssign( test && zrange,      -1.0, &smin );
  MaskedAssign( test && zrange, kInfinity, &smax );

  //
  // Step 2: find distances for intersections with side planes.
  //

#ifndef VECGEOM_PLANESHELL_DISABLE
  // If disttoplanes is such that smin < dist < smax, then distance=disttoplanes
  Float_t disttoplanes = unplaced.GetPlanes()->DistanceToIn<Backend>(point, smin, dir);

  //=== special cases
  // 1) point is inside (wrong-side) --> return -1
  Bool_t inside = (smin<-kHalfTolerance && disttoplanes<-kHalfTolerance);
  MaskedAssign( !done && inside, -1.0, &distance );  // wrong-side point
  done |= inside;

  // 2) point is outside plane-shell and moving away --> return infinity
  done |= disttoplanes == kInfinity;

  // 3) track misses the trapezoid: smin<0 && disttoplanes>smax
  done |= (smin < 0 && disttoplanes > smax);
  if (Backend::early_returns && IsFull(done)) return;

  // at this point we know there is a valid distance - start with the z-plane based one
  MaskedAssign( !done, smin, &distance);

  // and then take the one from the planar shell, if valid
  Bool_t hitplanarshell = (disttoplanes > smin) && (disttoplanes < smax);
  MaskedAssign( !done && hitplanarshell, disttoplanes, &distance);

  // finally check whether entry point candidate is globally valid
  // --> if it is, entry point should be on the surface.
  //     if entry point candidate is completely outside, track misses the solid --> return infinity
  Vector3D<Float_t> entry = point + distance*dir;
  Bool_t complIn, complOut;
  GenericKernelForContainsAndInside<Backend,true>( unplaced, entry, complIn, complOut );
  MaskedAssign( !done && complOut, kInfinity, &distance );

#else
  // here for VECGEOM_PLANESHELL_DISABLE
  TrapSidePlane const* fPlanes = unplaced.GetPlanes();

  // loop over side planes - find pdist,Comp for each side plane
  Float_t pdist[4], comp[4], vdist[4];
  // auto-vectorizable part of loop
  for (unsigned int i = 0; i < 4; ++i) {
    // Note: normal vector is pointing outside the volume (convention), therefore
    // pdist>0 if point is outside  and  pdist<0 means inside
    pdist[i] = fPlanes[i].fA * point.x() + fPlanes[i].fB * point.y()
      + fPlanes[i].fC * point.z() + fPlanes[i].fD;

    // Comp is projection of dir over the normal vector of side plane, hence
    // Comp > 0 if pointing ~same direction as normal and Comp<0 if ~opposite to normal
    comp[i] = fPlanes[i].fA * dir.x() + fPlanes[i].fB * dir.y() + fPlanes[i].fC * dir.z();

    vdist[i] = -pdist[i] / comp[i];
  }

  // wrong-side check: if (inside && smin<0) return -1
  Bool_t inside = smin < MakeMinusTolerant<true>(0.0)
               && smax > MakePlusTolerant<true>(0.0);

  for(int i=0; i<4; ++i) {
    inside &= pdist[i] < MakeMinusTolerant<true>(0.0);
  }
  MaskedAssign( inside, -1.0, &distance );
  done |= inside;
  if(Backend::early_returns && IsFull(done)) return;

  // this part does not auto-vectorize
  for(unsigned int i=0; i<4; ++i) {
    // if outside and moving away, return infinity
    Bool_t posPoint = pdist[i] > MakeMinusTolerant<true>(0.);
    Bool_t posDir = comp[i] > 0;

    // discard the ones moving away from this plane
    done |= (posPoint && posDir);

    // check if trajectory will intercept plane within current range (smin,smax)

    Bool_t interceptFromInside = (!posPoint && posDir);
    done |= ( interceptFromInside  && vdist[i]<smin );

    Bool_t interceptFromOutside = (posPoint && !posDir);
    done |= ( interceptFromOutside && vdist[i]>smax );
    if ( Backend::early_returns && IsFull(done) ) return;

    //   If dist is such that smin < dist < smax, then adjust either smin or smax.
    Bool_t validVdist = (vdist[i]>smin && vdist[i]<smax);
    MaskedAssign( !done && interceptFromInside  && validVdist, vdist[i], &smax );
    MaskedAssign( !done && interceptFromOutside && validVdist, vdist[i], &smin );
  }

  // Checks in non z plane intersections ensure smin<smax
  MaskedAssign( !done , smin, &distance );
#endif // end of #ifdef PLANESHELL
}


template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void TrapezoidImplementation<transCodeT, rotCodeT>::DistanceToOut(
    UnplacedTrapezoid const &unplaced,
    Vector3D<typename Backend::precision_v> const &point,
    Vector3D<typename Backend::precision_v> const &dir,
    typename Backend::precision_v const &/*stepMax*/,
    typename Backend::precision_v &distance) {

  typedef typename Backend::precision_v Float_t;
  typedef typename Backend::bool_v Bool_t;

  distance = kInfinity;  // default initialization
  Bool_t done(Backend::kFalse);

  // step 0: if point is outside any plane --> return -1
  Bool_t outside = Abs(point.z()) > MakePlusTolerant<true>(unplaced.GetDz());
  MaskedAssign( outside, -1.0, &distance );
  done |= outside;
  if(Backend::early_returns && IsFull(done)) return;
 
  //
  // Step 1: find range of distances along dir between Z-planes (smin, smax)
  //

  // convenience variables for direction pointing to +z or -z.
  // Note that both posZdir and NegZdir may be false, if dir.z() is zero!
  //Bool_t posZdir = dir.z() > 0.0;
  Bool_t negZdir = dir.z() <= -0.0;  // -0.0 is needed, see JIRA-150

  // TODO: consider use of copysign or some other standard function
  Float_t zdirSign = Backend::kOne;  // z-direction
  MaskedAssign( negZdir, -Backend::kOne, &zdirSign);

  Float_t max = zdirSign*unplaced.GetDz() - point.z();  // z-dist to farthest z-plane

  // // do we need this ????
  // // check if moving away towards +z
  // done |= (posZdir && max <= MakePlusTolerant<true>(0.));
  // // check if moving away towards -z
  // done |= (negZdir && max >= MakeMinusTolerant<true>(0.));

  // // if all particles moving away, we're done
  // MaskedAssign( done, Float_t(0.0), &distance );
  // if ( IsFull(done) ) return;

  // Step 1.b) general case: assign distance to z plane
  distance = max/( dir.z() + zdirSign * vecgeom::kEpsilon );

  //
  // Step 2: find distances for intersections with side planes.
  //

#ifndef VECGEOM_PLANESHELL_DISABLE
  Float_t disttoplanes = unplaced.GetPlanes()->DistanceToOut<Backend>(point, dir);

  // reconcile sides with endcaps
  MaskedAssign( disttoplanes<distance, disttoplanes, &distance );

#else
  //=== Here for VECGEOM_PLANESHELL_DISABLE
  TrapSidePlane const* fPlanes = unplaced.GetPlanes();

  // loop over side planes - find pdist,Comp for each side plane
  Float_t pdist[4], comp[4], vdist[4];
  for (unsigned int i = 0; i < 4; ++i) {
    // Note: normal vector is pointing outside the volume (convention), therefore
    // pdist>0 if point is outside  and  pdist<0 means inside
    pdist[i] = fPlanes[i].fA * point.x() + fPlanes[i].fB * point.y()
      + fPlanes[i].fC * point.z() + fPlanes[i].fD;

    // Comp is projection of dir over the normal vector of side plane, hence
    // Comp > 0 if pointing ~same direction as normal and Comp<0 if pointing ~opposite to normal
    comp[i] = fPlanes[i].fA * dir.x() + fPlanes[i].fB * dir.y() + fPlanes[i].fC * dir.z();

    vdist[i] = -pdist[i] / comp[i];
  }

  // early return if point is outside of plane
  for (unsigned int i=0; i<4; ++i) {
    done |= ( pdist[i] > MakePlusTolerant<true>(0.));
  }
  MaskedAssign( done, -1.0, &distance );
  if (Backend::early_returns && IsFull(done)) return;

  for (unsigned int i = 0; i < 4; ++i) {
    // if point is inside, pointing towards plane and vdist<distance, then distance=vdist
    Bool_t test = (vdist[i] >= MakeMinusTolerant<true>(0.)) && comp[i]>0.0 && vdist[i]<distance;
    MaskedAssign(!done && test, vdist[i], &distance);
  }
#endif
}


template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void TrapezoidImplementation<transCodeT, rotCodeT>::SafetyToIn(
  UnplacedTrapezoid const &unplaced,
    Transformation3D const &transformation,
    Vector3D<typename Backend::precision_v> const &masterPoint,
    typename Backend::precision_v &safety) {

  typedef typename Backend::precision_v Float_t;

  // convert from master to local coordinates
  Vector3D<Float_t> point = transformation.Transform<transCodeT, rotCodeT>(masterPoint);
  safety = Abs(point.z()) - unplaced.GetDz();

#ifndef VECGEOM_PLANESHELL_DISABLE
  // Get safety over side planes
  unplaced.GetPlanes()->SafetyToIn<Backend>(point, safety);
#else
  // Loop over side planes
  TrapSidePlane const* fPlanes = unplaced.GetPlanes();
  Float_t Dist[4];
  for (int i = 0; i < 4; ++i) {
    Dist[i] = fPlanes[i].fA * point.x() + fPlanes[i].fB * point.y()
      + fPlanes[i].fC * point.z() + fPlanes[i].fD;
  }

  for (int i = 0; i < 4; ++i) {
    MaskedAssign( Dist[i]>safety, Dist[i], &safety );
  }

  MaskedAssign( safety < 0, 0.0, &safety );
#endif
}

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void TrapezoidImplementation<transCodeT, rotCodeT>::SafetyToOut(
    UnplacedTrapezoid const &unplaced,
    Vector3D<typename Backend::precision_v> const &point,
    typename Backend::precision_v &safety) {

  typedef typename Backend::bool_v Bool_t;

  Bool_t done(false);

  // If point is outside (wrong-side) --> safety to negative value
  safety = unplaced.GetDz() - Abs(point.z());
  done |= (safety < kHalfTolerance);

  // If all test points are outside, we're done
  if ( Backend::early_returns && IsFull(done) ) return;

#ifndef VECGEOM_PLANESHELL_DISABLE
  // Get safety over side planes
  unplaced.GetPlanes()->SafetyToOut<Backend>(point, safety);
#else
  // Loop over side planes
  TrapSidePlane const* fPlanes = unplaced.GetPlanes();

  // auto-vectorizable loop
  typename Backend::precision_v Dist[4];
  for (int i = 0; i < 4; ++i) {
    Dist[i] = -(fPlanes[i].fA * point.x() + fPlanes[i].fB * point.y()
		+ fPlanes[i].fC * point.z() + fPlanes[i].fD);
  }

  // unvectorizable loop
  for (int i = 0; i < 4; ++i) {
    MaskedAssign( !done && Dist[i] < safety, Dist[i], &safety );
  }
#endif
}

} // End of inline VECGEOM_IMPL_NAMESPACE namespace
} // End global namespace

#endif // VECGEOM_VOLUMES_KERNEL_TRAPEZOIDIMPLEMENTATION_H_
