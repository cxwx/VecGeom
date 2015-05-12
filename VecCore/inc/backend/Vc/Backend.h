/// \file vc/backend.h

#ifndef VECCORE_BACKEND_VCBACKEND_H_
#define VECCORE_BACKEND_VCBACKEND_H_

//#include "VecCoreGlobal.h"
//#include "backend/scalar/Backend.h"

#include <Vc/Vc>

namespace VecCore {

namespace Backend {

namespace Vector {

//inline namespace VECGEOM_IMPL_NAMESPACE {

#define VECGEOM_INLINE
template< typename RealInput_t = double >
struct kVc {
  // The convention for type names is
  // _t is a type but it should be a primitive or otherwise scalar type
  // _v means that this type could potentially be a vector type ( but not necessarily )

  // define the base primitive type
  typedef RealInput_t                 Real_t;
  // we could define the other base primitive types
  typedef double                      Double_t;
  typedef double                      Float_t;

  // define vector types for wanted precision
  typedef typename Vc::Vector<Real_t>          Real_v;
  typedef typename Vc::Vector<Real_t>::Mask    Bool_v; // ( should we call this RealBool )

  // define vector types for standard double + float types
  typedef typename Vc::Vector<double>          Double_v;
  typedef typename Vc::Vector<float>           Float_v;
  typedef typename Vc::Vector<double>::Mask    DoubleBool_v;
  typedef typename Vc::Vector<float>::Mask     FloatBool_v;

  typedef typename Vc::int_v                   Int_v;
  typedef Int_v                       Inside_v;
  typedef typename Vc::int_v::Mask             IntBool_v;

  // alternative typedefs ( might supercede above typedefs )
  // We can no longer define Double_t like this this !!
  //typedef Vc::int_v                   Int_t;
  //typedef Vc::Vector<Precision>       Double_t;
  //typedef Vc::Vector<Precision>::Mask Bool_t;
  typedef Vc::Vector<Real_t>          Index_t;

  // numeric constants of this backend
  // can they be made constexpr ?
  const static Real_v kOne;
  const static Real_v kZero;
  const static Bool_v kTrue;
  const static Bool_v kFalse;

  // what about numeric constants for Double_v or Float_v
  // could by like
  static Double_v kDoubleOne;
  static Double_v kDoubleZero;
  static DoubleBool_v kDoubleTrue;
  static DoubleBool_v kDoubleFalse;

  // other properties of this backend
  // this is to be benchmarked ... could also make it a
  // configurable template parameter of this trait
  constexpr static bool early_returns = false;

  // the VectorSizes
  constexpr static int kRealVectorSize = Real_v::Size;
  constexpr static int kDoubleVectorSize = Double_v::Size;
  constexpr static int kFloatVectorSize = Float_v::Size;
  // ... in principle also kIntVectorSize ...
};

// Functions to extract individual components of vectors
// meant to provide a homogeneous way to provide operator[]
// to both complex and primitive types
// (since operator[] does not exist for primitive double, float, ... )
// function gives back component index of x and this component
// x can be modified

// Alternative Names: Operator[]
template <typename Type>
VECGEOM_INLINE
static
void
ComponentAssign( int index, typename Vc::Vector<Type>::EntryType what, typename Vc::Vector<Type> & to ) {
    to[index]=what;
}


template <typename Type>
VECGEOM_INLINE
static
typename Vc::Vector<Type>::EntryType const
GetComponent( typename Vc::Vector<Type> const & x, int index ) {
    return x[index];
}


// same for Mask
template <typename Type>
VECGEOM_INLINE
static
bool &
GetWritableComponent( typename Vc::Vector<Type>::Mask & x, int index ) {
    return x(index);
}

// might need to abstract on other things
// LoadFromArray; StoreToArray


template <typename Type>
VECGEOM_INLINE
static
void CondAssign(typename Vc::Vector<Type>::Mask const &cond,
                Vc::Vector<Type> const &thenval,
                Vc::Vector<Type> const &elseval,
                Vc::Vector<Type> *const output) {
  (*output)(cond) = thenval;
  (*output)(!cond) = elseval;
}

template <typename Type>
VECGEOM_INLINE
void CondAssign(typename Vc::Vector<Type>::Mask const &cond,
                Type const &thenval,
                Type const &elseval,
                Vc::Vector<Type> *const output) {
  (*output)(cond) = thenval;
  (*output)(!cond) = elseval;
}

template <typename Type>
VECGEOM_INLINE
void MaskedAssign(typename Vc::Vector<Type>::Mask const &cond,
                  Vc::Vector<Type> const &thenval,
                  Vc::Vector<Type> *const output) {
  (*output)(cond) = thenval;
}

template <typename Type>
VECGEOM_INLINE
void MaskedAssign(typename Vc::Vector<Type>::Mask const &cond,
                  Type const &thenval,
                  Vc::Vector<Type> *const output) {
  (*output)(cond) = thenval;
}

// special version of MaskedAssignment when output
// is of type int_v
//VECGEOM_INLINE
//void MaskedAssign(VcBool const &cond,
//                  const Inside_t thenval,
//                  VcInside *const output) {
//  (*output)(VcInside::Mask(cond)) = thenval;
//}

// returns if all lanes/slots in (vector) condition are true
template <typename Type>
VECGEOM_INLINE
bool IsFull(typename Vc::Vector<Type>::Mask const &cond) {
  return cond.isFull();
}

// returns if any lane/slot in (vector) condition is true
template <typename Type>
VECGEOM_INLINE
bool Any(typename Vc::Vector<Type>::Mask const &cond) {
  return !cond.isEmpty();
}

// returns if all lanes/slots in (vector) condition are false
template <typename Type>
VECGEOM_INLINE
bool IsEmpty(typename Vc::Vector<Type>::Mask const &cond) {
  return cond.isEmpty();
}

template <typename Type>
VECGEOM_INLINE
typename Vc::Vector<Type> Abs(typename Vc::Vector<Type> const &val) {
  return Vc::abs(val);
}

template <typename Type>
VECGEOM_INLINE
typename Vc::Vector<Type> Sqrt(typename Vc::Vector<Type> const &val) {
  return Vc::sqrt(val);
}

template <typename Type>
VECGEOM_INLINE
typename Vc::Vector<Type> ATan2(typename Vc::Vector<Type> const &y,
                  typename Vc::Vector<Type> const &x) {
  return Vc::atan2(y, x);
}

template <typename Type>
VECGEOM_INLINE
typename Vc::Vector<Type> sin(typename Vc::Vector<Type> const &x) {
  return Vc::sin(x);
}

template <typename Type>
VECGEOM_INLINE
typename Vc::Vector<Type> cos(typename Vc::Vector<Type> const &x) {
  return Vc::cos(x);
}

template <typename Type>
VECGEOM_INLINE
typename Vc::Vector<Type> tan(typename Vc::Vector<Type> const &radians) {
  // apparently Vc does not have a tan function
  // return Vc::tan(radians);
  // emulating it for the moment
  typename Vc::Vector<Type> s,c;
  Vc::sincos(radians,&s,&c);
  return s/c;
}

// ??????
template <typename Type>
VECGEOM_INLINE
typename Vc::Vector<Type> Pow(typename Vc::Vector<Type> const &x,
                              typename Vc::Vector<Type> & arg) {
    // What about a Vc version ?
    return std::pow(x,arg);
}

template <typename Type>
VECGEOM_INLINE
typename Vc::Vector<Type> Min(typename Vc::Vector<Type> const &val1,
                              typename Vc::Vector<Type> const &val2) {
  return Vc::min(val1, val2);
}

template <typename Type>
VECGEOM_INLINE
typename Vc::Vector<Type> Max(typename Vc::Vector<Type> const &val1,
                              typename Vc::Vector<Type> const &val2) {
  return Vc::max(val1, val2);
}


template <typename Type>
VECGEOM_INLINE
typename Vc::Vector<Type> Floor( typename Vc::Vector<Type> const &val ){
  return Vc::floor( val );
}

//} // End inline namespace

} // end Vector namespace

} // end Backend namespace

} // End global namespace


#endif // VECCORE_BACKEND_VCBACKEND_H_
