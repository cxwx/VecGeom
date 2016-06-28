#ifndef VOLUMES_UNPLACEDVOLUMEIMPLHELPER_H_
#define VOLUMES_UNPLACEDVOLUMEIMPLHELPER_H_

#include "base/Global.h"
#include "base/SOA3D.h"
#include "volumes/UnplacedVolume.h"
#include "management/VolumeFactory.h"

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2t(class, CommonUnplacedVolumeImplHelper, typename, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {

// kernels
template <class Implementation, typename Real_v>
VECGEOM_INLINE
static void DistanceToOutLoop(typename Implementation::UnplacedStruct_t const *shapestruct, const size_t offset,
                              const size_t size, SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                              double const *step_max, double *output)
{
  using vecCore::FromPtr;
  for (decltype(points.size()) i(offset); i < size; i += vecCore::VectorSize<Real_v>()) {
    Vector3D<Real_v> point(FromPtr<Real_v>(points.x() + i), FromPtr<Real_v>(points.y() + i),
                           FromPtr<Real_v>(points.z() + i));
    Vector3D<Real_v> dir(FromPtr<Real_v>(directions.x() + i), FromPtr<Real_v>(directions.y() + i),
                         FromPtr<Real_v>(directions.z() + i));
    Real_v stepMax_v = FromPtr<Real_v>(&step_max[i]);
    Real_v result;
    Implementation::template DistanceToOut<Real_v>(*shapestruct, point, dir, stepMax_v, result);
    vecCore::Store(result, &output[i]);
  }
}

template <class Implementation, typename Real_v>
VECGEOM_INLINE
static void SafetyToOutLoop(typename Implementation::UnplacedStruct_t const *shapestruct, const size_t offset,
                            const size_t size, SOA3D<Precision> const &points, double *output)
{
  using vecCore::FromPtr;
  const decltype(points.size()) len(size);
  for (decltype(points.size()) i(offset); i < len; i += vecCore::VectorSize<Real_v>()) {
    Vector3D<Real_v> point(FromPtr<Real_v>(points.x() + i), FromPtr<Real_v>(points.y() + i),
                           FromPtr<Real_v>(points.z() + i));
    Real_v result(kInfinity);
    Implementation::template SafetyToOut<Real_v>(*shapestruct, point, result);
    vecCore::Store(result, &output[i]);
  }
}

template <class Implementation, class BaseUnplVol = VUnplacedVolume>
class CommonUnplacedVolumeImplHelper : public BaseUnplVol {

public:
  using UnplacedStruct_t = typename Implementation::UnplacedStruct_t;
  using UnplacedVolume_t = typename Implementation::UnplacedVolume_t;

  // bring in constructors
  using BaseUnplVol::BaseUnplVol;
  // bring in some members from base (to avoid name hiding)
  using BaseUnplVol::SafetyToOut;
  using BaseUnplVol::SafetyToIn;
  using BaseUnplVol::DistanceToOut;
  using BaseUnplVol::DistanceToIn;

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  virtual Precision DistanceToOut(Vector3D<Precision> const &p, Vector3D<Precision> const &d,
                                  Precision step_max = kInfinity) const override
  {
#ifndef VECGEOM_NVCC
    assert(d.IsNormalized() && " direction not normalized in call to  DistanceToOut ");
#endif
    Precision output;
    Implementation::template DistanceToOut(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), p, d, step_max,
                                           output);

//#ifdef VECGEOM_DISTANCE_DEBUG
//    DistanceComparator::CompareDistanceToOut( this, output, point, direction, stepMax );
//#endif

// detect -inf responses which are often an indication for a real bug
#ifndef VECGEOM_NVCC
    assert(!((output < 0.) && std::isinf((double)output)));
#endif
    return output;
  }

  // the extended DistanceToOut interface
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  virtual Precision DistanceToOut(Vector3D<Precision> const &p, Vector3D<Precision> const &d,
                                  Vector3D<Precision> &normal, bool &convex,
                                  Precision step_max = kInfinity) const override
  {
    (void)p;
    (void)d;
    (void)normal;
    (void)convex;
    (void)step_max;
    assert(false);
    return 0.;
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  virtual bool Contains(Vector3D<Precision> const &p) const override
  {
    bool output(false);
    Implementation::Contains(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), p, output);
    return output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  virtual EnumInside Inside(Vector3D<Precision> const &p) const override
  {
    using Inside_t = vecCore::Int32_s;
    Inside_t output(0);
    Implementation::Inside(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), p, output);
    return (EnumInside)output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision DistanceToIn(Vector3D<Precision> const &p, Vector3D<Precision> const &d,
                                 const Precision step_max = kInfinity) const override
  {
#ifndef VECGEOM_NVCC
    assert(d.IsNormalized() && " direction not normalized in call to  DistanceToOut ");
#endif
    Precision output(kInfinity);
    Implementation::DistanceToIn(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), p, d, step_max, output);
    return output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision SafetyToOut(Vector3D<Precision> const &p) const override
  {
    Precision output(kInfinity);
    Implementation::SafetyToOut(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), p, output);
    return output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision SafetyToIn(Vector3D<Precision> const &p) const override
  {
    Precision output(kInfinity);
    Implementation::SafetyToIn(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), p, output);
    return output;
  }

  virtual int memory_size() const override { return sizeof(*this); }
};

// helper that dispatches the vector interface to a loop over the SIMD interface
template <class Implementation, class BaseUnplVol = VUnplacedVolume>
class SIMDUnplacedVolumeImplHelper : public CommonUnplacedVolumeImplHelper<Implementation, BaseUnplVol> {
public:
  using Real_v           = vecgeom::VectorBackend::Real_v;
  using UnplacedVolume_t = typename Implementation::UnplacedVolume_t;
  using Common_t         = CommonUnplacedVolumeImplHelper<Implementation, BaseUnplVol>;

  // bring in constructor
  using Common_t::Common_t;

  // the explicit SIMD interface
  VECGEOM_CUDA_HEADER_BOTH
  virtual Real_v DistanceToOutVec(Vector3D<Real_v> const &p, Vector3D<Real_v> const &d,
                                  Real_v const &step_max) const override
  {
    Real_v output;
    Implementation::template DistanceToOut<Real_v>(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), p, d,
                                                   step_max, output);
    return output;
  }

  // the explicit SIMD interface
  VECGEOM_CUDA_HEADER_BOTH
  virtual Real_v DistanceToInVec(Vector3D<Real_v> const &p, Vector3D<Real_v> const &d,
                                 Real_v const &step_max) const override
  {
    Real_v output(kInfinity);
    Implementation::template DistanceToIn<Real_v>(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), p, d,
                                                  step_max, output);
    return output;
  }

  // the explicit SIMD interface
  VECGEOM_CUDA_HEADER_BOTH
  virtual Real_v SafetyToOutVec(Vector3D<Real_v> const &p) const override
  {
    Real_v output(kInfinity);
    Implementation::template SafetyToOut<Real_v>(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), p, output);
    return output;
  }

  // the explicit SIMD interface
  VECGEOM_CUDA_HEADER_BOTH
  virtual Real_v SafetyToInVec(Vector3D<Real_v> const &p) const override
  {
    Real_v output(kInfinity);
    Implementation::template SafetyToIn<Real_v>(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), p, output);
    return output;
  }

  using UnplacedStruct_t = typename Implementation::UnplacedStruct_t;
  using Common_t::DistanceToOut;
  using Common_t::SafetyToOut;

  virtual void DistanceToOut(SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                             Precision const *const step_max, Precision *const output) const override
  {
    auto offset = points.size() - points.size() % vecCore::VectorSize<VectorBackend::Real_v>();
    auto &shape = ((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct();
    // vector loop treatment
    DistanceToOutLoop<Implementation, VectorBackend::Real_v>(&shape, 0, offset, points, directions, step_max, output);
    // tail treatment
    DistanceToOutLoop<Implementation, ScalarBackend::Real_v>(&shape, offset, points.size(), points, directions,
                                                             step_max, output);
  }

  virtual void SafetyToOut(SOA3D<Precision> const &points, Precision *const output) const override
  {
    auto offset = points.size() - points.size() % vecCore::VectorSize<VectorBackend::Real_v>();
    auto &shape = ((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct();
    // vector loop treatment
    SafetyToOutLoop<Implementation, VectorBackend::Real_v>(&shape, 0, offset, points, output);
    // tail treatment
    SafetyToOutLoop<Implementation, ScalarBackend::Real_v>(&shape, offset, points.size(), points, output);
  }
};

// helper that dispatches the vector interface to a loop over the scalar interface
template <class Implementation, class BaseUnplVol = VUnplacedVolume>
class LoopUnplacedVolumeImplHelper : public CommonUnplacedVolumeImplHelper<Implementation, BaseUnplVol> {
public:
  using Real_v           = vecgeom::VectorBackend::Real_v;
  using Real_s           = vecgeom::ScalarBackend::Real_v;
  using UnplacedVolume_t = typename Implementation::UnplacedVolume_t;
  using Common_t         = CommonUnplacedVolumeImplHelper<Implementation, BaseUnplVol>;

  // constructors
  using Common_t::Common_t;
  using Common_t::SafetyToOut;
  using Common_t::DistanceToOut;
  using Common_t::SafetyToIn;
  using Common_t::DistanceToIn;

  // the explicit SIMD interface
  VECGEOM_CUDA_HEADER_BOTH
  virtual Real_v DistanceToOutVec(Vector3D<Real_v> const &p, Vector3D<Real_v> const &d,
                                  Real_v const &step_max) const override
  {
    // implementation of a vector interface in terms of a scalar interface
    Real_v output(kInfinity);
    using vecCore::LaneAt;
    for (size_t i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
      Vector3D<Real_s> ps(LaneAt(p.x(), i), LaneAt(p.y(), i), LaneAt(p.z(), i)); // scalar vector
      Vector3D<Real_s> ds(LaneAt(d.x(), i), LaneAt(d.y(), i), LaneAt(d.z(), i)); // scalar direction;
      Real_s result;
      Implementation::template DistanceToOut<Real_s>(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), ps, ds,
                                                     LaneAt(step_max, i), result);
      vecCore::AssignLane(output, i, result);
    }
    return output;
  }

  // the explicit SIMD interface
  VECGEOM_CUDA_HEADER_BOTH
  virtual Real_v DistanceToInVec(Vector3D<Real_v> const &p, Vector3D<Real_v> const &d,
                                 Real_v const &step_max) const override
  {
    // implementation of a vector interface in terms of a scalar interface
    Real_v output(kInfinity);
    using vecCore::LaneAt;
    for (size_t i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
      Vector3D<Real_s> ps(LaneAt(p.x(), i), LaneAt(p.y(), i), LaneAt(p.z(), i)); // scalar vector
      Vector3D<Real_s> ds(LaneAt(d.x(), i), LaneAt(d.y(), i), LaneAt(d.z(), i)); // scalar direction;
      Real_s tmp(-1.);
      Implementation::template DistanceToIn<Real_s>(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), ps, ds,
                                                    LaneAt(step_max, i), tmp);
      vecCore::AssignLane(output, i, tmp);
    }
    return output;
  }

  // the explicit SIMD interface
  VECGEOM_CUDA_HEADER_BOTH
  virtual Real_v SafetyToOutVec(Vector3D<Real_v> const &p) const override
  {
    // implementation of a vector interface in terms of a scalar interface
    Real_v output(-1.);
    using vecCore::LaneAt;
    for (size_t i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
      Vector3D<Real_s> ps(LaneAt(p.x(), i), LaneAt(p.y(), i), LaneAt(p.z(), i)); // scalar vector
      Real_s tmp(kInfinity);
      Implementation::template SafetyToOut<Real_s>(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), ps, tmp);
      vecCore::AssignLane(output, i, tmp);
    }
    return output;
  }

  // the explicit SIMD interface
  VECGEOM_CUDA_HEADER_BOTH
  virtual Real_v SafetyToInVec(Vector3D<Real_v> const &p) const override
  {
    // implementation of a vector interface in terms of a scalar interface
    Real_v output(kInfinity);
    using vecCore::LaneAt;
    for (size_t i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
      Vector3D<Real_s> ps(LaneAt(p.x(), i), LaneAt(p.y(), i), LaneAt(p.z(), i)); // scalar vector
      Real_s tmp;
      Implementation::template SafetyToIn<Real_s>(((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct(), ps, tmp);
      vecCore::AssignLane(output, i, tmp);
    }
    return output;
  }

  using UnplacedStruct_t = typename Implementation::UnplacedStruct_t;

  virtual void DistanceToOut(SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                             Precision const *const step_max, Precision *const output) const override
  {
    auto &shape = ((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct();
    DistanceToOutLoop<Implementation, Precision>(&shape, 0, points.size(), points, directions, step_max, output);
  }

  virtual void SafetyToOut(SOA3D<Precision> const &points, Precision *const output) const override
  {
    auto &shape = ((UnplacedVolume_t *)this)->UnplacedVolume_t::GetStruct();
    SafetyToOutLoop<Implementation, Precision>(&shape, 0, points.size(), points, output);
  }
};
}
} // end namespace

#endif /* VOLUMES_UnplacedVolumeImplHelper_H_ */
