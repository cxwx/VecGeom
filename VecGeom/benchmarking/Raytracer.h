/// \file Raytracer.h
/// \author Andrei Gheata (andrei.gheata@cern.ch)

#ifndef VECGEOM_RAYTRACER_H_
#define VECGEOM_RAYTRACER_H_

#include <list>
#include <VecGeom/base/Global.h>
#include <VecGeom/base/Vector3D.h>
#include <VecGeom/benchmarking/Color.h>

#ifdef VECGEOM_ENABLE_CUDA
#include <VecGeom/backend/cuda/Interface.h>
#endif

namespace vecgeom {

enum ERTmodel { kRTxray = 0, kRTspecular, kRTtransparent, kRTdiffuse };
enum ERTView { kRTVparallel = 0, kRTVperspective };

VECGEOM_DEVICE_FORWARD_DECLARE(class VPlacedVolume;);

inline namespace VECGEOM_IMPL_NAMESPACE {

class NavigationState;
class VPlacedVolume;

struct Ray_t {
  using VPlacedVolumePtr_t = VPlacedVolume const *;

  Vector3D<double> fPos;
  Vector3D<double> fDir;
  NavigationState *fCrtState  = nullptr;
  NavigationState *fNextState = nullptr;
  VPlacedVolumePtr_t fVolume  = nullptr;
  int fNcrossed               = 0;
  Color_t fColor;
  bool fDone = false;
};

struct RaytracerData_t {

  using VPlacedVolumePtr_t = VPlacedVolume const *;

  double fScale     = 0;                 ///< Scaling from pixels to world coordinates
  double fShininess = 1.;                ///< Shininess exponent in the specular model
  double fZoom      = 1.;                ///< Zoom with respect to the default view
  Vector3D<double> fSourceDir;           ///< Light source direction
  Vector3D<double> fScreenPos;           ///< Screen position
  Vector3D<double> fStart;               ///< Eye position in perspectove mode
  Vector3D<double> fDir;                 ///< Start direction of all rays in parallel view mode
  Vector3D<double> fUp;                  ///< Up vector in the shooting rectangle plane
  Vector3D<double> fRight;               ///< Right vector in the shooting rectangle plane
  Vector3D<double> fLeftC;               ///< left-down corner of the ray shooting rectangle
  int fVerbosity      = 0;               ///< Verbosity level
  int fNrays          = 0;               ///< Number of rays left to propagate
  int fSize_px        = 1024;            ///< Image pixel size in x
  int fSize_py        = 1024;            ///< Image pixel size in y
  int fVisDepth       = 1;               ///< Visible geometry depth
  int fMaxDepth       = 0;               ///< Maximum geometry depth
  Color_t fLightColor = 0xFFFFFFFF;      ///< Light color
  Color_t fObjColor   = 0x0000FFFF;      ///< Object color
  ERTmodel fModel     = kRTxray;         ///< Selected RT model
  ERTView fView       = kRTVperspective; ///< View type

  VPlacedVolumePtr_t fWorld = nullptr; ///< World volume
  NavigationState *fVPstate = nullptr; ///< Navigation state corresponding to the viewpoint

  VECCORE_ATT_HOST_DEVICE
  void Print();
};

/// \brief Raytracing a logical volume content using a given model
///
/// In order to run a benchmark, a world volume must be provided to the
/// raytracer.
namespace Raytracer {

using VPlacedVolumePtr_t = VPlacedVolume const *;

/// \param world Mother volume containing daughters that will be benchmarked.
/// \param screen_position position of the screen in world coordinates, rays are starting from a plane passing through
/// it, normal vector pointing to the origin of the world reference frame \param up_vector the projection of this
/// vector on the camera plane determines the 'up' direction of the image \param img_size_px image size on X in pixels
/// \param img_size_py image size on Y in pixels
VECCORE_ATT_HOST_DEVICE
void InitializeModel(VPlacedVolumePtr_t world, RaytracerData_t &data);

VECCORE_ATT_HOST_DEVICE
void ApplyRTmodel(Ray_t &ray, double step, RaytracerData_t const &rtdata);

/// \brief Entry point to propagate all rays
VECCORE_ATT_HOST_DEVICE
void PropagateRays(RaytracerData_t &data);

VECCORE_ATT_HOST_DEVICE
Color_t RaytraceOne(int px, int py, RaytracerData_t const &rtdata);

// Navigation methods (just temporary here)
VECCORE_ATT_HOST_DEVICE
VPlacedVolumePtr_t LocateGlobalPoint(VPlacedVolume const *vol, Vector3D<Precision> const &point, NavigationState &path,
                                     bool top);
VECCORE_ATT_HOST_DEVICE
VPlacedVolumePtr_t LocateGlobalPointExclVolume(VPlacedVolume const *vol, VPlacedVolume const *excludedvolume,
                                               Vector3D<Precision> const &point, NavigationState &path, bool top);
VECCORE_ATT_HOST_DEVICE
VPlacedVolumePtr_t RelocatePointFromPathForceDifferent(Vector3D<Precision> const &localpoint, NavigationState &path);
VECCORE_ATT_HOST_DEVICE
double ComputeStepAndPropagatedState(Vector3D<Precision> const &globalpoint, Vector3D<Precision> const &globaldir,
                                     Precision step_limit, NavigationState const &in_state, NavigationState &out_state);

} // End namespace Raytracer

} // End namespace VECGEOM_IMPL_NAMESPACE

void write_ppm(std::string filename, float *buffer, int px, int py);
void RenderCPU(VPlacedVolume const *const world, int px, int py, int maxdepth);

} // End namespace vecgeom

#endif // VECGEOM_RAYTRACER_H_
