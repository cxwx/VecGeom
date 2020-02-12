/// \file Raytracer.cpp
/// \author Andrei Gheata (andrei.gheata@cern.ch)

#include "VecGeom/benchmarking/Raytracer.h"

#include "VecGeom/base/Stopwatch.h"

#ifdef VECGEOM_CUDA_INTERFACE
#include "VecGeom/backend/cuda/Backend.h"
#include "VecGeom/management/CudaManager.h"
#endif

#include <random>
#include <sstream>
#include <utility>

namespace vecgeom {
/*
Raytracer::Raytracer(VPlacedVolume const *world, Vector3D<double> const &source_position,  Vector3D<double> const
&up_vector, int img_size_px, int img_size_py) : fSize_px(img_size_px), fSize_py(img_size_py), fStart(camera_position),
fUp(up_vector)
{
  SetWorld(world);
}
*/
Raytracer::~Raytracer()
{
  delete[] fRays;
  delete[] fNavStates;
}

void Raytracer::SetWorld(VPlacedVolume_t world)
{
  using namespace vecCore::math;

  // fVolumes.clear();
  fWorld = world;
  if (!world) return;
  GenerateVolumePointers(fWorld);
  if (fVerbosity > 2) {
    // printf("Found %lu volumes in world volume to be used for benchmarking.\n", fVolumes.size());
  }
  // adjust up vector, image scaling
  Vector3D<double> aMin, aMax, vcenter, vsize;
  fWorld->Extent(aMin, aMax);
  vcenter = 0.5 * (aMin + aMax);
  vsize   = 0.5 * (aMax - aMin);

  double imgRadius = vsize.Mag();
  assert(fSize_px * fSize_py > 0 && "SetWorld: image size not set");

  // Make sure the image fits the parrallel world view, leaving 20% margin
  fScale = 0.8 * imgRadius / Min(fSize_px, fSize_py);

  // Project up vector on the source plane
  fDir = vcenter - fStart;
  fDir.Normalize();
  fRight = Vector3D<double>::Cross(fUp, fDir);
  fRight.Normalize();
  fUp = Vector3D<double>::Cross(fDir, fRight);
  fUp.Normalize();
  fLeftC = fStart - 0.5 * fScale * (fSize_px * fRight + fSize_py * fUp);

  // Calculate light source position for the specular model
  if (fModel == kRTspecular) {
    if (fSourceDir.Mag2() < kTolerance) {
      // Light position on top-left
      Vector3D<double> ptstart = fStart - fRight;
      Vector3D<double> ptend   = fDir;
      fSourceDir               = ptend - ptstart;
      fSourceDir.Normalize();
    }
  }
  // Allocate rays
  fNrays                = fSize_px * fSize_py;
  fRays                 = new Ray_t[fNrays];
  int maxdepth          = GeoManager::Instance().getMaxDepth();
  size_t nav_state_size = NavigationState::SizeOfInstance(maxdepth);
  fNavStates            = new char[2 * nav_state_size * fNrays];
  for (auto i = 0; i < fNrays; ++i) {
    fRays[i].fCrtState = NavigationState::MakeInstanceAt(maxdepth, (void *)(fNavStates + 2 * i * nav_state_size));
    fRays[i].fNextState =
        NavigationState::MakeInstanceAt(maxdepth, (void *)(fNavStates + (2 * i + 1) * nav_state_size));
  }
}

void Raytracer::GenerateVolumePointers(VPlacedVolume_t vol)
{
  for (auto vd : vol->GetDaughters()) {

    fVolumes.push_back(vd);
    // can check now the property of the conversions of *i
    GenerateVolumePointers(vd);
  }
}

void Raytracer::StartRay(int iray)
{
  constexpr int kMaxTries = 10;
  constexpr double kPush  = 1.e-8;
  int px                  = iray / fSize_py;
  int py                  = iray - px * fSize_py;
  Vector3D<double> start  = fLeftC + fScale * (px * fRight + py * fUp);
  // Locate starting point
  Ray_t &ray  = fRays[px * fSize_py + py];
  ray.fVolume = GlobalLocator::LocateGlobalPoint(fWorld, ray.fPos, *ray.fCrtState, true);
  // Special case when starting point is outside the setup
  int itry = 0;
  while (!ray.fVolume && itry < kMaxTries) {
    auto snext = fWorld->DistanceToIn(ray.fPos, ray.fDir);
    ray.fDone  = snext == kInfLength;
    if (ray.fDone) return;
    // Propagate to the world volume (but do not increment the boundary count)
    ray.fPos += (snext + kPush) * ray.fDir;
    ray.fVolume = GlobalLocator::LocateGlobalPoint(fWorld, ray.fPos, *ray.fCrtState, true);
  }
  *ray.fNextState = *ray.fCrtState;
  ray.fDone       = ray.fVolume == nullptr;
}

void Raytracer::ApplyRTmodel(Ray_t &ray, double step)
{
  if (fModel == kRTspecular) { // specular reflection
    // Calculate normal at the hit point
    Transformation3D m;
    ray.fNextState->TopMatrix(m);
    auto localpoint = m.Transform(ray.fPos);
    Vector3D<double> norm, lnorm;
    bool valid =
        (ray.fVolume != nullptr) && ray.fVolume->GetLogicalVolume()->GetUnplacedVolume()->Normal(localpoint, lnorm);
    if (valid) {
      m.InverseTransformDirection(lnorm, norm);
      Vector3D<double> refl = 2 * norm.Dot(fSourceDir) - fSourceDir;
      refl.Normalize();
      double calf = -fDir.Dot(refl);
      if (calf < 0) calf = 0;
      calf                   = vecCore::math::Pow(calf, fShininess);
      Color_t specular_color = fLightColor;
      specular_color *= calf;
      auto object_color = fObjColor;
      object_color *= (1. - calf);
      ray.fColor = specular_color + object_color;
    }
    ray.fDone = true;
  }
}

void Raytracer::PropagateRays()
{
  for (int iray = 0; iray < fNrays; ++iray) {
    StartRay(iray);
  }

  while (PropagateAllOneStep()) {
  }
}

int Raytracer::PropagateAllOneStep()
{
  int count_alive = 0;
  // dummy strategy: loop over all rays and propagate the ones not done
  // This will have tails since not all rays finish at the same time. The
  // Alternative is to store rays alive in a container/queue
  for (int iray = 0; iray < fNrays; ++iray) {
    if (fRays[iray].fDone) continue;
    count_alive += PropagateOneStep(iray);
  }
  return count_alive;
}

int Raytracer::PropagateOneStep(int iray)
{
  constexpr double kPush  = 1.e-8;
  constexpr int kMaxTries = 10;
  int px                  = iray / fSize_py;
  int py                  = iray - px * fSize_py;
  Vector3D<double> start  = fLeftC + fScale * (px * fRight + py * fUp);
  Ray_t &ray              = fRays[px * fSize_py + py];

  auto nav     = dynamic_cast<NewSimpleNavigator<false> *>(NewSimpleNavigator<false>::Instance());
  auto nextvol = ray.fVolume;
  double snext = kInfLength;
  int nsmall   = 0;
  while (nextvol == ray.fVolume && nsmall < kMaxTries) {
    snext   = nav->ComputeStepAndPropagatedState(ray.fPos, ray.fDir, kInfLength, *ray.fCrtState, *ray.fNextState);
    nextvol = ray.fNextState->Top();
    ray.fPos += (snext + kPush) * ray.fDir;
    nsmall++;
  }
  if (nsmall == kMaxTries) {
    ray.fDone  = true;
    ray.fColor = 0;
    return 0;
  }
  // Apply the selected RT model
  ray.fNcrossed++;
  ray.fVolume = nextvol;
  if (ray.fVolume == nullptr) ray.fDone = true;
  ApplyRTmodel(ray, snext);
  return (ray.fDone) ? 0 : 1;
}

#ifdef VECGEOM_CUDA_INTERFACE
void Raytracer::GetVolumePointers(std::list<DevicePtr<cuda::VPlacedVolume>> &volumesGpu)
{
  CudaManager::Instance().LoadGeometry(GetWorld());
  CudaManager::Instance().Synchronize();
  for (std::list<VolumePointers>::const_iterator v = fVolumes.begin(); v != fVolumes.end(); ++v) {
    volumesGpu.push_back(CudaManager::Instance().LookupPlaced(v->Specialized()));
  }
}
#endif

} // End namespace vecgeom
