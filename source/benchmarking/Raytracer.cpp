/// \file Raytracer.cpp
/// \author Andrei Gheata (andrei.gheata@cern.ch)

#include "VecGeom/benchmarking/Raytracer.h"

#include <VecGeom/base/Transformation3D.h>
#include <VecGeom/base/Stopwatch.h>

#include <VecGeom/navigation/NavigationState.h>
#include <VecGeom/volumes/PlacedVolume.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/navigation/NewSimpleNavigator.h>
#include <VecGeom/navigation/SimpleABBoxNavigator.h>
#include <VecGeom/navigation/SimpleABBoxLevelLocator.h>
#include <VecGeom/navigation/HybridNavigator2.h>

#ifdef VECGEOM_CUDA_INTERFACE
#include "VecGeom/backend/cuda/Backend.h"
#include "VecGeom/management/CudaManager.h"
#endif

#include <random>
#include <sstream>
#include <utility>

namespace vecgeom {

Raytracer::Raytracer(VPlacedVolumePtr_t world, Vector3D<double> const &screen_position,
                     Vector3D<double> const &up_vector, int img_size_px, int img_size_py, ERTmodel model, ERTView view)
    : fScreenPos(screen_position), fUp(up_vector), fSize_px(img_size_px), fSize_py(img_size_py), fModel(model),
      fView(view)
{
  SetWorld(world);
}

Raytracer::~Raytracer()
{
  delete[] fRays;
  delete[] fNavStates;
  NavigationState::ReleaseInstance(fVPstate);
}

void Raytracer::SetWorld(VPlacedVolumePtr_t world)
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
  fScale = 1.2 * imgRadius / Min(fSize_px, fSize_py);
  if (fView == kRTVperspective) fScale *= 0.5;

  // Project up vector on the source plane
  fDir   = vcenter - fScreenPos;
  fStart = fScreenPos - fDir;
  fDir.Normalize();
  fRight = Vector3D<double>::Cross(fDir, fUp);
  fRight.Normalize();
  fUp = Vector3D<double>::Cross(fRight, fDir);
  fUp.Normalize();
  fLeftC = fScreenPos - 0.5 * fScale * (fSize_px * fRight + fSize_py * fUp);

  // Light position on top-left
  fSourceDir = fDir + fUp + fRight;
  fSourceDir.Normalize();

  // Create navigators (only for CPU case)
  CreateNavigators();

  // Setup viewpoint state
  int maxdepth = GeoManager::Instance().getMaxDepth();
  fVPstate     = NavigationState::MakeInstance(maxdepth);
  GlobalLocator::LocateGlobalPoint(fWorld, fStart, *fVPstate, true);

  // Allocate rays
  fNrays                = fSize_px * fSize_py;
  fRays                 = new Ray_t[fNrays];
  size_t nav_state_size = NavigationState::SizeOfInstance(maxdepth);
  fNavStates            = new char[2 * nav_state_size * fNrays];
  for (auto i = 0; i < fNrays; ++i) {
    fRays[i].fCrtState = NavigationState::MakeInstanceAt(maxdepth, (void *)(fNavStates + 2 * i * nav_state_size));
    fRays[i].fNextState =
        NavigationState::MakeInstanceAt(maxdepth, (void *)(fNavStates + (2 * i + 1) * nav_state_size));
  }
}

void Raytracer::CreateNavigators()
{
  // Create all navigators.
  for (auto &lvol : vecgeom::GeoManager::Instance().GetLogicalVolumesMap()) {
    if (lvol.second->GetDaughtersp()->size() < 4) {
      lvol.second->SetNavigator(vecgeom::NewSimpleNavigator<>::Instance());
    }
    if (lvol.second->GetDaughtersp()->size() >= 5) {
      lvol.second->SetNavigator(vecgeom::SimpleABBoxNavigator<>::Instance());
    }
    if (lvol.second->GetDaughtersp()->size() >= 10) {
      lvol.second->SetNavigator(vecgeom::HybridNavigator<>::Instance());
      vecgeom::HybridManager2::Instance().InitStructure((lvol.second));
    }
    lvol.second->SetLevelLocator(vecgeom::SimpleABBoxLevelLocator::GetInstance());
  }
}

void Raytracer::GenerateVolumePointers(VPlacedVolumePtr_t vol)
{
  for (auto vd : vol->GetDaughters()) {

    fVolumes.push_back(vd);
    // can check now the property of the conversions of *i
    GenerateVolumePointers(vd);
  }
}

void Raytracer::StartRay(int iray)
{
  constexpr int kMaxTries       = 10;
  constexpr double kPush        = 1.e-8;
  int px                        = iray / fSize_py;
  int py                        = iray - px * fSize_py;
  Vector3D<double> pos_onscreen = fLeftC + fScale * (px * fRight + py * fUp);
  Vector3D<double> start        = (fView == kRTVperspective) ? fStart : pos_onscreen;
  // Locate starting point
  Ray_t &ray = fRays[px * fSize_py + py];
  ray.fPos   = start;
  ray.fDir   = (fView == kRTVperspective) ? pos_onscreen - fStart : fDir;
  ray.fDir.Normalize();
  ray.fColor  = 0xFFFFFFFF;
  ray.fVolume = (fView == kRTVperspective) ? fVPstate->Top()
                                           : GlobalLocator::LocateGlobalPoint(fWorld, ray.fPos, *ray.fCrtState, true);
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
  int depth = ray.fNextState->GetLevel();
  if (fModel == kRTspecular) { // specular reflection
    // Calculate normal at the hit point
    bool valid = ray.fVolume != nullptr && depth >= fVisDepth;
    if (valid) {
      Transformation3D m;
      ray.fNextState->TopMatrix(m);
      auto localpoint = m.Transform(ray.fPos);
      Vector3D<double> norm, lnorm;
      ray.fVolume->GetLogicalVolume()->GetUnplacedVolume()->Normal(localpoint, lnorm);
      m.InverseTransformDirection(lnorm, norm);
      Vector3D<double> refl = fSourceDir - 2 * norm.Dot(fSourceDir) * norm;
      refl.Normalize();
      double calf = -fDir.Dot(refl);
      // if (calf < 0) calf = 0;
      // calf                   = vecCore::math::Pow(calf, fShininess);
      auto specular_color = fLightColor;
      specular_color.MultiplyLightChannel(1. + 0.7 * calf);
      auto object_color = fObjColor;
      object_color.MultiplyLightChannel(1. + 0.7 * calf);
      ray.fColor = specular_color + object_color;
      ray.fDone  = true;
      // std::cout << "calf = " << calf << "red=" << (int)ray.fColor.fComp.red << " green=" <<
      // (int)ray.fColor.fComp.green
      //          << " blue=" << (int)ray.fColor.fComp.blue << " alpha=" << (int)ray.fColor.fComp.alpha << std::endl;
    }
  } else if (fModel == kRTtransparent) { // everything transparent 100% except volumes at visible depth
    bool valid = ray.fVolume != nullptr && depth == fVisDepth;
    if (valid) {
      float transparency = 0.85;
      auto object_color  = fObjColor;
      object_color *= (1 - transparency);
      ray.fColor += object_color;
    }
  }
  if (ray.fVolume == nullptr) ray.fDone = true;
}

void Raytracer::PropagateRays()
{
  for (int iray = 0; iray < fNrays; ++iray) {
    StartRay(iray);
  }
  int remaining;
  while ((remaining = PropagateAllOneStep())) {
    std::cout << "remaining " << remaining << std::endl;
  }

  // Write out the image
  fprintf(stderr, "P3\n%d %d\n255\n", fSize_px, fSize_py);
  for (int j = 0; j < fSize_py; j++) {
    for (int i = 0; i < fSize_px; i++) {
      size_t pixel_index = i * fSize_py + j;
      auto col           = fRays[pixel_index].fColor;
      int red            = 255 * col.Red();
      int green          = 255 * col.Green();
      int blue           = 255 * col.Blue();
      fprintf(stderr, "%d %d %d\n", red, green, blue);
    }
  }
}

int Raytracer::PropagateAllOneStep()
{
  int count_alive = 0;
  // dummy strategy: loop over all rays and propagate the ones not done
  // This will have tails since not all rays finish at the same time. The
  // Alternative is to store rays alive in a container/queue
  for (int iray = 0; iray < fNrays; ++iray) {
    if (fRays[iray].fDone) {
      // std::cout << "ray " << iray << ": red=" << (int)fRays[iray].fColor.fComp.red << " green=" <<
      // (int)fRays[iray].fColor.fComp.green
      //          << " blue=" << (int)fRays[iray].fColor.fComp.blue << " alpha=" << (int)fRays[iray].fColor.fComp.alpha
      //          << std::endl;
      continue;
    }
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
  // Vector3D<double> start  = fLeftC + fScale * (px * fRight + py * fUp);
  Ray_t &ray              = fRays[px * fSize_py + py];

  auto nav = ray.fVolume->GetLogicalVolume()->GetNavigator();
  // auto nav     = static_cast<NewSimpleNavigator<false> *>(NewSimpleNavigator<false>::Instance());
  auto nextvol = ray.fVolume;
  double snext = kInfLength;
  int nsmall   = 0;

  while (nextvol == ray.fVolume && nsmall < kMaxTries) {
    // if (px > 500 && px < 502 && py > 500 && py < 502) {
    //  std::cout << "ray " << iray << "\n";
    //}
    snext   = nav->ComputeStepAndPropagatedState(ray.fPos, ray.fDir, kInfLength, *ray.fCrtState, *ray.fNextState);
    nextvol = ray.fNextState->Top();
    ray.fPos += (snext + kPush) * ray.fDir;
    nsmall++;
  }
  if (nsmall == kMaxTries) {
    std::cout << "error for ray " << iray << std::endl;
    ray.fDone  = true;
    ray.fColor = 0;
    return 0;
  }

  if (px > fSize_px / 2 + 200 && px < fSize_px / 2 + 202 && py > fSize_py / 2 && py < fSize_py / 2 + 2) {
    ray.fColor = 0xFF0000FF;
    std::cout << "ray " << iray << " currently in: ";
    ray.fCrtState->Print();
    std::cout << " entering at step = " << snext << ": ";
    ray.fNextState->Print();

    // ray.fDone = true;
    // return 0;
  }

  // Apply the selected RT model
  ray.fNcrossed++;
  ray.fVolume = nextvol;
  if (ray.fVolume == nullptr) ray.fDone = true;
  if (nextvol) ApplyRTmodel(ray, snext);
  auto tmpstate  = ray.fCrtState;
  ray.fCrtState  = ray.fNextState;
  ray.fNextState = tmpstate;

  return (ray.fDone) ? 0 : 1;
}

#ifdef VECGEOM_CUDA_INTERFACE
void Raytracer::GetVolumePointers(std::list<DevicePtr<cuda::VPlacedVolume>> &volumesGpu)
{
  CudaManager::Instance().LoadGeometry(GetWorld());
  CudaManager::Instance().Synchronize();
  for (auto v : fVolumes) {
    volumesGpu.push_back(CudaManager::Instance().LookupPlaced(v));
  }
}
#endif

} // End namespace vecgeom
