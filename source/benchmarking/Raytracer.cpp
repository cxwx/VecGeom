/// \file Raytracer.cpp
/// \author Andrei Gheata (andrei.gheata@cern.ch)

#include "VecGeom/benchmarking/Raytracer.h"

#include <VecGeom/base/Transformation3D.h>
#include <VecGeom/base/Stopwatch.h>

#include <VecGeom/navigation/NavigationState.h>
#include <VecGeom/volumes/PlacedVolume.h>
#include <VecGeom/management/GeoManager.h>

/*
#include <VecGeom/navigation/NewSimpleNavigator.h>
#include <VecGeom/navigation/SimpleABBoxNavigator.h>
#include <VecGeom/navigation/SimpleABBoxLevelLocator.h>
#include <VecGeom/navigation/HybridNavigator2.h>
*/

#ifdef VECGEOM_CUDA_INTERFACE
#include "VecGeom/backend/cuda/Backend.h"
#include "VecGeom/management/CudaManager.h"
#endif

#include <random>
#include <sstream>
#include <fstream>
#include <utility>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

namespace Raytracer {

void InitializeModel(VPlacedVolumePtr_t world, RaytracerData_t &rtdata)
{
  using namespace vecCore::math;

  if (!world) return;
  rtdata.fWorld = world;

  // adjust up vector, image scaling
  Vector3D<double> aMin, aMax, vcenter, vsize;
  world->Extent(aMin, aMax);
  vcenter = 0.5 * (aMin + aMax);
  vsize   = 0.5 * (aMax - aMin);

  double imgRadius = vsize.Mag();
  // std::cout << *(fWorld->GetLogicalVolume()->GetUnplacedVolume()) << std::endl;
  // std::cout << "vcenter =  " << vcenter << "   vsize = " << vsize << "  ingRadius = " << imgRadius << std::endl;
  assert(rtdata.fSize_px * rtdata.fSize_py > 0 && "SetWorld: image size not set");

  // Make sure the image fits the parrallel world view, leaving 20% margin
  constexpr double d0 = 1.;
  double dd           = (vcenter - rtdata.fScreenPos).Mag();
  rtdata.fScale       = 2.1 * imgRadius / Min(rtdata.fSize_px, rtdata.fSize_py) / rtdata.fZoom;
  if (rtdata.fView == kRTVperspective) rtdata.fScale *= d0 / (dd + d0);

  // Project up vector on the source plane
  rtdata.fDir = vcenter - rtdata.fScreenPos;
  rtdata.fDir.Normalize();
  rtdata.fStart = rtdata.fScreenPos - d0 * rtdata.fDir;
  rtdata.fRight = Vector3D<double>::Cross(rtdata.fDir, rtdata.fUp);
  rtdata.fRight.Normalize();
  rtdata.fUp = Vector3D<double>::Cross(rtdata.fRight, rtdata.fDir);
  rtdata.fUp.Normalize();
  rtdata.fLeftC =
      rtdata.fScreenPos - 0.5 * rtdata.fScale * (rtdata.fSize_px * rtdata.fRight + rtdata.fSize_py * rtdata.fUp);

  // Light position on top-left
  rtdata.fSourceDir = rtdata.fDir + rtdata.fUp + rtdata.fRight;
  rtdata.fSourceDir.Normalize();

  // Create navigators (only for CPU case)
  // CreateNavigators();

  // Allocate rays
  rtdata.fNrays = rtdata.fSize_px * rtdata.fSize_py;
}

Color_t RaytraceOne(int px, int py, RaytracerData_t const &rtdata)
{
  constexpr int kMaxTries      = 10;
  constexpr double kPush       = 1.e-8;
  constexpr size_t buffer_size = 2048;             // should contain current/next states
  int maxdepth                 = rtdata.fMaxDepth; /*GeoManager::Instance().getMaxDepth();*/
  size_t navstate_size         = NavigationState::SizeOfInstance(maxdepth);
  // make sure the buffer state is large enough to contain 2 navigation states
  assert(buffer_size > 2 * navstate_size);

  char navstates_buffer[buffer_size];

  Ray_t ray;
  ray.fCrtState  = NavigationState::MakeInstanceAt(maxdepth, (void *)(navstates_buffer));
  ray.fNextState = NavigationState::MakeInstanceAt(maxdepth, (void *)(navstates_buffer + navstate_size));

  Vector3D<double> pos_onscreen = rtdata.fLeftC + rtdata.fScale * (px * rtdata.fRight + py * rtdata.fUp);
  Vector3D<double> start        = (rtdata.fView == kRTVperspective) ? rtdata.fStart : pos_onscreen;
  ray.fPos                      = start;
  ray.fDir                      = (rtdata.fView == kRTVperspective) ? pos_onscreen - rtdata.fStart : rtdata.fDir;
  ray.fDir.Normalize();
  ray.fColor  = 0xFFFFFFFF; // white
  ray.fVolume = (rtdata.fView == kRTVperspective) ? rtdata.fVPstate->Top()
                                                  : LocateGlobalPoint(rtdata.fWorld, ray.fPos, *ray.fCrtState, true);
  int itry = 0;
  while (!ray.fVolume && itry < kMaxTries) {
    auto snext = rtdata.fWorld->DistanceToIn(ray.fPos, ray.fDir);
    ray.fDone  = snext == kInfLength;
    if (ray.fDone) return ray.fColor;
    // Propagate to the world volume (but do not increment the boundary count)
    ray.fPos += (snext + kPush) * ray.fDir;
    ray.fVolume = LocateGlobalPoint(rtdata.fWorld, ray.fPos, *ray.fCrtState, true);
  }
  ray.fDone = ray.fVolume == nullptr;
  if (ray.fDone) return ray.fColor;
  *ray.fNextState = *ray.fCrtState;

  // Now propagate ray
  while (!ray.fDone) {
    auto nextvol = ray.fVolume;
    double snext = kInfLength;
    int nsmall   = 0;

    while (nextvol == ray.fVolume && nsmall < kMaxTries) {
      snext   = ComputeStepAndPropagatedState(ray.fPos, ray.fDir, kInfLength, *ray.fCrtState, *ray.fNextState);
      nextvol = ray.fNextState->Top();
      ray.fPos += (snext + kPush) * ray.fDir;
      nsmall++;
    }
    if (nsmall == kMaxTries) {
      // std::cout << "error for ray (" << px << ", " << py << ")\n";
      ray.fDone  = true;
      ray.fColor = 0;
      return ray.fColor;
    }

    // Apply the selected RT model
    ray.fNcrossed++;
    ray.fVolume = nextvol;
    if (ray.fVolume == nullptr) ray.fDone = true;
    if (nextvol) ApplyRTmodel(ray, snext, rtdata);
    auto tmpstate  = ray.fCrtState;
    ray.fCrtState  = ray.fNextState;
    ray.fNextState = tmpstate;
  }

  return ray.fColor;
}

void ApplyRTmodel(Ray_t &ray, double step, RaytracerData_t const &rtdata)
{
  int depth = ray.fNextState->GetLevel();
  if (rtdata.fModel == kRTspecular) { // specular reflection
    // Calculate normal at the hit point
    bool valid = ray.fVolume != nullptr && depth >= rtdata.fVisDepth;
    if (valid) {
      Transformation3D m;
      ray.fNextState->TopMatrix(m);
      auto localpoint = m.Transform(ray.fPos);
      Vector3D<double> norm, lnorm;
      ray.fVolume->GetLogicalVolume()->GetUnplacedVolume()->Normal(localpoint, lnorm);
      m.InverseTransformDirection(lnorm, norm);
      Vector3D<double> refl = ray.fDir - 2 * norm.Dot(ray.fDir) * norm;
      refl.Normalize();
      double calf = -rtdata.fSourceDir.Dot(refl);
      // if (calf < 0) calf = 0;
      // calf                   = vecCore::math::Pow(calf, fShininess);
      auto specular_color = rtdata.fLightColor;
      specular_color.MultiplyLightChannel(1. + 0.5 * calf);
      auto object_color = rtdata.fObjColor;
      object_color.MultiplyLightChannel(1. + 0.5 * calf);
      ray.fColor = specular_color + object_color;
      ray.fDone  = true;
      // std::cout << "calf = " << calf << "red=" << (int)ray.fColor.fComp.red << " green=" <<
      // (int)ray.fColor.fComp.green
      //          << " blue=" << (int)ray.fColor.fComp.blue << " alpha=" << (int)ray.fColor.fComp.alpha << std::endl;
    }
  } else if (rtdata.fModel == kRTtransparent) { // everything transparent 100% except volumes at visible depth
    bool valid = ray.fVolume != nullptr && depth == rtdata.fVisDepth;
    if (valid) {
      float transparency = 0.85;
      auto object_color  = rtdata.fObjColor;
      object_color *= (1 - transparency);
      ray.fColor += object_color;
    }
  }
  if (ray.fVolume == nullptr) ray.fDone = true;
}

void PropagateRays(RaytracerData_t &rtdata)
{
  // Propagate all rays and write out the image on the CPU
#ifndef VECGEOM_ENABLE_CUDA
  // Setup viewpoint state
  constexpr size_t buffer_size = 1024;             // should contain state for the viewpoint
  int maxdepth                 = rtdata.fMaxDepth; /*GeoManager::Instance().getMaxDepth();*/
  assert(buffer_size > NavigationState::SizeOfInstance(maxdepth));

  char vpstate_buffer[buffer_size];
  auto vpstate = NavigationState::MakeInstanceAt(maxdepth, (void *)(vpstate_buffer));
  // GlobalLocator::LocateGlobalPoint(fWorld, fStart, *vpstate, true);
  LocateGlobalPoint(rtdata.fWorld, rtdata.fStart, *vpstate, true);
  rtdata.fVPstate = vpstate;

  float *buffer = (float *)malloc(4 * rtdata.fNrays * sizeof(float));
  size_t n10    = 0.1 * rtdata.fNrays;
  size_t icrt   = 0;
  // fprintf(stderr, "P3\n%d %d\n255\n", fSize_px, fSize_py);
  for (int py = 0; py < rtdata.fSize_py; py++) {
    for (int px = 0; px < rtdata.fSize_px; px++) {
      if ((icrt % n10) == 0) printf("%lu %%\n", 10 * icrt / n10);
      int pixel_index = 4 * (py * rtdata.fSize_px + px);

      auto pixel_color = RaytraceOne(px, py, rtdata);

      buffer[pixel_index + 0] = pixel_color.Red();
      buffer[pixel_index + 1] = pixel_color.Green();
      buffer[pixel_index + 2] = pixel_color.Blue();
      buffer[pixel_index + 3] = 1.0f;
      icrt++;
    }
  }
  write_ppm("output.ppm", buffer, rtdata.fSize_px, rtdata.fSize_py);
  free(buffer);
#endif
}

///< Explicit navigation functions, we should be using the navigator functionality when it works
VPlacedVolume const *LocateGlobalPoint(VPlacedVolume const *vol, Vector3D<Precision> const &point,
                                       NavigationState &path, bool top)
{
  VPlacedVolume const *candvolume = vol;
  Vector3D<Precision> currentpoint(point);
  if (top) {
    assert(vol != nullptr);
    if (!vol->UnplacedContains(point)) return nullptr;
  }
  path.Push(candvolume);
  LogicalVolume const *lvol         = candvolume->GetLogicalVolume();
  Vector<Daughter> const *daughters = lvol->GetDaughtersp();

  bool godeeper = true;
  while (daughters->size() > 0 && godeeper) {
    for (size_t i = 0; i < daughters->size() && godeeper; ++i) {
      VPlacedVolume const *nextvolume = (*daughters)[i];
      Vector3D<Precision> transformedpoint;
      if (nextvolume->Contains(currentpoint, transformedpoint)) {
        path.Push(nextvolume);
        currentpoint = transformedpoint;
        candvolume   = nextvolume;
        daughters    = candvolume->GetLogicalVolume()->GetDaughtersp();
        break;
      }
    }
    godeeper = false;
  }
  return candvolume;
}

VPlacedVolume const *LocateGlobalPointExclVolume(VPlacedVolume const *vol, VPlacedVolume const *excludedvolume,
                                                 Vector3D<Precision> const &point, NavigationState &path, bool top)
{
  VPlacedVolume const *candvolume = vol;
  Vector3D<Precision> currentpoint(point);
  if (top) {
    assert(vol != nullptr);
    candvolume = (vol->UnplacedContains(point)) ? vol : nullptr;
  }
  if (candvolume) {
    path.Push(candvolume);
    LogicalVolume const *lvol         = candvolume->GetLogicalVolume();
    Vector<Daughter> const *daughters = lvol->GetDaughtersp();

    bool godeeper = true;
    while (daughters->size() > 0 && godeeper) {
      godeeper = false;
      // returns nextvolume; and transformedpoint; modified path
      for (size_t i = 0; i < daughters->size(); ++i) {
        VPlacedVolume const *nextvolume = (*daughters)[i];
        if (nextvolume != excludedvolume) {
          Vector3D<Precision> transformedpoint;
          if (nextvolume->Contains(currentpoint, transformedpoint)) {
            path.Push(nextvolume);
            currentpoint = transformedpoint;
            candvolume   = nextvolume;
            daughters    = candvolume->GetLogicalVolume()->GetDaughtersp();
            godeeper     = true;
            break;
          }
        } // end if excludedvolume
      }
    }
  }
  return candvolume;
}

VPlacedVolume const *RelocatePointFromPathForceDifferent(Vector3D<Precision> const &localpoint, NavigationState &path)
{
  VPlacedVolume const *currentmother = path.Top();
  VPlacedVolume const *entryvol      = currentmother;
  if (currentmother != nullptr) {
    Vector3D<Precision> tmp = localpoint;
    while (currentmother) {
      if (currentmother == entryvol || currentmother->GetLogicalVolume()->GetUnplacedVolume()->IsAssembly() ||
          !currentmother->UnplacedContains(tmp)) {
        path.Pop();
        Vector3D<Precision> pointhigherup = currentmother->GetTransformation()->InverseTransform(tmp);
        tmp                               = pointhigherup;
        currentmother                     = path.Top();
      } else {
        break;
      }
    }

    if (currentmother) {
      path.Pop();
      return LocateGlobalPointExclVolume(currentmother, entryvol, tmp, path, false);
    }
  }
  return currentmother;
}

double ComputeStepAndPropagatedState(Vector3D<Precision> const &globalpoint, Vector3D<Precision> const &globaldir,
                                     Precision step_limit, NavigationState const &in_state, NavigationState &out_state)
{
  // calculate local point/dir from global point/dir
  // call the static function for this provided/specialized by the Impl
  Vector3D<Precision> localpoint;
  Vector3D<Precision> localdir;
  // Impl::DoGlobalToLocalTransformation(in_state, globalpoint, globaldir, localpoint, localdir);
  Transformation3D m;
  in_state.TopMatrix(m);
  localpoint = m.Transform(globalpoint);
  localdir   = m.TransformDirection(globaldir);

  Precision step                    = step_limit;
  VPlacedVolume const *hitcandidate = nullptr;
  auto pvol                         = in_state.Top();
  auto lvol                         = pvol->GetLogicalVolume();

  // need to calc DistanceToOut first
  // step = Impl::TreatDistanceToMother(pvol, localpoint, localdir, step_limit);
  step = pvol->DistanceToOut(localpoint, localdir, step_limit);
  if (step < 0) step = 0;
  // "suck in" algorithm from Impl and treat hit detection in local coordinates for daughters
  //((Impl *)this)
  //    ->Impl::CheckDaughterIntersections(lvol, localpoint, localdir, &in_state, &out_state, step, hitcandidate);
  auto *daughters = lvol->GetDaughtersp();
  auto ndaughters = daughters->size();
  for (decltype(ndaughters) d = 0; d < ndaughters; ++d) {
    auto daughter    = daughters->operator[](d);
    double ddistance = daughter->DistanceToIn(localpoint, localdir, step);

    // if distance is negative; we are inside that daughter and should relocate
    // unless distance is minus infinity
    const bool valid =
        (ddistance < step && !IsInf(ddistance)) && !((ddistance <= 0.) && in_state.GetLastExited() == daughter);
    hitcandidate = valid ? daughter : hitcandidate;
    step         = valid ? ddistance : step;
  }

  // fix state
  bool done;
  // step = Impl::PrepareOutState(in_state, out_state, step, step_limit, hitcandidate, done);
  // now we have the candidates and we prepare the out_state
  in_state.CopyTo(&out_state);
  done = false;
  if (step == kInfLength && step_limit > 0.) {
    step = vecgeom::kTolerance;
    out_state.SetBoundaryState(true);
    do {
      out_state.Pop();
    } while (out_state.Top()->GetLogicalVolume()->GetUnplacedVolume()->IsAssembly());
    done = true;
  } else {
    // is geometry further away than physics step?
    // this is a physics step
    if (step > step_limit) {
      // don't need to do anything
      step = step_limit;
      out_state.SetBoundaryState(false);
    } else {
      // otherwise it is a geometry step
      out_state.SetBoundaryState(true);
      if (hitcandidate) out_state.Push(hitcandidate);

      if (step < 0.) {
        // std::cerr << "WARNING: STEP NEGATIVE; NEXTVOLUME " << nexthitvolume << std::endl;
        // InspectEnvironmentForPointAndDirection( globalpoint, globaldir, currentstate );
        step = 0.;
      }
    }
  }

  if (done) {
    if (out_state.Top() != nullptr) {
      assert(!out_state.Top()->GetLogicalVolume()->GetUnplacedVolume()->IsAssembly());
    }
    return step;
  }
  // step was physics limited
  if (!out_state.IsOnBoundary()) return step;

  // otherwise if necessary do a relocation
  // try relocation to refine out_state to correct location after the boundary

  // ((Impl *)this)->Impl::Relocate(MovePointAfterBoundary(localpoint, localdir, step), in_state, out_state);
  localpoint += (step + 1.E-6) * localdir;

  if (out_state.Top() == in_state.Top()) {
    RelocatePointFromPathForceDifferent(localpoint, out_state);
  } else {
    // continue directly further down ( next volume should have been stored in out_state already )
    VPlacedVolume const *nextvol = out_state.Top();
    out_state.Pop();
    LocateGlobalPoint(nextvol, nextvol->GetTransformation()->Transform(localpoint), out_state, false);
  }

  if (out_state.Top() != nullptr) {
    while (out_state.Top()->IsAssembly()) {
      out_state.Pop();
    }
    assert(!out_state.Top()->GetLogicalVolume()->GetUnplacedVolume()->IsAssembly());
  }
  return step;
}

/*
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
*/
} // End namespace Raytracer
} // End namespace VECGEOM_IMPL_NAMESPACE

#ifndef VECGEOM_CUDA_INTERFACE
void write_ppm(std::string filename, float *buffer, int px, int py)
{
  std::ofstream image(filename);

  image << "P3\n" << px << " " << py << "\n255\n";

  for (int j = py - 1; j >= 0; j--) {
    for (int i = 0; i < px; i++) {
      int idx = 4 * (j * px + i);
      int ir  = int(255.99 * buffer[idx + 0]);
      int ig  = int(255.99 * buffer[idx + 1]);
      int ib  = int(255.99 * buffer[idx + 2]);
      image << ir << " " << ig << " " << ib << "\n";
    }
  }
}
#endif

/// COPY OF GPU CODE COMES BELOW
/// This is needed to have the code in the cxx namespace,
/// as the CPU code in the cuda namespace does not work.

VPlacedVolume const *LocateGlobalPoint(VPlacedVolume const *vol, Vector3D<Precision> const &point,
                                       NavigationState &path, bool top)
{
  VPlacedVolume const *candvolume = vol;
  Vector3D<Precision> currentpoint(point);
  if (top) {
    assert(vol != nullptr);
    if (!vol->UnplacedContains(point)) return nullptr;
  }
  path.Push(candvolume);
  LogicalVolume const *lvol         = candvolume->GetLogicalVolume();
  Vector<Daughter> const *daughters = lvol->GetDaughtersp();

  bool godeeper = true;
  while (daughters->size() > 0 && godeeper) {
    for (size_t i = 0; i < daughters->size() && godeeper; ++i) {
      VPlacedVolume const *nextvolume = (*daughters)[i];
      Vector3D<Precision> transformedpoint;
      if (nextvolume->Contains(currentpoint, transformedpoint)) {
        path.Push(nextvolume);
        currentpoint = transformedpoint;
        candvolume   = nextvolume;
        daughters    = candvolume->GetLogicalVolume()->GetDaughtersp();
        break;
      }
    }
    godeeper = false;
  }
  return candvolume;
}

VPlacedVolume const *LocateGlobalPointExclVolume(VPlacedVolume const *vol, VPlacedVolume const *excludedvolume,
                                                 Vector3D<Precision> const &point, NavigationState &path, bool top)
{
  VPlacedVolume const *candvolume = vol;
  Vector3D<Precision> currentpoint(point);
  if (top) {
    assert(vol != nullptr);
    candvolume = (vol->UnplacedContains(point)) ? vol : nullptr;
  }
  if (candvolume) {
    path.Push(candvolume);
    LogicalVolume const *lvol         = candvolume->GetLogicalVolume();
    Vector<Daughter> const *daughters = lvol->GetDaughtersp();

    bool godeeper = true;
    while (daughters->size() > 0 && godeeper) {
      godeeper = false;
      // returns nextvolume; and transformedpoint; modified path
      for (size_t i = 0; i < daughters->size(); ++i) {
        VPlacedVolume const *nextvolume = (*daughters)[i];
        if (nextvolume != excludedvolume) {
          Vector3D<Precision> transformedpoint;
          if (nextvolume->Contains(currentpoint, transformedpoint)) {
            path.Push(nextvolume);
            currentpoint = transformedpoint;
            candvolume   = nextvolume;
            daughters    = candvolume->GetLogicalVolume()->GetDaughtersp();
            godeeper     = true;
            break;
          }
        } // end if excludedvolume
      }
    }
  }
  return candvolume;
}

VPlacedVolume const *RelocatePointFromPathForceDifferent(Vector3D<Precision> const &localpoint, NavigationState &path)
{
  VPlacedVolume const *currentmother = path.Top();
  VPlacedVolume const *entryvol      = currentmother;
  if (currentmother != nullptr) {
    Vector3D<Precision> tmp = localpoint;
    while (currentmother) {
      if (currentmother == entryvol || currentmother->GetLogicalVolume()->GetUnplacedVolume()->IsAssembly() ||
          !currentmother->UnplacedContains(tmp)) {
        path.Pop();
        Vector3D<Precision> pointhigherup = currentmother->GetTransformation()->InverseTransform(tmp);
        tmp                               = pointhigherup;
        currentmother                     = path.Top();
      } else {
        break;
      }
    }

    if (currentmother) {
      path.Pop();
      return LocateGlobalPointExclVolume(currentmother, entryvol, tmp, path, false);
    }
  }
  return currentmother;
}

Precision ComputeStepAndPropagatedState(Vector3D<Precision> const &globalpoint, Vector3D<Precision> const &globaldir,
                                        Precision step_limit, NavigationState const &in_state,
                                        NavigationState &out_state)
{
  // calculate local point/dir from global point/dir
  // call the static function for this provided/specialized by the Impl
  Vector3D<Precision> localpoint;
  Vector3D<Precision> localdir;
  // Impl::DoGlobalToLocalTransformation(in_state, globalpoint, globaldir, localpoint, localdir);
  Transformation3D m;
  in_state.TopMatrix(m);
  localpoint = m.Transform(globalpoint);
  localdir   = m.TransformDirection(globaldir);

  Precision step                    = step_limit;
  VPlacedVolume const *hitcandidate = nullptr;
  auto pvol                         = in_state.Top();
  auto lvol                         = pvol->GetLogicalVolume();

  // need to calc DistanceToOut first
  // step = Impl::TreatDistanceToMother(pvol, localpoint, localdir, step_limit);
  step = pvol->DistanceToOut(localpoint, localdir, step_limit);
  if (step < 0) step = 0;
  // "suck in" algorithm from Impl and treat hit detection in local coordinates for daughters
  //((Impl *)this)
  //    ->Impl::CheckDaughterIntersections(lvol, localpoint, localdir, &in_state, &out_state, step, hitcandidate);
  auto *daughters = lvol->GetDaughtersp();
  auto ndaughters = daughters->size();
  for (decltype(ndaughters) d = 0; d < ndaughters; ++d) {
    auto daughter       = daughters->operator[](d);
    Precision ddistance = daughter->DistanceToIn(localpoint, localdir, step);

    // if distance is negative; we are inside that daughter and should relocate
    // unless distance is minus infinity
    const bool valid =
        (ddistance < step && !IsInf(ddistance)) && !((ddistance <= 0.) && in_state.GetLastExited() == daughter);
    hitcandidate = valid ? daughter : hitcandidate;
    step         = valid ? ddistance : step;
  }

  // fix state
  bool done;
  // step = Impl::PrepareOutState(in_state, out_state, step, step_limit, hitcandidate, done);
  // now we have the candidates and we prepare the out_state
  in_state.CopyTo(&out_state);
  done = false;
  if (step == kInfLength && step_limit > 0.) {
    step = vecgeom::kTolerance;
    out_state.SetBoundaryState(true);
    do {
      out_state.Pop();
    } while (out_state.Top()->GetLogicalVolume()->GetUnplacedVolume()->IsAssembly());
    done = true;
  } else {
    // is geometry further away than physics step?
    // this is a physics step
    if (step > step_limit) {
      // don't need to do anything
      step = step_limit;
      out_state.SetBoundaryState(false);
    } else {
      // otherwise it is a geometry step
      out_state.SetBoundaryState(true);
      if (hitcandidate) out_state.Push(hitcandidate);

      if (step < 0.) {
        // std::cerr << "WARNING: STEP NEGATIVE; NEXTVOLUME " << nexthitvolume << std::endl;
        // InspectEnvironmentForPointAndDirection( globalpoint, globaldir, currentstate );
        step = 0.;
      }
    }
  }

  if (done) {
    if (out_state.Top() != nullptr) {
      assert(!out_state.Top()->GetLogicalVolume()->GetUnplacedVolume()->IsAssembly());
    }
    return step;
  }
  // step was physics limited
  if (!out_state.IsOnBoundary()) return step;

  // otherwise if necessary do a relocation
  // try relocation to refine out_state to correct location after the boundary

  // ((Impl *)this)->Impl::Relocate(MovePointAfterBoundary(localpoint, localdir, step), in_state, out_state);
  localpoint += (step + 1.E-6) * localdir;

  if (out_state.Top() == in_state.Top()) {
    RelocatePointFromPathForceDifferent(localpoint, out_state);
  } else {
    // continue directly further down ( next volume should have been stored in out_state already )
    VPlacedVolume const *nextvol = out_state.Top();
    out_state.Pop();
    LocateGlobalPoint(nextvol, nextvol->GetTransformation()->Transform(localpoint), out_state, false);
  }

  if (out_state.Top() != nullptr) {
    while (out_state.Top()->IsAssembly()) {
      out_state.Pop();
    }
    assert(!out_state.Top()->GetLogicalVolume()->GetUnplacedVolume()->IsAssembly());
  }
  return step;
}

Color_t raytrace(VPlacedVolume const *const world, Vector3D<Precision> origin, Vector3D<Precision> dir, int maxdepth)
{
  Ray_t ray;
  char cstate[512];
  char nstate[512];

  ray.fPos       = origin;
  ray.fDir       = dir;
  ray.fColor     = 0xFFFFFFFF;
  ray.fCrtState  = NavigationState::MakeInstanceAt(maxdepth, (void *)(cstate));
  ray.fNextState = NavigationState::MakeInstanceAt(maxdepth, (void *)(nstate));

  ray.fVolume = LocateGlobalPoint(world, ray.fPos, *ray.fCrtState, true);

  do {
    auto t = world->DistanceToIn(ray.fPos, ray.fDir);

    if (t == kInfLength) return ray.fColor;

    ray.fPos += (t + 1e-7) * ray.fDir;
    ray.fVolume = LocateGlobalPoint(world, ray.fPos, *ray.fCrtState, true);

  } while (!ray.fVolume);

  if (!ray.fVolume) return ray.fColor;

  *ray.fNextState = *ray.fCrtState;

  while (ray.fVolume) {
    auto t = ComputeStepAndPropagatedState(ray.fPos, ray.fDir, kInfLength, *ray.fCrtState, *ray.fNextState);

    if (t == kInfLength || !ray.fNextState->Top()) break;

    ray.fPos += (t + 1e-7) * ray.fDir;
    ray.fVolume    = ray.fNextState->Top();
    *ray.fCrtState = *ray.fNextState;

    if (ray.fVolume->GetDaughters().size() == 0) ray.fColor += 0x0000ff15;
  }

  return ray.fColor;
}

void RenderCPU(VPlacedVolume const *const world, int px, int py, int maxdepth)
{
  float *buffer = (float *)malloc(4 * px * py * sizeof(float));

  for (int i = 0; i < px; ++i) {
    for (int j = 0; j < py; ++j) {
      int pixel_index = 4 * (j * px + i);

      float u = float(i) / float(px);
      float v = float(j) / float(py);

      // model view hard-coded for debugging
      // traceML size is 2200,2200,6200, centered at 0,0,0
      Vector3D<Precision> origin    = {0, -7000, 0};
      Vector3D<Precision> direction = {v - 0.5, 1.9, 2 * u - 1};
      direction.Normalize();

      Color_t pixel_color = raytrace(world, origin, direction, maxdepth);

      buffer[pixel_index + 0] = pixel_color.Red();
      buffer[pixel_index + 1] = pixel_color.Green();
      buffer[pixel_index + 2] = pixel_color.Blue();
      buffer[pixel_index + 3] = 1.0f;
    }
  }

  write_ppm("output.ppm", buffer, px, py);

  free(buffer);
}
} // End namespace vecgeom
