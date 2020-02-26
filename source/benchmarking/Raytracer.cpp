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
  //CreateNavigators();

  // Setup viewpoint state
  int maxdepth = GeoManager::Instance().getMaxDepth();
  fVPstate     = NavigationState::MakeInstance(maxdepth);
  //GlobalLocator::LocateGlobalPoint(fWorld, fStart, *fVPstate, true);
  LocateGlobalPoint(fWorld, fStart, *fVPstate, true);

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

VPlacedVolume const *Raytracer::LocateGlobalPoint(VPlacedVolume const *vol, Vector3D<Precision> const &point,
                                                  NavigationState &path, bool top) const
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

VPlacedVolume const *Raytracer::LocateGlobalPointExclVolume(VPlacedVolume const *vol, VPlacedVolume const *excludedvolume,
                                                            Vector3D<Precision> const &point, NavigationState &path, bool top) const
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

VPlacedVolume const *Raytracer::RelocatePointFromPathForceDifferent(Vector3D<Precision> const &localpoint,
                                                                    NavigationState &path) const
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

double Raytracer::ComputeStepAndPropagatedState(Vector3D<Precision> const &globalpoint,
                                                Vector3D<Precision> const &globaldir, Precision step_limit,
                                                NavigationState const &in_state,
                                                NavigationState &out_state) const
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
    auto daughter = daughters->operator[](d);
    double ddistance = daughter->DistanceToIn(localpoint, localdir, step);

    // if distance is negative; we are inside that daughter and should relocate
    // unless distance is minus infinity
    const bool valid = (ddistance < step && !IsInf(ddistance)) &&
                       !((ddistance <= 0.) && in_state.GetLastExited() == daughter);
    hitcandidate = valid ? daughter : hitcandidate;
    step         = valid ? ddistance : step;
  }

  // fix state
  bool done;
  //step = Impl::PrepareOutState(in_state, out_state, step, step_limit, hitcandidate, done);
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
*/

void Raytracer::GenerateVolumePointers(VPlacedVolumePtr_t vol)
{
  for (auto vd : vol->GetDaughters()) {

    fVolumes.push_back(vd);
    // can check now the property of the conversions of *i
    GenerateVolumePointers(vd);
  }
}

int Raytracer::RaytraceOne(int px, int py, VPlacedVolumePtr_t world, ERTmodel model, ERTView view,
                  Vector3D<double> start, Vector3D<double> dir, Vector3D<double> leftc, Vector3D<double> up, Vector3D<double> right, double scale,
                  Vector3D<double> source_dir, NavigationState *vpstate)
{
  constexpr size_t buffer_size = 2048; // should contain current/next states
  int maxdepth = GeoManager::Instance().getMaxDepth();
  size_t navstate_size = NavigationState::SizeOfInstance(maxdepth);

  char navstates_buffer[buffer_size];

  Ray_t ray;
  ray.fCrtState = NavigationState::MakeInstanceAt(maxdepth, (void *)(navstates_buffer));
  ray.fNextState = NavigationState::MakeInstanceAt(maxdepth, (void *)(navstates_buffer + navstate_size));


  return 0;
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
                                           : /*GlobalLocator::*/LocateGlobalPoint(fWorld, ray.fPos, *ray.fCrtState, true);
  // Special case when starting point is outside the setup
  int itry = 0;
  while (!ray.fVolume && itry < kMaxTries) {
    auto snext = fWorld->DistanceToIn(ray.fPos, ray.fDir);
    ray.fDone  = snext == kInfLength;
    if (ray.fDone) return;
    // Propagate to the world volume (but do not increment the boundary count)
    ray.fPos += (snext + kPush) * ray.fDir;
    ray.fVolume = /*GlobalLocator::*/LocateGlobalPoint(fWorld, ray.fPos, *ray.fCrtState, true);
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

  // auto nav = ray.fVolume->GetLogicalVolume()->GetNavigator();
  // auto nav     = static_cast<NewSimpleNavigator<false> *>(NewSimpleNavigator<false>::Instance());
  auto nextvol = ray.fVolume;
  double snext = kInfLength;
  int nsmall   = 0;

  while (nextvol == ray.fVolume && nsmall < kMaxTries) {
    // if (px > 500 && px < 502 && py > 500 && py < 502) {
    //  std::cout << "ray " << iray << "\n";
    //}
    snext   = /*nav->*/ComputeStepAndPropagatedState(ray.fPos, ray.fDir, kInfLength, *ray.fCrtState, *ray.fNextState);
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
