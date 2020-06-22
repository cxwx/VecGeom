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

/**
 * @brief Rounds up an address to the aligned value
 * @param buf Buffer address to align
 */
VECCORE_ATT_HOST_DEVICE
static char *round_up_align(char *buf)
{
  long remainder = ((long)buf) % 64;
  if (remainder == 0) return buf;
  return (buf + 64 - remainder);
}

/**
 * @brief Rounds up a value to upper aligned version
 * @param buf Buffer address to align
 */
VECCORE_ATT_HOST_DEVICE
static size_t round_up_align(size_t value)
{
  size_t remainder = ((size_t)value) % 64;
  if (remainder == 0) return value;
  return (value + 64 - remainder);
}

inline namespace VECGEOM_IMPL_NAMESPACE {

Ray_t::Ray_t(void *addr, int maxdepth) : fMaxDepth(maxdepth)
{
  char *path_addr = round_up_align((char *)addr + sizeof(Ray_t));
  // Geometry paths follow
  fCrtState       = NavigationState::MakeInstanceAt(maxdepth, path_addr);
  path_addr      += round_up_align(NavigationState::SizeOfInstance(maxdepth));
  fNextState      = NavigationState::MakeInstanceAt(maxdepth, path_addr);
}

void Ray_t::FixGPUpointers()
{
  // Rays are copied from host to device in a buffer, containing still the CPU pointers, which are invalid
  // Since the states are following the regular Ray_t data, we just need to recompute the pointers (when on the device)
  char *path_addr = round_up_align((char *)this + sizeof(Ray_t));
  fCrtState = reinterpret_cast<NavigationState*>(path_addr);
  path_addr      += round_up_align(NavigationState::SizeOfInstance(fMaxDepth));
  fNextState = reinterpret_cast<NavigationState*>(path_addr);
}

size_t Ray_t::SizeOfInstance(int maxdepth)
{
  size_t size = sizeof(Ray_t) + 2 * round_up_align(NavigationState::SizeOfInstance(maxdepth));
  return size;
}

void RaytracerData_t::Print()
{
  printf("  screen_pos(%g, %g, %g) screen_size(%d, %d)\n", fScreenPos[0], fScreenPos[1], fScreenPos[2], fSize_px, fSize_py);
  printf("  light_dir(%g, %g, %g) light_color(0x%08x) obj_color(0x%08x)\n", fSourceDir[0], fSourceDir[1], fSourceDir[2], fLightColor.fColor, fObjColor.fColor);
  printf("  zoom_factor(%g) visible_depth(%d/%d) rt_model(%d) rt_view(%d)\n", fZoom, fVisDepth, fMaxDepth, (int)fModel, (int)fView);
  printf("  viewpoint_state: ");
  if (fVPstate) fVPstate->Print();
  else printf("nullptr\n");
}

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

Color_t RaytraceOne(int px, int py, RaytracerData_t const &rtdata, void *input_buffer)
{
  constexpr int kMaxTries      = 10;
  constexpr double kPush       = 1.e-8;  
  int maxdepth                 = rtdata.fMaxDepth;

  //if (px == rtdata.fSize_px/2 && py == rtdata.fSize_py/2) {
  //  printf("px=%d  py=%d\n", px, py);
  //}

  size_t statesize = NavigationState::SizeOfInstance(maxdepth);
  size_t raysize   = Ray_t::SizeOfInstance(maxdepth);
  char *raybuff    = (char*)input_buffer + statesize;

  int ray_index = py * rtdata.fSize_px + px;
  Ray_t *ray = (Ray_t*)(raybuff + ray_index * raysize);

  Vector3D<double> pos_onscreen = rtdata.fLeftC + rtdata.fScale * (px * rtdata.fRight + py * rtdata.fUp);
  Vector3D<double> start        = (rtdata.fView == kRTVperspective) ? rtdata.fStart : pos_onscreen;
  ray->fPos                     = start;
  ray->fDir                     = (rtdata.fView == kRTVperspective) ? pos_onscreen - rtdata.fStart : rtdata.fDir;
  ray->fDir.Normalize();
  ray->fColor  = 0xFFFFFFFF; // white
  ray->fVolume = (rtdata.fView == kRTVperspective) ? rtdata.fVPstate->Top()
                                                  : LocateGlobalPoint(rtdata.fWorld, ray->fPos, *ray->fCrtState, true);
  int itry = 0;
  while (!ray->fVolume && itry < kMaxTries) {
    auto snext = rtdata.fWorld->DistanceToIn(ray->fPos, ray->fDir);
    ray->fDone  = snext == kInfLength;
    if (ray->fDone) return ray->fColor;
    // Propagate to the world volume (but do not increment the boundary count)
    ray->fPos += (snext + kPush) * ray->fDir;
    ray->fVolume = LocateGlobalPoint(rtdata.fWorld, ray->fPos, *ray->fCrtState, true);
  }
  ray->fDone = ray->fVolume == nullptr;
  if (ray->fDone) return ray->fColor;
  *ray->fNextState = *ray->fCrtState;

  // Now propagate ray
  while (!ray->fDone) {
    auto nextvol = ray->fVolume;
    double snext = kInfLength;
    int nsmall   = 0;

    while (nextvol == ray->fVolume && nsmall < kMaxTries) {
      snext   = ComputeStepAndPropagatedState(ray->fPos, ray->fDir, kInfLength, *ray->fCrtState, *ray->fNextState);
      nextvol = ray->fNextState->Top();
      ray->fPos += (snext + kPush) * ray->fDir;
      nsmall++;
    }
    if (nsmall == kMaxTries) {
      // std::cout << "error for ray (" << px << ", " << py << ")\n";
      ray->fDone  = true;
      ray->fColor = 0;
      return ray->fColor;
    }

    // Apply the selected RT model
    ray->fNcrossed++;
    ray->fVolume = nextvol;
    if (ray->fVolume == nullptr) ray->fDone = true;
    if (nextvol) ApplyRTmodel(*ray, snext, rtdata);
    auto tmpstate  = ray->fCrtState;
    ray->fCrtState  = ray->fNextState;
    ray->fNextState = tmpstate;
  }

  return ray->fColor;
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

void PropagateRays(RaytracerData_t &rtdata, void *input_buffer, void *output_buffer)
{
  // Propagate all rays and write out the image on the CPU
#ifdef VECGEOM_ENABLE_CUDA
  printf("Cannot run CPU propagation kernel on GPU\n");
  return;
#endif

  // The viewpoint state is the first in the buffer
  NavigationState *vpstate = reinterpret_cast<NavigationState*>(input_buffer);
  rtdata.fVPstate = vpstate;

  auto buffer = (unsigned char *)output_buffer;
  size_t n10    = 0.1 * rtdata.fNrays;
  size_t icrt   = 0;
  // fprintf(stderr, "P3\n%d %d\n255\n", fSize_px, fSize_py);
  for (int py = 0; py < rtdata.fSize_py; py++) {
    for (int px = 0; px < rtdata.fSize_px; px++) {
      if ((icrt % n10) == 0) printf("%lu %%\n", 10 * icrt / n10);
      int pixel_index = 4 * (py * rtdata.fSize_px + px);

      auto pixel_color = RaytraceOne(px, py, rtdata, input_buffer);

      buffer[pixel_index + 0] = pixel_color.fComp.red;
      buffer[pixel_index + 1] = pixel_color.fComp.green;
      buffer[pixel_index + 2] = pixel_color.fComp.blue;
      buffer[pixel_index + 3] = 255;
      icrt++;
    }
  }
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
*/

} // End namespace Raytracer
} // End namespace VECGEOM_IMPL_NAMESPACE

#ifndef VECGEOM_CUDA_INTERFACE
void write_ppm(std::string filename, unsigned char *buffer, int px, int py)
{
  std::ofstream image(filename);

  image << "P3\n" << px << " " << py << "\n255\n";

  for (int j = py - 1; j >= 0; j--) {
    for (int i = 0; i < px; i++) {
      int idx = 4 * (j * px + i);
      image << (int)buffer[idx + 0] << " " << (int)buffer[idx + 1] << " " << (int)buffer[idx + 2] << "\n";
    }
  }
}
#endif

} // End namespace vecgeom
