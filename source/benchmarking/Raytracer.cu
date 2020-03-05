/// \file Raytracer.cu
/// \author Guilherme Amadio

#include "VecGeom/benchmarking/Raytracer.h"

#include <VecGeom/base/Transformation3D.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/management/CudaManager.h>
#include <VecGeom/navigation/NavigationState.h>
#include <VecGeom/volumes/PlacedVolume.h>

#include <cassert>
#include <cstdio>

using namespace vecgeom;

void check_cuda(cudaError_t result, char const *const func, const char *const file, int const line)
{
  if (result) {
    fprintf(stderr, "CUDA error = %s at %s:%d\n", cudaGetErrorString(result), file, line);
    cudaDeviceReset();
    exit(1);
  }
}

#define checkCudaErrors(val) check_cuda((val), #val, __FILE__, __LINE__)

__device__
VPlacedVolume const*
LocateGlobalPoint(VPlacedVolume const *vol, Vector3D<Precision> const &point, NavigationState &path, bool top)
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

__device__
VPlacedVolume const*
LocateGlobalPointExclVolume(VPlacedVolume const *vol, VPlacedVolume const *excludedvolume,
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

__device__
VPlacedVolume const*
RelocatePointFromPathForceDifferent(Vector3D<Precision> const &localpoint, NavigationState &path)
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

__device__
Precision
ComputeStepAndPropagatedState(Vector3D<Precision> const &globalpoint, Vector3D<Precision> const &globaldir,
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

__device__
Color_t raytrace(VPlacedVolume const* const world, Vector3D<Precision> origin, Vector3D<Precision> dir, int maxdepth)
{
  Ray_t ray;
  char cstate[512];
  char nstate[512];

  ray.fPos = origin;
  ray.fDir = dir;
  ray.fColor = 0xFFFFFFFF;
  ray.fCrtState  = NavigationState::MakeInstanceAt(maxdepth, (void *)(cstate));
  ray.fNextState = NavigationState::MakeInstanceAt(maxdepth, (void *)(nstate));

  ray.fVolume = LocateGlobalPoint(world, ray.fPos, *ray.fCrtState, true);

  do {
    auto t = world->DistanceToIn(ray.fPos, ray.fDir);

    if (t == kInfLength)
      return ray.fColor;

    ray.fPos += (t + 1e-7) * ray.fDir;
    ray.fVolume = LocateGlobalPoint(world, ray.fPos, *ray.fCrtState, true);

  } while (!ray.fVolume);

  if (!ray.fVolume)
    return ray.fColor;

  *ray.fNextState = *ray.fCrtState;

  while (ray.fVolume) {
    auto t = ComputeStepAndPropagatedState(ray.fPos, ray.fDir, kInfLength, *ray.fCrtState, *ray.fNextState);

    if (t == kInfLength || !ray.fNextState->Top())
      break;

    ray.fPos += (t + 1e-7) * ray.fDir;
    ray.fVolume = ray.fNextState->Top();
    *ray.fCrtState = *ray.fNextState;

    if (ray.fVolume->GetDaughters().size() == 0)
      ray.fColor += 0x0000ff15;
  }

  return ray.fColor;
}

#if 0
void RenderCPU(VPlacedVolume const *const world, int px, int py, int maxdepth)
{
  float *buffer = (float*) malloc(4 * px * py * sizeof(float));

  for (int i = 0; i < px; ++i) {
    for (int j = 0; j < py; ++j) {
      int pixel_index = 4 * (j * px + i);

      float u = float(i) / float(px);
      float v = float(j) / float(py);

      // model view hard-coded for debugging
      // traceML size is 2200,2200,6200, centered at 0,0,0
      Vector3D<Precision> origin = {0, -7000, 0};
      Vector3D<Precision> direction = {v - 0.5, 1.9, 2*u - 1};
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
#endif

__global__
void RenderKernel(cuda::VPlacedVolume const *const world, int px, int py, float *buffer, int maxdepth)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;

  if ((i >= px) || (j >= py)) return;

  int pixel_index = 4 * (j * px + i);

  float u = float(i) / float(px);
  float v = float(j) / float(py);

  // model view hard-coded for debugging
  // traceML size is 2200,2200,6200, centered at 0,0,0
  Vector3D<Precision> origin = {0, -7000, 0};
  Vector3D<Precision> direction = {v - 0.5, 1.9, 2*u - 1};
  direction.Normalize();

  Color_t pixel_color = raytrace(world, origin, direction, maxdepth);

  buffer[pixel_index + 0] = pixel_color.Red();
  buffer[pixel_index + 1] = pixel_color.Green();
  buffer[pixel_index + 2] = pixel_color.Blue();
  buffer[pixel_index + 3] = 1.0f;
}

void RenderGPU(cuda::VPlacedVolume const *const world, int px, int py, int maxdepth)
{
  float *buffer = nullptr;
  checkCudaErrors(cudaMallocManaged((void **)&buffer, 4 * sizeof(float) * px * py));

  vecgeom::cxx::CudaManager::Instance().LoadGeometry((vecgeom::cxx::VPlacedVolume*) world);
  vecgeom::cxx::CudaManager::Instance().Synchronize();

  auto gpu_world = vecgeom::cxx::CudaManager::Instance().world_gpu();

  if (!gpu_world) {
    fprintf(stderr, "GPU world volume is a null pointer\n");
    exit(1);
  }

  dim3 blocks(px / 8 + 1, py / 8 + 1), threads(8, 8);
  RenderKernel<<<blocks, threads>>>(gpu_world, px, py, buffer, maxdepth);

  checkCudaErrors(cudaGetLastError());
  checkCudaErrors(cudaDeviceSynchronize());

  write_ppm("output.ppm", buffer, px, py);

  checkCudaErrors(cudaFree(buffer));
}
