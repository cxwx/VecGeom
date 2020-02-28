/// \file Raytracer.cu
/// \author Guilherme Amadio

#include "VecGeom/benchmarking/Raytracer.h"

#include <VecGeom/base/Transformation3D.h>
#include <VecGeom/base/Stopwatch.h>

#include <VecGeom/navigation/NavigationState.h>
#include <VecGeom/volumes/PlacedVolume.h>
#include <VecGeom/management/GeoManager.h>

#ifdef VECGEOM_CUDA_INTERFACE
#include "VecGeom/backend/cuda/Backend.h"
#include "VecGeom/management/CudaManager.h"
#endif

using namespace vecgeom;

#define checkCudaErrors(val) check_cuda((val), #val, __FILE__, __LINE__)

void check_cuda(cudaError_t result, char const *const func, const char *const file, int const line)
{
  if (result) {
    fprintf(stderr, "CUDA error = %s at %s:%d\n", cudaGetErrorString(result), file, line);
    cudaDeviceReset();
    exit(1);
  }
}

__device__
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

__device__
VPlacedVolume const *LocateGlobalPointExclVolume(VPlacedVolume const *vol,
                                                            VPlacedVolume const *excludedvolume,
                                                            Vector3D<Precision> const &point, NavigationState &path,
                                                            bool top)
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
VPlacedVolume const *RelocatePointFromPathForceDifferent(Vector3D<Precision> const &localpoint,
                                                                    NavigationState &path)
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
double ComputeStepAndPropagatedState(Vector3D<Precision> const &globalpoint,
                                                Vector3D<Precision> const &globaldir, Precision step_limit,
                                                NavigationState const &in_state, NavigationState &out_state)
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

__device__
void ApplyRTmodel(Ray_t &ray, double step)
{
	int fVisDepth       = 4;               ///< Visible geometry depth
	Color_t fLightColor = 0xFF0000FF;      ///< Light color
	Color_t fObjColor   = 0x0000FFFF;      ///< Object color
	ERTmodel fModel     = kRTspecular;  ///< Selected RT model
	Vector3D<double> LightDir(1,1,1);   ///< Light source direction
	LightDir.Normalize();

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
      Vector3D<double> refl = LightDir - 2 * norm.Dot(LightDir) * norm;
      refl.Normalize();
      double calf = -ray.fDir.Dot(refl);
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

__device__
Color_t trace_ray(VPlacedVolume const* const world, Vector3D<double> origin, Vector3D<double> dir)
{
  constexpr int kMaxTries      = 10;
  constexpr double kPush       = 1.e-8;

  int maxdepth = 8;
  // size_t navstate_size = NavigationState::SizeOfInstance(maxdepth);

  constexpr size_t N = 512; // should contain current/next states
  char cstate[N];
  char nstate[N];

  Ray_t ray;

  ray.fPos = origin;
  ray.fDir = dir;
  ray.fColor  = 0xffffffff; //0xAA33FFFF; // purple
  ray.fCrtState  = NavigationState::MakeInstanceAt(maxdepth, (void *)(cstate));
  ray.fNextState = NavigationState::MakeInstanceAt(maxdepth, (void *)(nstate));
  ray.fVolume = LocateGlobalPoint(world, ray.fPos, *ray.fCrtState, true);

  int itry    = 0;
  while (!ray.fVolume && itry < kMaxTries) {
    auto snext = world->DistanceToIn(ray.fPos, ray.fDir);
    ray.fDone  = snext == kInfLength;
    if (ray.fDone) return ray.fColor;
    // Propagate to the world volume (but do not increment the boundary count)
    ray.fPos += (snext + kPush) * ray.fDir;
    ray.fVolume = LocateGlobalPoint(world, ray.fPos, *ray.fCrtState, true);
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
      //std::cout << "error for ray (" << px << ", " << py << ")\n";
      ray.fDone  = true;
      ray.fColor = 0;
      return ray.fColor;
    }

    // Apply the selected RT model
    ray.fNcrossed++;
    ray.fVolume = nextvol;
    if (ray.fVolume == nullptr) ray.fDone = true;
    if (nextvol) ApplyRTmodel(ray, snext);
    auto tmpstate  = ray.fCrtState;
    ray.fCrtState  = ray.fNextState;
    ray.fNextState = tmpstate;
  }

  return ray.fColor;
}

__global__ void Render(cuda::VPlacedVolume const *const world, int px, int py, float *buffer)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;

  if ((i >= px) || (j >= py)) return;

  int pixel_index = 4 * (j * px + i);

  float u = float(i) / float(px);
  float v = float(j) / float(py);

  // hard-coded orthographic for debugging
  // traceML size is 2200,2200,6200, centered at 0,0,0

  Vector3D<double> origin = {2199 * (u-0.5), 2199*(v-0.5), -3201}; 
  Vector3D<double> direction = {0, 0, 1};

  Color_t pixel_color = trace_ray(world, origin, direction);

  buffer[pixel_index + 0] = pixel_color.Red();
  buffer[pixel_index + 1] = pixel_color.Green();
  buffer[pixel_index + 2] = pixel_color.Blue();
  buffer[pixel_index + 3] = 1.0f;
}

void RayTraceCUDA(cxx::VPlacedVolume const* const world, int px, int py,
                  float a, float b, float c, float d, float e, float f)
{
	float *buffer = nullptr;
	checkCudaErrors(cudaMallocManaged((void **)&buffer, 4 * sizeof(float) * px * py));

	cxx::CudaManager::Instance().LoadGeometry(world);
	cxx::CudaManager::Instance().Synchronize();

	auto gpu_world = cxx::CudaManager::Instance().world_gpu();

	if (!gpu_world) {
		fprintf(stderr, "GPU world volume is a null pointer\n");
		exit(1);
	}

	dim3 blocks(px / 8 + 1, py / 8 + 1), threads(8, 8);
	Render<<<blocks, threads>>>(gpu_world, px, py, buffer);

	std::ofstream image("out.ppm");

	image << "P3\n" << px << " " << py << "\n255\n";
	for (int j = py - 1; j >= 0; j--) {
		for (int i = 0; i < px; i++) {
			size_t pixel_index = 4 * (j * px + i);
			int ir             = int(255.99 * buffer[pixel_index + 0]);
			int ig             = int(255.99 * buffer[pixel_index + 1]);
			int ib             = int(255.99 * buffer[pixel_index + 2]);
			image << ir << " " << ig << " " << ib << "\n";
		}
	}

	checkCudaErrors(cudaFree(buffer));
}
