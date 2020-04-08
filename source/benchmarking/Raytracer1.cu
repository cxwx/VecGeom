/// \file Raytracer.cu
/// \author Guilherme Amadio

#include "VecGeom/benchmarking/Raytracer.h"

#include <VecGeom/base/Transformation3D.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/management/CudaManager.h>
#include <VecGeom/navigation/NavigationState.h>
#include <VecGeom/volumes/PlacedVolume.h>
#include <VecGeom/benchmarking/Raytracer.h>

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
void RenderKernel(RaytracerData_t rtdata, char *input_buffer, unsigned char *output_buffer)
{
  int px = threadIdx.x + blockIdx.x * blockDim.x;
  int py = threadIdx.y + blockIdx.y * blockDim.y;

  if ((px >= rtdata.fSize_px) || (py >= rtdata.fSize_py)) return;

  int pixel_index = 4 * (py * rtdata.fSize_px + px);

  // model view hard-coded for debugging
  // traceML size is 2200,2200,6200, centered at 0,0,0
  //printf("Creating instance of GPU ray-tracer\n");
  //rtdata.Print();

  Color_t pixel_color = Raytracer::RaytraceOne(px, py, rtdata, input_buffer);

  output_buffer[pixel_index + 0] = pixel_color.fComp.red;
  output_buffer[pixel_index + 1] = pixel_color.fComp.green;
  output_buffer[pixel_index + 2] = pixel_color.fComp.blue;
  output_buffer[pixel_index + 3] = 255;
}

void RenderGPU(cuda::VPlacedVolume const *const world, int px, int py, int maxdepth)
{
  using Vector3 = cuda::Vector3D<double>;

  // Create the raytracer model. The parameters should come as arguments from the caller
  RaytracerData_t rtdata;
  rtdata.fScreenPos.Set(0, -7000, 0);
  rtdata.fUp.Set(1, 0, 0);
  rtdata.fZoom       = 1.;
  rtdata.fModel      = kRTspecular;
  rtdata.fView       = kRTVperspective;
  rtdata.fSize_px    = px;
  rtdata.fSize_py    = py;
  rtdata.fLightColor = 0xFF0000FF;
  rtdata.fObjColor   = 0x0000FFFF;
  rtdata.fVisDepth   = 4;
  rtdata.fMaxDepth   = maxdepth;

  Raytracer::InitializeModel(world, rtdata);

  // Allocate ray data and output data on the device
  size_t statesize = NavigationState::SizeOfInstance(maxdepth);
  size_t raysize = Ray_t::SizeOfInstance(maxdepth);

  printf("=== Allocating %.3f MB of ray data on the device\n", (float)rtdata.fNrays * raysize / 1048576);
  char *input_buffer = nullptr;  
  checkCudaErrors(cudaMallocManaged((void **)&input_buffer, statesize + rtdata.fNrays * raysize));

  unsigned char *output_buffer = nullptr;
  checkCudaErrors(cudaMallocManaged((void **)&output_buffer, 4 * sizeof(unsigned char) * px * py));

  // Load and synchronize the geometry on the GPU
  vecgeom::cxx::CudaManager::Instance().LoadGeometry((vecgeom::cxx::VPlacedVolume*) world);
  vecgeom::cxx::CudaManager::Instance().Synchronize();
  
  auto gpu_world = vecgeom::cxx::CudaManager::Instance().world_gpu();
  assert(gpu_world && "GPU world volume is a null pointer");
  rtdata.fWorld   = gpu_world;

  // Initialize the navigation state for the view point
  auto vpstate = NavigationState::MakeInstanceAt(rtdata.fMaxDepth, (void *)(input_buffer));
  Raytracer::LocateGlobalPoint(rtdata.fWorld, rtdata.fStart, *vpstate, true)
  rtdata.fVPstate = vpstate;

  rtdata.Print();

  dim3 blocks(px / 8 + 1, py / 8 + 1), threads(8, 8);
  RenderKernel<<<blocks, threads>>>(gpu_world, rtdata, input_buffer, output_buffer);

  checkCudaErrors(cudaGetLastError());
  checkCudaErrors(cudaDeviceSynchronize());

  write_ppm("output.ppm", buffer, px, py);

  checkCudaErrors(cudaFree(buffer));
}
