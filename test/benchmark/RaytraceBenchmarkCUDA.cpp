#include <VecGeom/base/Vector3D.h>
#include <VecGeom/navigation/GlobalLocator.h>
#include <VecGeom/navigation/NavigationState.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/management/CudaManager.h>
#include <VecGeom/benchmarking/Raytracer.h>
#include <VecGeom/management/GeoManager.h>

#ifdef VECGEOM_GDML
#include "Frontend.h"
#endif

using namespace vecgeom;

#ifdef VECGEOM_ENABLE_CUDA
void RayTrace(vecgeom::cuda::VPlacedVolume const* const world, int px, int py, int maxdepth, bool use_cuda);
#endif

int main(int argc, char *argv[])
{
#if !defined(VECGEOM_GDML) || !defined(VECGEOM_ENABLE_CUDA)
   fprintf(stderr, "VecGeom must have GDML and CUDA enabled to run this example.\n");
   return 1;
#else
   int nx = 2048;
   int ny = 1024;

   if (!vgdml::Frontend::Load(argv[1])) {
      fprintf(stderr, "Failed to open %s\n", argv[1]);
      return 2;
   }

   auto world = GeoManager::Instance().GetWorld();

   if (!world) {
      fprintf(stderr, "World volume is a null pointer\n");
      return 3;
   }

   bool use_cuda = true;
   int maxdepth = GeoManager::Instance().getMaxDepth();

   if (maxdepth <= 0)
      fprintf(stderr, "Error, geometry depth is reported as %d\n", maxdepth);
   else
      fprintf(stdout, "geometry depth is %d\n", maxdepth);

   RayTrace((vecgeom::cuda::VPlacedVolume*)world, nx, ny, maxdepth, use_cuda);

   return 0;
#endif
}
