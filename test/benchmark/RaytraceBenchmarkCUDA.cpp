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

void RayTraceCUDA(VPlacedVolume const* const world, int px, int py,
                  float a, float b, float c, float d, float e, float f);

int main(int argc, char *argv[])
{
#if !defined(VECGEOM_GDML) || !defined(VECGEOM_ENABLE_CUDA)
   fprintf(stderr, "VecGeom must have GDML and CUDA enabled to run this example.\n");
   return 1;
#else
   int nx = 1024;
   int ny = 1024;

   if (!vgdml::Frontend::Load(argv[1])) {
      fprintf(stderr, "Failed to open %s\n", argv[1]);
      return 2;
   }

   auto world = GeoManager::Instance().GetWorld();

   Vector3D<double> aMin, aMax;
	world->Extent(aMin, aMax);

   if (!world) {
      fprintf(stderr, "World volume is a null pointer\n");
      return 3;
   }

   Vector3D<double> wmin, wmax, wcenter, wsize;

   world->Extent(wmin, wmax);
   wcenter = 0.5 * (wmax + wmin);
   wsize   =       (wmax - wmin);

   std::cout << "world min: " << wmin << ", max: " << wmax << "\n";
   std::cout << "world center: " << wcenter << ", size: " << wsize << "\n";

   RayTraceCUDA(world, nx, ny, aMin[0], aMin[1], aMin[2],
                               aMax[0], aMax[1], aMax[2]);

   return 0;
#endif
}
