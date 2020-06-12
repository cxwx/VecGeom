/*
 * RaytraceBenchmark.cpp
 *
 *  Created on: Feb 13, 2020
 *      Author: andrei.gheata@cern.ch
 */

#include <VecGeom/base/Vector3D.h>
#include <VecGeom/benchmarking/Raytracer.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/navigation/NavigationState.h>
#include "ArgParser.h"

#ifdef VECGEOM_GDML
#include "Frontend.h"
#endif

using namespace vecgeom;

int RaytraceBenchmarkCPU(vecgeom::cxx::RaytracerData_t &rtdata);
#ifdef VECGEOM_ENABLE_CUDA
int RaytraceBenchmarkGPU(vecgeom::cuda::VPlacedVolume const* const world,
                         int px, int py, int model, int view, double zoom,
			 double screenx, double screeny, double screenz,
			 double upx, double upy, double upz,
			 int lightcol, int objcol, int maxdepth, int vdepth);
#endif

int main(int argc, char *argv[])
{
#ifndef VECGEOM_GDML
  std::cout << "### VecGeom must be compiled with GDML support to run this.\n";
  return 1;
#endif

  OPTION_STRING(gdml_name, "default.gdml");
  // image size in pixels
  OPTION_INT(px, 1200);
  OPTION_INT(py, 800);

  // RT model as in { kRTxray = 0, kRTspecular, kRTtransparent, kRTdiffuse };
  OPTION_INT(model, 2);

  // RT view as in { kRTVparallel = 0, kRTVperspective };
  OPTION_INT(view, 1);

  // zoom w.r.t to the default view mode
  OPTION_DOUBLE(zoom, 1.5);

  // Screen position in world coordinates
  OPTION_DOUBLE(screenx, -10000);
  OPTION_DOUBLE(screeny, -5000);
  OPTION_DOUBLE(screenz, 0);

  // Up vector (no need to be normalized)
  OPTION_DOUBLE(upx, 0);
  OPTION_DOUBLE(upy, 1);
  OPTION_DOUBLE(upz, 0);
  vecgeom::Vector3D<double> up(upx, upy, upz);

  // Light color, object color (no color per volume yet) - in RGBA chars compressed into an unsigned integer
  OPTION_INT(lightcol, 0xFF0000FF); // red
  OPTION_INT(objcol, 0x0000FFFF);   // blue
  OPTION_INT(vdepth, 4);            // visible depth

  OPTION_INT(on_gpu, 0);            // run on GPU

// Try to open the input file
#ifdef VECGEOM_GDML
  bool load = vgdml::Frontend::Load(gdml_name.c_str());
  if (!load) return 2;
#endif

  auto world = GeoManager::Instance().GetWorld();
  if (!world) return 3;
  RaytracerData_t rtdata;

  rtdata.fScreenPos.Set(screenx, screeny, screenz);
  rtdata.fUp.Set(upx, upy, upz);
  rtdata.fZoom       = zoom;
  rtdata.fModel      = (ERTmodel)model;
  rtdata.fView       = (ERTView)view;
  rtdata.fSize_px    = px;
  rtdata.fSize_py    = py;
  rtdata.fLightColor = lightcol;
  rtdata.fObjColor   = objcol;
  rtdata.fVisDepth   = vdepth;
  rtdata.fMaxDepth   = GeoManager::Instance().getMaxDepth();

  Raytracer::InitializeModel(world, rtdata);
  rtdata.Print();

  auto ierr = 0;
  if (on_gpu) {
#ifdef VECGEOM_ENABLE_CUDA
    ierr = RaytraceBenchmarkGPU((const vecgeom::cuda::VPlacedVolume *)world, px, py, model, view, zoom,
                                screenx, screeny, screenz, upx, upy, upz,
                                lightcol, objcol, rtdata.fMaxDepth, vdepth);
#else
    std::cout << "=== Cannot run RaytracerBenchmark on GPU since VecGeom CUDA support not compiled.\n";
    return 1;
#endif
  } else {
   ierr = RaytraceBenchmarkCPU(rtdata);
  }
  if (ierr) std::cout << "TestNavIndex FAILED\n";

  return ierr;
}

int RaytraceBenchmarkCPU(vecgeom::cxx::RaytracerData_t &rtdata)
{
  // Allocate and initialize all rays on the host
  size_t statesize = NavigationState::SizeOfInstance(rtdata.fMaxDepth);
  size_t raysize = Ray_t::SizeOfInstance(rtdata.fMaxDepth);
  printf("=== Allocating %.3f MB of ray data on the host\n", (float)rtdata.fNrays * raysize / 1048576);
  char *input_buffer = new char[statesize + rtdata.fNrays * raysize];
  char *output_buffer = new char[4 * rtdata.fNrays * sizeof(char)];

  // Initialize the navigation state for the view point
  auto vpstate = NavigationState::MakeInstanceAt(rtdata.fMaxDepth, (void *)(input_buffer));
  Raytracer::LocateGlobalPoint(rtdata.fWorld, rtdata.fStart, *vpstate, true);

  rtdata.fVPstate = vpstate;

  // Construct rays in place
  char *raybuff = input_buffer + statesize;
  for (int iray = 0; iray < rtdata.fNrays; ++iray)
    Ray_t::MakeInstanceAt(raybuff + iray * raysize, rtdata.fMaxDepth);

  // Run the CPU propagation kernel
  Raytracer::PropagateRays(rtdata, input_buffer, output_buffer);

  // Write the output
  write_ppm("output.ppm", (unsigned char*)output_buffer, rtdata.fSize_px, rtdata.fSize_py);

  return 0;
}
