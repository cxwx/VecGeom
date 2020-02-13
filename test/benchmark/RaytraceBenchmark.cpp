/*
 * RaytraceBenchmark.cpp
 *
 *  Created on: Feb 13, 2020
 *      Author: swenzel
 */

#include <VecGeom/base/Vector3D.h>
#include <VecGeom/benchmarking/Raytracer.h>
#include <VecGeom/management/GeoManager.h>
#include "ArgParser.h"

#ifdef VECGEOM_GDML
#include "Frontend.h"
#endif

using namespace vecgeom;

int main(int argc, char *argv[])
{
#ifndef VECGEOM_GDML
  std::cout << "### VecGeom must be compiled with GDML support to run this.\n";
  return 1;
#endif

  OPTION_STRING(gdml_name, "default.gdml");
  // image size in pixels
  OPTION_INT(px, 1024);
  OPTION_INT(py, 1024);

  // RT model as in { kRTxray = 0, kRTspecular, kRTtransparent, kRTdiffuse };
  OPTION_INT(model, 1);

  // Light source position in world coordinates
  OPTION_DOUBLE(srcx, -9999);
  OPTION_DOUBLE(srcy, -9999);
  OPTION_DOUBLE(srcz, -9999);
  vecgeom::Vector3D<double> src_pos(srcx, srcy, srcz);

  // Up vector (no need to be normalized)
  OPTION_DOUBLE(upx, 0);
  OPTION_DOUBLE(upy, 0);
  OPTION_DOUBLE(upz, 1);
  vecgeom::Vector3D<double> up(upx, upy, upz);

  // Light color, object color (no color per volume yet) - in RGBA chars compressed into an unsigned integer
  OPTION_INT(lightcol, 0xFFFFFFFF); // white
  OPTION_INT(objcol, 0x0000FFFF);   // blue

// Try to open the input file
#ifdef VECGEOM_GDML
  bool load = vgdml::Frontend::Load(gdml_name.c_str());
  if (!load) return 2;
#endif

  auto world = vecgeom::GeoManager::Instance().GetWorld();
  if (!world) return 3;
  vecgeom::Raytracer raytracer(world, src_pos, up, px, py, (vecgeom::ERTmodel)model);
  raytracer.SetLightColor(lightcol);
  raytracer.SetObjColor(objcol);

  raytracer.PropagateRays();
  return 0;
}
