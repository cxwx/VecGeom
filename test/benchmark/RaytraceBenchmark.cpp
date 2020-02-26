/*
 * RaytraceBenchmark.cpp
 *
 *  Created on: Feb 13, 2020
 *      Author: andrei.gheata@cern.ch
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

  // RT view as in { kRTVparallel = 0, kRTVperspective };
  OPTION_INT(view, 1);

  // Screen position in world coordinates
  OPTION_DOUBLE(screenx, -2000);
  OPTION_DOUBLE(screeny, -1500);
  OPTION_DOUBLE(screenz, -600);
  vecgeom::Vector3D<double> screen_pos(screenx, screeny, screenz);

  // Up vector (no need to be normalized)
  OPTION_DOUBLE(upx, 0);
  OPTION_DOUBLE(upy, 1);
  OPTION_DOUBLE(upz, 0);
  vecgeom::Vector3D<double> up(upx, upy, upz);

  // Light source direction
  OPTION_DOUBLE(ldirx, 0);
  OPTION_DOUBLE(ldiry, 0);
  OPTION_DOUBLE(ldirz, 0);
  vecgeom::Vector3D<double> ldir(ldirx, ldiry, ldirz);

  // Light color, object color (no color per volume yet) - in RGBA chars compressed into an unsigned integer
  OPTION_INT(lightcol, 0xFF0000FF); // red
  OPTION_INT(objcol, 0x0000FFFF);   // blue
  OPTION_INT(vdepth, 1);            // visible depth

// Try to open the input file
#ifdef VECGEOM_GDML
  bool load = vgdml::Frontend::Load(gdml_name.c_str());
  if (!load) return 2;
#endif

  auto world = vecgeom::GeoManager::Instance().GetWorld();
  if (!world) return 3;
  vecgeom::Raytracer raytracer(world, screen_pos, up, px, py, (vecgeom::ERTmodel)model, (vecgeom::ERTView)view);
  raytracer.SetLightColor(lightcol);
  raytracer.SetObjColor(objcol);
  // raytracer.SetLightSourceDir(ldir);
  raytracer.SetVisDepth(vdepth);

  raytracer.PropagateRays();
  return 0;
}
