/// \file shape_testParallelepiped.h
/// \author: Mihaela Gheata (mihaela.gheata@cern.ch)

#include "../benchmark/ArgParser.h"
#include "ShapeTester.h"
#include "volumes/Parallelepiped.h"
typedef vecgeom::SimpleParallelepiped Para_t;

int main(int argc, char *argv[])
{
  if (argc == 1) {
    std::cout << "Usage: shape_testParallelepiped -test <#>:\n"
                 "       0 - alpha, theta, phi\n"
                 "       1 - alpha=0 theta=0 phi=0 (box)\n";
  }

  OPTION_INT(npoints, 10000);
  OPTION_BOOL(debug, false);
  OPTION_BOOL(stat, false);
  OPTION_BOOL(usolids, false);
  OPTION_INT(type, 0);

  if (usolids) {
    std::cerr << "\n*** ERROR: '-usolids true' is not valid for SExtru shape!\n Aborting...\n\n";
    return 1;
  }

  using namespace vecgeom;

  Para_t *solid   = 0;
  Precision dx    = 10.;
  Precision dy    = 7;
  Precision dz    = 15;
  Precision alpha = 30.;
  Precision theta = 30.;
  Precision phi   = 45.;

  switch (type) {
  case 0:
    // Non-zero alpha, theta, phi
    std::cout << "Testing general parallelepiped\n";
    solid = new Para_t("test_VecGeomPara", dx, dy, dz, alpha, theta, phi);
    break;
  case 1:
    // Parallelepiped degenerated to box
    std::cout << "Testing box-like parallelepiped\n";
    solid = new Para_t("test_VecGeomPara", dx, dy, dz, 0, 0, 0);
    break;
  default:
    std::cout << "Unknown test case.\n";
  }

  ShapeTester<vecgeom::VPlacedVolume> tester;
  tester.setConventionsMode(usolids);
  tester.setDebug(debug);
  tester.setStat(stat);
  tester.SetMaxPoints(npoints);
  tester.SetSolidTolerance(1.e-9);
  tester.SetTestBoundaryErrors(true);
  int errCode = tester.Run(solid);

  std::cout << "Final Error count for Shape *** " << solid->GetName() << "*** = " << errCode << " ("
            << (tester.getConventionsMode() ? "USolids" : "VecGeom") << " conventions)\n";
  std::cout << "=========================================================" << std::endl;

  if (solid) delete solid;
  return 0;
}
