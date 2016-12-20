#include "../benchmark/ArgParser.h"
#include "ShapeTester.h"
#include "volumes/Orb.h"

typedef vecgeom::SimpleOrb Orb_t;

int main(int argc, char *argv[])
{
  OPTION_INT(npoints, 10000);
  OPTION_BOOL(debug, false);
  OPTION_BOOL(stat, false);
  OPTION_BOOL(usolids, false);

  if (usolids) {
    std::cerr << "\n*** ERROR: '-usolids true' is not valid for SExtru shape!\n Aborting...\n\n";
    return 1;
  }

  Orb_t const *orb = new Orb_t("test_orb", 35);
  orb->Print();

  ShapeTester<vecgeom::VPlacedVolume> tester;
  tester.setConventionsMode(usolids);
  tester.setDebug(debug);
  tester.setStat(stat);
  tester.SetMaxPoints(npoints);
  tester.SetSolidTolerance(1.e-9);
  tester.SetTestBoundaryErrors(true);
  int errCode = tester.Run(orb);

  std::cout << "Final Error count for Shape *** " << orb->GetName() << "*** = " << errCode << " ("
            << (tester.getConventionsMode() ? "USolids" : "VecGeom") << " conventions)\n";
  std::cout << "=========================================================" << std::endl;

  if (orb) delete orb;
  return 0;
}
