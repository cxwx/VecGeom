#include "ShapeTester.h"
#include "VUSolid.hh"
#include "UPolycone.hh"
#include "volumes/Polycone.h"

#include "base/Vector3D.h"

#ifdef VECGEOM_ROOT
#include "TApplication.h"
#endif
#include "stdlib.h"

typedef vecgeom::SimplePolycone Poly_t;

int main(int argc, char *argv[])
{

  double Z_ValP[15];
  Z_ValP[0]  = -1520;
  Z_ValP[1]  = -804;
  Z_ValP[2]  = -804;
  Z_ValP[3]  = -515.345;
  Z_ValP[4]  = -515.345;
  Z_ValP[5]  = -177;
  Z_ValP[6]  = -177;
  Z_ValP[7]  = 149.561;
  Z_ValP[8]  = 149.561;
  Z_ValP[9]  = 575;
  Z_ValP[10] = 575;
  Z_ValP[11] = 982.812;
  Z_ValP[12] = 982.812;
  Z_ValP[13] = 1166.7;
  Z_ValP[14] = 1524;
  double R_MinP[15];
  R_MinP[0]  = 1238;
  R_MinP[1]  = 1238;
  R_MinP[2]  = 1238;
  R_MinP[3]  = 1238;
  R_MinP[4]  = 1238;
  R_MinP[5]  = 1238;
  R_MinP[6]  = 1238;
  R_MinP[7]  = 1238;
  R_MinP[8]  = 1238;
  R_MinP[9]  = 1238;
  R_MinP[10] = 1238;
  R_MinP[11] = 1238;
  R_MinP[12] = 1238;
  R_MinP[13] = 1238;
  R_MinP[14] = 1455.22;
  double R_MaxP[15];

  R_MaxP[0]  = 1555.01;
  R_MaxP[1]  = 1555.01;
  R_MaxP[2]  = 1538.05;
  R_MaxP[3]  = 1538.05;
  R_MaxP[4]  = 1523.26;
  R_MaxP[5]  = 1523.26;
  R_MaxP[6]  = 1506.24;
  R_MaxP[7]  = 1506.24;
  R_MaxP[8]  = 1488.52;
  R_MaxP[9]  = 1488.52;
  R_MaxP[10] = 1471.28;
  R_MaxP[11] = 1471.28;
  R_MaxP[12] = 1455.22;
  R_MaxP[13] = 1455.22;
  R_MaxP[14] = 1455.22;

  int Nz         = 15;
  VUSolid *poly2 = new Poly_t("Test", 0.,              /* initial phi starting angle */
                              UUtils::kTwoPi,          /* total phi angle */
                              Nz,                      /* number corners in r,z space */
                              Z_ValP, R_MinP, R_MaxP); /* r coordinate of these corners */

  ShapeTester tester;

  if (argc > 1) {
    if (strcmp(argv[1], "vis") == 0) {
#ifdef VECGEOM_ROOT
      TApplication theApp("App", 0, 0);
      tester.Run(poly2);
      theApp.Run();
#endif
    }
  } else {
    tester.Run(poly2);
  }

  delete poly2;

  return 0;
}
