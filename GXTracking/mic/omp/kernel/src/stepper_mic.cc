#include "GPTypeDef.h"
#include "GXTrack.h"
#include "GXFieldMap.h"
#include "GXMagneticField.h"
#include "GXEquationOfMotion.h"
#include "GXClassicalRK4.h"
#include "GXCashKarpRKF45.h"
#include "GXNystromRK4.h"

#include <omp.h>
#include "stdio.h"
//-----------------------------------------------------------------------------
//  MIC wrapper for the coprocessor
//-----------------------------------------------------------------------------

GLOBALFUNC
void rk4_mic(GXFieldMap *magMap, GXTrack* track, unsigned int ntracks)
{
  GXMagneticField magField;    
  GXMagneticField_Constructor(&magField,magMap);

  GXEquationOfMotion equaOfMotion;
  GXEquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);

  GXClassicalRK4 rk4;
  GXClassicalRK4_Constructor(&rk4,&equaOfMotion);

  const G4int nvar = 6;
  G4double y[nvar], dydx[nvar], yOut[nvar], yErr[nvar];
  G4double hstep,terror,dchord;

  unsigned int tid;

  //  int nthreads = omp_get_max_threads();
  //  int nthreads = omp_get_num_threads();
  //  #pragma omp parallel for num_threads(nthreads)

  //  for (tid = 0; tid < ntracks ; tid++) {

  y[0]  = track->x;
  y[1]  = track->y;
  y[2]  = track->z;
  y[3]  = track->px;
  y[4]  = track->py;
  y[5]  = track->pz;
  hstep = track->s;
  
  GXClassicalRK4_RightHandSide(&rk4,y,dydx);
  GXClassicalRK4_Stepper(&rk4,y,dydx,hstep,yOut,yErr);
  
  terror = GXClassicalRK4_TruncatedError(&rk4,hstep,yOut,yErr);
  dchord = GXClassicalRK4_DistChord(&rk4);
  
  //output
  track->x    = yOut[0];
  track->y    = yOut[1];
  track->z    = yOut[2];
  track->px   = yOut[3];
  track->py   = yOut[4];
  track->pz   = yOut[5];
  track->s    = terror+dchord; //temporary
  //  }
}

//---------------------------------------------
// CashKarpRKF45 
//---------------------------------------------

GLOBALFUNC
void rkf45_mic( GXFieldMap *magMap, GXTrack *track, unsigned int nTrackSize)
{
  GXMagneticField magField;    
  GXMagneticField_Constructor(&magField,magMap);

  GXEquationOfMotion equaOfMotion;
  GXEquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);

  GXCashKarpRKF45 rkf45;
  GXCashKarpRKF45_Constructor(&rkf45,&equaOfMotion);

  const G4int nvar = 6;
  G4double y[nvar], dydx[nvar], yOut[nvar], yErr[nvar];
  G4double hstep,terror,dchord;

  //  for (int tid = 0 ; tid < nTrackSize ; tid++) {

  //input
  y[0]  = track->x;
  y[1]  = track->y;
  y[2]  = track->z;
  y[3]  = track->px;
  y[4]  = track->py;
  y[5]  = track->pz;
  hstep = track->s;
  
  GXCashKarpRKF45_RightHandSide(&rkf45,y,dydx);
  GXCashKarpRKF45_Stepper(&rkf45,y,dydx,hstep,yOut,yErr);
  
  terror = GXCashKarpRKF45_TruncatedError(&rkf45,hstep,yOut,yErr);
  dchord = GXCashKarpRKF45_DistChord(&rkf45);
  
  //output
  track->x    = yOut[0];
  track->y    = yOut[1];
  track->z    = yOut[2];
  track->px   = yOut[3];
  track->py   = yOut[4];
  track->pz   = yOut[5];
  track->s    = terror + dchord; //temporary
}

//---------------------------------------------
// NystromRK4
//---------------------------------------------

GLOBALFUNC
void nrk4_mic( GXFieldMap *magMap, GXTrack *track, unsigned int nTrackSize)
{
  GXMagneticField magField;    
  GXMagneticField_Constructor(&magField,magMap);

  GXEquationOfMotion equaOfMotion;
  GXEquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);

  GXNystromRK4 nrk4;
  GXNystromRK4_Constructor(&nrk4,&equaOfMotion);

  const G4int nvar = 6;
  G4double y[nvar], dydx[nvar], yOut[nvar], yErr[nvar];
  G4double hstep,terror,dchord;

  //  for (int tid = 0 ; tid < nTrackSize ; tid++) {

  //input
  y[0] = track->x;
  y[1] = track->y;
  y[2] = track->z;
  y[3] = track->px;
  y[4] = track->py;
  y[5] = track->pz;
  hstep = track->s;
  
  GXNystromRK4_RightHandSide(&nrk4,y,dydx);
  GXNystromRK4_Stepper(&nrk4,y,dydx,hstep,yOut,yErr);
  
  terror = GXNystromRK4_TruncatedError(&nrk4,hstep,yOut,yErr);
  dchord = GXNystromRK4_DistChord(&nrk4);
  
  //output
  track->x    = yOut[0];
  track->y    = yOut[1];
  track->z    = yOut[2];
  track->px   = yOut[3];
  track->py   = yOut[4];
  track->pz   = yOut[5];
  track->s    = terror + dchord; //temporary
  //  }
}

void rk4_cpu(GXFieldMap *magMap, GXTrack* track, unsigned int ntracks)
{
  GXMagneticField magField;    
  GXMagneticField_Constructor(&magField,magMap);

  GXEquationOfMotion equaOfMotion;
  GXEquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);

  GXClassicalRK4 rk4;
  GXClassicalRK4_Constructor(&rk4,&equaOfMotion);

  const G4int nvar = 6;
  G4double y[nvar], dydx[nvar], yOut[nvar], yErr[nvar];
  G4double hstep,terror,dchord;

  unsigned int tid;

  //  int nthreads = omp_get_max_threads();
  //  int nthreads = omp_get_num_threads();
  //  #pragma omp parallel for num_threads(nthreads)

  for (tid = 0; tid < ntracks ; tid++) {

    y[0] = track[tid].x;
    y[1] = track[tid].y;
    y[2] = track[tid].z;
    y[3] = track[tid].px;
    y[4] = track[tid].py;
    y[5] = track[tid].pz;
    hstep = track[tid].s;

    GXClassicalRK4_RightHandSide(&rk4,y,dydx);
    GXClassicalRK4_Stepper(&rk4,y,dydx,hstep,yOut,yErr);

    terror = GXClassicalRK4_TruncatedError(&rk4,hstep,yOut,yErr);
    dchord = GXClassicalRK4_DistChord(&rk4);

    //output
    track[tid].x    = yOut[0];
    track[tid].y    = yOut[1];
    track[tid].z    = yOut[2];
    track[tid].px   = yOut[3];
    track[tid].py   = yOut[4];
    track[tid].pz   = yOut[5];
    track[tid].s    = terror+dchord; //temporary
  }
}

//---------------------------------------------
// CashKarpRKF45 
//---------------------------------------------

void rkf45_cpu( GXFieldMap *magMap, GXTrack *track, unsigned int nTrackSize)
{
  GXMagneticField magField;    
  GXMagneticField_Constructor(&magField,magMap);

  GXEquationOfMotion equaOfMotion;
  GXEquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);

  GXCashKarpRKF45 rkf45;
  GXCashKarpRKF45_Constructor(&rkf45,&equaOfMotion);

  const G4int nvar = 6;
  G4double y[nvar], dydx[nvar], yOut[nvar], yErr[nvar];
  G4double hstep,terror,dchord;

  for (int tid = 0 ; tid < nTrackSize ; tid++) {

    //input
    y[0] = track[tid].x;
    y[1] = track[tid].y;
    y[2] = track[tid].z;
    y[3] = track[tid].px;
    y[4] = track[tid].py;
    y[5] = track[tid].pz;
    hstep = track[tid].s;

    GXCashKarpRKF45_RightHandSide(&rkf45,y,dydx);
    GXCashKarpRKF45_Stepper(&rkf45,y,dydx,hstep,yOut,yErr);

    terror = GXCashKarpRKF45_TruncatedError(&rkf45,hstep,yOut,yErr);
    dchord = GXCashKarpRKF45_DistChord(&rkf45);

    //output
    track[tid].x    = yOut[0];
    track[tid].y    = yOut[1];
    track[tid].z    = yOut[2];
    track[tid].px   = yOut[3];
    track[tid].py   = yOut[4];
    track[tid].pz   = yOut[5];
    track[tid].s    = terror + dchord; //temporary
  }
}

//---------------------------------------------
// NystromRK4
//---------------------------------------------

void nrk4_cpu( GXFieldMap *magMap, GXTrack *track, unsigned int nTrackSize)
{
  GXMagneticField magField;    
  GXMagneticField_Constructor(&magField,magMap);

  GXEquationOfMotion equaOfMotion;
  GXEquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);

  GXNystromRK4 nrk4;
  GXNystromRK4_Constructor(&nrk4,&equaOfMotion);

  const G4int nvar = 6;
  G4double y[nvar], dydx[nvar], yOut[nvar], yErr[nvar];
  G4double hstep,terror,dchord;

  for (int tid = 0 ; tid < nTrackSize ; tid++) {

    //input
    y[0] = track[tid].x;
    y[1] = track[tid].y;
    y[2] = track[tid].z;
    y[3] = track[tid].px;
    y[4] = track[tid].py;
    y[5] = track[tid].pz;
    hstep = track[tid].s;

    GXNystromRK4_RightHandSide(&nrk4,y,dydx);
    GXNystromRK4_Stepper(&nrk4,y,dydx,hstep,yOut,yErr);

    terror = GXNystromRK4_TruncatedError(&nrk4,hstep,yOut,yErr);
    dchord = GXNystromRK4_DistChord(&nrk4);

    //output
    track[tid].x    = yOut[0];
    track[tid].y    = yOut[1];
    track[tid].z    = yOut[2];
    track[tid].px   = yOut[3];
    track[tid].py   = yOut[4];
    track[tid].pz   = yOut[5];
    track[tid].s    = terror + dchord; //temporary
  }
}
//-----------------------------------------------------------------------------
// Common
//-----------------------------------------------------------------------------
#include "GPThreeVector.cu"
#include "GPThreeVectorList.cu"
#include "GPRotationMatrix.cu"
#include "GPUtils.cu"
#include "GPVParticleChange.cu"
#include "GPPhysicsTable.cu"
#include "GPPhysics2DVector.cu"
#include "GXPhysicsTable.cu"
#include "GXPhysics2DVector.cu"

//-----------------------------------------------------------------------------
// Material
//-----------------------------------------------------------------------------
#include "GPElement.cu"
#include "GPMaterial.cu"
#include "GPStep.cu"
#include "GPStepPoint.cu"
#include "GPAtomicShells.cu"
#include "GPSandiaTable.cu"
#include "GPIonisParamMat.cu"
//-----------------------------------------------------------------------------
// Electron processes/models
//-----------------------------------------------------------------------------
#include "GXeBremsstrahlung.cu"
#include "GXeIonisation.cu"
#include "GXeMultipleScattering.cu"

#include "GPeBremsstrahlung.cu"
#include "GPeIonisation.cu"
#include "GPeMultipleScattering.cu"
#include "GPSeltzerBergerRelModel.cu"
#include "GPMollerBhabhaModel.cu"
#include "GPUniversalFluctuation.cu"
#include "GPUrbanMscModel95.cu"
//-----------------------------------------------------------------------------
// Photon processes/models
//-----------------------------------------------------------------------------
#include "GPPhotonModel.cu"
#include "GPPhotonProcess.cu"

//-----------------------------------------------------------------------------
// Manager
//-----------------------------------------------------------------------------
#include "GXProcessManager.cu"
#include "GXTrackingManager.cu"
#include "GXSteppingManager.cu"
//-----------------------------------------------------------------------------
// GPGeometry
//-----------------------------------------------------------------------------
#include "GPAffineTransform.cu"
#include "GPAuxiliaryNavServices.cu"
#include "GPCombinedNavigation.cu"
#include "GPLineSection.cu"
#include "GPLogicalVolume.cu"
#include "GPNavigationLevel.cu"
#include "GPNavigator.cu"
#include "GPNavigationHistory.cu"
#include "GPNormalNavigation.cu"
#include "GPSmartVoxelHeader.cu"
#include "GPSmartVoxelNode.cu"
#include "GPSmartVoxelProxy.cu"
#include "GPTouchableHistory.cu"
#include "GPVoxelLimits.cu"
#include "GPVoxelNavigation.cu"
#include "GPVPhysicalVolume.cu"
#include "GPBox.cu"
#include "GPCons.cu"
#include "GPOrb.cu"
#include "GPTrd.cu"
#include "GPTubs.cu"
#include "GPVSolid.cu"
#include "GPChargeState.cu"
#include "GPVoxelHeader.cu"

//-----------------------------------------------------------------------------
// Magnetic Field and Transporation
//-----------------------------------------------------------------------------
#include "GXFieldTrack.cu"
#include "GXFieldManager.cu"
#include "GXMagneticField.cu"
#include "GXClassicalRK4.cu"
#include "GXCashKarpRKF45.cu"
#include "GXNystromRK4.cu"
#include "GXEquationOfMotion.cu"
#include "GXMagIntegratorDriver.cu"
#include "GXChordFinder.cu"
#include "GXPropagatorInField.cu"
#include "GXMultiLevelLocator.cu"
#include "GXTransportation.cu"

//-----------------------------------------------------------------------------
// FrameWork and Application
//-----------------------------------------------------------------------------
#include "GXRunManager.cc"
#include "GXUserGeometry.cc"
#include "GXSimpleEcal.cc"