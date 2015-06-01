//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: geomdefs.hh,v 1.8 2007-05-18 07:21:39 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Constants, typedefs, enums for Geometry Section
//
// History:
// 30.06.95 P.Kent

//enum ESolid { kBox, kTrd };
//enum ESolid { kBox, kOrb, kTubs, kCons, kPolyCone, kTrd };

#ifndef GeomDefs_hh
#define GeomDefs_hh

//#include "globals.hh"
#include "GPTypeDef.h"
#include "GPConstants.h"

// `Infinity' - Distance returned for no intersection etc.
CONSTTYPE static const G4double kInfinity = 9.0E99;

//G4GeometryTolerance
CONSTTYPE static const G4double kCarTolerance = 1E-9;//*mm;
CONSTTYPE static const G4double kAngTolerance = 1E-9;//*rad;
CONSTTYPE static const G4double kRadTolerance = 1E-9;//*mm;

//geometry/management/include/meshdefs.hh
CONSTTYPE const G4double kMeshAngleDefault=(pi/4); // Angle for mesh `wedges' in rads
CONSTTYPE const G4int kMinMeshSections=3;  // Min wedges+1 to make
CONSTTYPE const G4int kMaxMeshSections=37; // max wedges+1 to make

// Minimum cosine of angle between surface normal & track direction
// for exiting normal optimisation
CONSTTYPE static const G4double kMinExitingNormalCosine = 1E-3;

// Define axes for function params etc.
// X/Y/ZAxis = Normal Catesian axes
// Rho = Radial axis in cylindrical polar
// Radial3D = Radial axis in spherical polar
// Phi = Phi axis in cylindrical polar
// kUndefined = Not defined axis
enum EAxis {kXAxis,kYAxis,kZAxis,kRho,kRadial3D,kPhi,kUndefined};

// G4VSolid::Inside return codes
// kSurface => within tolerance of exact surface
enum EInside {kOutside,kSurface,kInside};

// kNormal = (G4PVPlacement) Conventional positioning
// kReplica = (G4PVReplica)  Consumed parameterised case
//                           => Distances & location computed with
//                              simple formulae & MOTHER volume(s)
//                              must also be checked
// kParameterised = (G4PVParameterised) General parameterised volume
//                           => Distance & location computed to volumes
//                              after setup/modification via user object
enum EVolume {kNormal,kReplica,kParameterised};

// Default max size of Navigation history
//
//static const G4int kHistoryMax    = 16;
CONSTTYPE static const G4int kHistoryMax    = 8;

// History increase stride of Navigation history
//
//static const G4int kHistoryStride = 16;
CONSTTYPE static const G4int kHistoryStride = 8;

// Voxel stack depth maximum [no resizing]
//
CONSTTYPE static const G4int kNavigatorVoxelStackMax = 3;
// Voxel stack depth maximum [no resizing]
//
//static const G4int max_depth = 10;
CONSTTYPE static const G4int max_depth = 5;

#endif /* GeomDefs_hh */