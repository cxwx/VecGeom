// -*- C++ -*-
// $Id:$
// ----------------------------------------------------------------------
// HEP coherent Physical Constants
//
// This file has been provided by Geant4 (simulation toolkit for HEP).
//
// The basic units are :
//  		millimeter  
// 		nanosecond  
// 		Mega electron Volt  
// 		positon charge 
// 		degree Kelvin
//              amount of substance (mole)
//              luminous intensity (candela)
// 		radian  
//              steradian 
//
// Below is a non exhaustive list of Physical CONSTANTS,
// computed in the Internal HEP System Of Units.
//
// Most of them are extracted from the Particle Data Book :
//        Phys. Rev. D  volume 50 3-1 (1994) page 1233
// 
//        ...with a meaningful (?) name ...
//
// You can add your own constants.
//
// Author: M.Maire
//
// History:
//
// 23.02.96 Created
// 26.03.96 Added constants for standard conditions of temperature
//          and pressure; also added Gas threshold.
// 29.04.08   use PDG 2006 values
// 03.11.08   use PDG 2008 values

#ifndef HEP_PHYSICAL_CONSTANTS_H
#define HEP_PHYSICAL_CONSTANTS_H

#include "GPTypeDef.h"
#include "GPSystemOfUnits.h"


//
//
//
CONSTTYPE const G4double     pi  = 3.14159265358979323846;
CONSTTYPE const G4double  twopi  = 2*pi;
//CONSTTYPE const G4double halfpi  = pi/2;
//CONSTTYPE const G4double     pi2 = pi*pi;

//
// 
//
CONSTTYPE const G4double Avogadro = 6.02214179e+23/mole;

//
// c   = 299.792458 mm/ns
// c^2 = 898.7404 (mm/ns)^2 
//
CONSTTYPE const G4double c_light   = 2.99792458e+8 * m/s;
CONSTTYPE const G4double c_squared = c_light * c_light;

//
// h     = 4.13566e-12 MeV*ns
// hbar  = 6.58212e-13 MeV*ns
// hbarc = 197.32705e-12 MeV*mm
//
CONSTTYPE const G4double h_Planck      = 6.62606896e-34 * joule*s;
CONSTTYPE const G4double hbar_Planck   = h_Planck/twopi;
CONSTTYPE const G4double hbarc         = hbar_Planck * c_light;
//CONSTTYPE const G4double hbarc_squared = hbarc * hbarc;

//
//
//
CONSTTYPE const G4double electron_charge = - eplus; // see SystemOfUnits.h
CONSTTYPE const G4double e_squared = eplus * eplus;

//
// amu_c2 - atomic equivalent mass unit
//        - AKA, unified atomic mass unit (u)
// amu    - atomic mass unit
//
CONSTTYPE const G4double electron_mass_c2 = 0.510998910 * MeV;
CONSTTYPE const G4double   proton_mass_c2 = 938.272013 * MeV;
CONSTTYPE const G4double  neutron_mass_c2 = 939.56536 * MeV;
CONSTTYPE const G4double           amu_c2 = 931.494028 * MeV;
CONSTTYPE const G4double              amu = amu_c2/c_squared;

//
// permeability of free space mu0    = 2.01334e-16 Mev*(ns*eplus)^2/mm
// permittivity of free space epsil0 = 5.52636e+10 eplus^2/(MeV*mm)
//
CONSTTYPE const G4double mu0      = 4*pi*1.e-7 * henry/m;
CONSTTYPE const G4double epsilon0 = 1./(c_squared*mu0);

//
// electromagnetic coupling = 1.43996e-12 MeV*mm/(eplus^2)
//
CONSTTYPE const G4double elm_coupling           = e_squared/(4*pi*epsilon0);
CONSTTYPE const G4double fine_structure_const   = elm_coupling/hbarc;
CONSTTYPE const G4double classic_electr_radius  = elm_coupling/electron_mass_c2;
CONSTTYPE const G4double electron_Compton_length = hbarc/electron_mass_c2;
//CONSTTYPE const G4double Bohr_radius = electron_Compton_length/fine_structure_const;

CONSTTYPE const G4double alpha_rcl2 = fine_structure_const
                                   *classic_electr_radius
                                   *classic_electr_radius;

CONSTTYPE const G4double twopi_mc2_rcl2 = twopi*electron_mass_c2
                                             *classic_electr_radius
                                             *classic_electr_radius;
//
//
//
CONSTTYPE const G4double k_Boltzmann = 8.617343e-11 * MeV/kelvin;

//
//
//
//CONSTTYPE const G4double STP_Temperature = 273.15*kelvin;
//CONSTTYPE const G4double STP_Pressure    = 1.*atmosphere;
//CONSTTYPE const G4double kGasThreshold   = 10.*mg/cm3;

//
//
//
//CONSTTYPE const G4double universe_mean_density = 1.e-25*g/cm3;


#endif /* HEP_PHYSICAL_CONSTANTS_H */