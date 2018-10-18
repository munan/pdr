//#include "../copyright.h"
/*==============================================================================
 * FILE: cool_tigress.c
 *
 * PURPOSE: heating and cooling function for TIGRESS simulations.
 *
 * REFERENCES: Gong, Ostriker and wolfire (2017)
 *
 * METHODS:
 * Estimate chemical abundances from gas and radiation parameters based on the
 * equilibrium chemistry in Gong, Ostriker and Wolfire (2017) network. Then
 * calculate:
 * Heating: cosmic-ray (CR) heating, photo-electric heating on dust (PE), UV
 * pumping of H2
 * Cooling: Ly-alpha, OI, C+, CI, CO, recombination of e- on PAHs
 * 
 * NOMENCLATURE:
 *
 *   Written by Munan Gong at 2018-10-15
 *
 * GLOBAL DEFAULT PARAMETERS
 * xCstd = 1.6e-4 (carbon abundance at solar metallicity)
 * xOstd = 3.2e-4 (oxygen abundance at solar metallicity)
 * xHetot = 0.1 (helium abundance)
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * cooling(nH, T, dvdr, Z, xi_CR, G_PE, G_CI, G_CO(, G_H2))
 * heating(nH, T, Z, G_PE)
 *
 * Z = Z_d = Z_g is the same for gas and dust metallicity
 *
 * CONTAINS PRIVATE FUNCTIONS:
 * chemical abundances:
 * fH2_CR(nH, T, Z_d, xi_CR)
 * fH2_CR_FUV(nH, T, Z_d, xi_CR, G_H2)
 * fCplus(x_e, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI)
 * fHplus(x_e, x_H2, x_Cplus, nH, T, Z_d, xi_CR, G_PE)
 * fe_e(x_e, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI)
 * fe(x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI)
 * fCO(x_H2, x_Cplus, nH, Z_d, Z_g, xi_CR, G_CO)
 *
 * Heating:
 * heatingPE(x_e, nH, T, Z_d, G_PE)
 * heatingCR(x_e, x_HI, x_H2, nH, xi_CR)
 * heatingH2pump(x_HI, x_H2, nH, T, G_H2)
 * 
 * Cooling:
 * coolingLya(x_HI, x_e, nH, T)
 * coolingOI(x_OI, x_HI, x_H2, x_e, nH, T)
 * coolingCII(x_Cplus, x_HI, x_H2, x_e, nH, T)
 * coolingCI(x_CI, x_HI, x_H2, x_e, nH, T)
 * coolingCO(x_CO, x_HI, x_H2, x_e, nH, T, dvdr)
 * coolingRec(x_e, nH, T, Z_d, G_PE)
 *
 *============================================================================*/

#include <math.h>
#define Real double //TODO:replace this by include def.in file in athena

static const Real xCstd=1.6e-4, xOstd=3.2e-4, xHetot=0.1;

/*----------------------------------------------------------------------------*/
/* PRIVATE FUCNTIONS                                                          */
/*----------------------------------------------------------------------------*/
//H2 abundance consider only CR destruction (no FUV)
static Real fH2_CR(const Real nH, const Real T, const Real Z_d,
                   const Real xi_CR);
static Real fH2_CR_FUV(const Real nH, const Real T, const Real Z_d,
                       const Real xi_CR, const Real G_H2);

/*----------------------------------------------------------------------------*/
/* PUBLIC FUCNTIONS                                                           */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* IMPLEMENTATION of FUCNTIONS                                                */
/*----------------------------------------------------------------------------*/

Real fH2_CR(const Real nH, const Real T, const Real Z_d, const Real xi_CR) {
  Real kgr = 3.0e-17*Z_d;
  Real R = kgr * nH / (2. * xi_CR);
  Real b = - (1.65*1.5 + 2.*R);
  Real a = 1.65*0.7;
  Real c = R;
  Real x_H2 = (-b - sqrt(b*b - 4*a*c) )/(2.*a);
  return x_H2;
}

Real fH2_CR_FUV(const Real nH, const Real T, const Real Z_d, const Real xi_CR, 
                const Real G_H2) {
  Real kgr = 3.0e-17*Z_d;
  Real R = kgr * nH / (2. * xi_CR);
  Real a = 1.65*0.7;
  Real c = R;
  Real k_FUV = 5.7e-11 * G_H2;
  Real b = - (1.65*1.5 + 2.*R + k_FUV/(2.*xi_CR));
  Real x_H2 = (-b - sqrt(b*b - 4*a*c) )/(2.*a);
  return x_H2;
}
