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
