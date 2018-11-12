//#include "../copyright.h"
//TODO: uncomment the copyright include
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
 * heating(x_e, x_HI, x_H2, nH, T, Z, xi_CR, G_PE, G_H2)
 * cooling(x_e, x_HI, x_H2, x_Cplus, x_CI, x_CO, x_OI, 
           nH, T, dvdr, Z, G_PE)
 *
 * Notes:
 * - Input parameters nH, T, dvdr are in CGS units, radiation field G_i in
 * units of Draine 1978 radiation field, xi_CR in [s^{-1}H^{-1}] (per second per
 * H atom).
 * - Output heating and cooling rates in [ergs s^{-1} cm^{-3} ] units.
 * TODO: multiply
 * by density, since the thermo.cpp functions output in [ergs/s/H] units.
 * - Z = Z_d = Z_g is the same for gas and dust metallicity
 *
 * CONTAINS PRIVATE FUNCTIONS:
 * chemical abundances:
 * fH2(nH, T, Z_d, xi_CR, G_H2)
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
 * coolingLya(x_e, x_HI, nH, T)
 * coolingOI(x_e, x_OI, x_HI, x_H2, nH, T)
 * coolingCII(x_e, x_Cplus, x_HI, x_H2, nH, T)
 * coolingCI(x_e, x_CI, x_HI, x_H2, nH, T)
 * coolingCO(x_e, x_CO, x_HI, x_H2, nH, T, dvdr)
 * coolingRec(x_e, nH, T, Z_d, G_PE)
 *
 * Notes:
 * The heating and cooling rates above are in [ergs s^{-1} H^{-1}] units. 
 * In order to get the heating and cooling rates for TIGRESS in
 * [ergs s^{-1} cm^{-3}] units, they are multiplied by density nH in functions
 * heating() and cooling().
 *
 * helper functions for cooling:
 * cooling2Level_()
 * cooling3Level_()
 * linearInterpIndex_()
 * linearInterp_()
 * LP1Di_()
 * LP2Di_()
 *
 *============================================================================*/

#include <math.h>
#define Real double //TODO:replace this by include def.in file in athena

/*elemental abundances*/
static const Real xCstd=1.6e-4, xOstd=3.2e-4, xHetot=0.1;
//TODO: check constants for heating and cooling, delete the unused ones.
/*physical constants*/
static const Real eV_ = 1.602e-12;
static const Real kb_ = 1.381e-16;
static const Real ca_ = 2.27e-4;
static const Real TCMB_ = 2.73;
/*ortho to para ratio of H2*/
static const Real o2p_ = 3.;
static const Real fo_ = 0.75;
static const Real fp_ = 0.25;
static const Real sigmaPE_ = 1.0e-21;
static const Real sigmaISRF_ = 3.0e-22;
static const Real sigmad10_ = 2.0e-25;
static const Real alpha_GD_ = 3.2e-34;
/*----C+, 2 level system---*/
static const Real A10CII_ = 2.3e-6;
static const Real E10CII_ = 1.26e-14;
static const Real g0CII_ = 2.;
static const Real g1CII_ = 4.;
/*----HI+, 2 level system---*/
static const Real A10HI_ = 6.265e8;
static const Real E10HI_ = 1.634e-11;
static const Real g0HI_ = 1.;
static const Real g1HI_ = 3.;
/*----CI, 3 level system---*/
static const Real g0CI_ = 1;
static const Real g1CI_ = 3;
static const Real g2CI_ = 5;
static const Real A10CI_ = 7.880e-08;
static const Real A20CI_ = 1.810e-14;
static const Real A21CI_ = 2.650e-07;
static const Real E10CI_ = 3.261e-15;
static const Real E20CI_ = 8.624e-15;
static const Real E21CI_ = 5.363e-15;
/*----OI, 3 level system---*/
static const Real g0OI_ = 5;
static const Real g1OI_ = 3;
static const Real g2OI_ = 1;
static const Real A10OI_ = 8.910e-05;
static const Real A20OI_ = 1.340e-10;
static const Real A21OI_ = 1.750e-05;
static const Real E10OI_ = 3.144e-14;
static const Real E20OI_ = 4.509e-14;
static const Real E21OI_ = 1.365e-14;

/*-----CO cooling table data, from Omukai+2010-----*/
static const Real TCO_[lenTCO_] = {10,	20,	30,	50,	80,	100,
                                   300,	600,	1000,	1500,	2000};
static const Real NeffCO_[lenNeffCO_] = {14.0, 14.5, 15.0, 15.5, 16.0, 16.5,
                                         17.0, 17.5, 18.0, 18.5, 19.0};
static const Real L0CO_[lenTCO_] = {24.77,	24.38, 24.21,	24.03, 23.89, 23.82,
  23.34238089,  22.99832519,  22.75384686,  22.56640625, 22.43740866};
static const Real LLTECO_[lenNeffCO_*lenTCO_] = {
21.08, 20.35, 19.94, 19.45, 19.01, 18.80, 17.81, 17.23, 16.86, 16.66, 16.55,
21.09, 20.35, 19.95, 19.45, 19.01, 18.80, 17.81, 17.23, 16.86, 16.66, 16.55,
21.11, 20.37, 19.96, 19.46, 19.01, 18.80, 17.81, 17.23, 16.86, 16.66, 16.55,
21.18, 20.40, 19.98, 19.47, 19.02, 18.81, 17.82, 17.23, 16.87, 16.66, 16.55,
21.37, 20.51, 20.05, 19.52, 19.05, 18.83, 17.82, 17.23, 16.87, 16.66, 16.55,
21.67, 20.73, 20.23, 19.64, 19.13, 18.90, 17.85, 17.25, 16.88, 16.67, 16.56,
22.04, 21.05, 20.52, 19.87, 19.32, 19.06, 17.92, 17.28, 16.90, 16.69, 16.58,
22.44, 21.42, 20.86, 20.19, 19.60, 19.33, 18.08, 17.38, 16.97, 16.75, 16.63,
22.87, 21.82, 21.24, 20.55, 19.95, 19.66, 18.34, 17.59, 17.15, 16.91, 16.78,
23.30, 22.23, 21.65, 20.94, 20.32, 20.03, 18.67, 17.89, 17.48, 17.26, 17.12,
23.76, 22.66, 22.06, 21.35, 20.71, 20.42, 19.03, 18.26, 17.93, 17.74, 17.61
};
static const Real nhalfCO_[lenNeffCO_*lenTCO_] = {
  3.29, 3.49 ,3.67  ,3.97,  4.30, 4.46, 5.17, 5.47, 5.53, 5.30, 4.70, 
  3.27, 3.48 ,3.66  ,3.96,  4.30, 4.45, 5.16, 5.47, 5.53, 5.30, 4.70, 
  3.22, 3.45 ,3.64  ,3.94,  4.29, 4.45, 5.16, 5.47, 5.53, 5.30, 4.70, 
  3.07, 3.34 ,3.56  ,3.89,  4.26, 4.42, 5.15, 5.46, 5.52, 5.30, 4.70, 
  2.72, 3.09 ,3.35  ,3.74,  4.16, 4.34, 5.13, 5.45, 5.51, 5.29, 4.68, 
  2.24, 2.65 ,2.95  ,3.42,  3.92, 4.14, 5.06, 5.41, 5.48, 5.26, 4.64, 
  1.74, 2.15 ,2.47  ,2.95,  3.49, 3.74, 4.86, 5.30, 5.39, 5.17, 4.53, 
  1.24, 1.65 ,1.97  ,2.45,  3.00, 3.25, 4.47, 5.02, 5.16, 4.94, 4.27, 
 0.742, 1.15 ,1.47  ,1.95,  2.50, 2.75, 3.98, 4.57, 4.73, 4.52, 3.84, 
 0.242, 0.652,0.966 ,1.45,  2.00, 2.25, 3.48, 4.07, 4.24, 4.03, 3.35, 
-0.258, 0.152,0.466 ,0.95,	 1.50, 1.75, 2.98, 3.57, 3.74, 3.53, 2.85
};
static const Real alphaCO_[lenNeffCO_*lenTCO_] = {
0.439, 0.409, 0.392, 0.370, 0.361, 0.357, 0.385, 0.437, 0.428, 0.354, 0.322,	
0.436, 0.407, 0.391, 0.368, 0.359, 0.356, 0.385, 0.437, 0.427, 0.354, 0.322,	
0.428, 0.401, 0.385, 0.364, 0.356, 0.352, 0.383, 0.436, 0.427, 0.352, 0.320,	
0.416, 0.388, 0.373, 0.353, 0.347, 0.345, 0.380, 0.434, 0.425, 0.349, 0.316,	
0.416, 0.378, 0.360, 0.338, 0.332, 0.330, 0.371, 0.429, 0.421, 0.341, 0.307,	
0.450, 0.396, 0.367, 0.334, 0.322, 0.317, 0.355, 0.419, 0.414, 0.329, 0.292,	
0.492, 0.435, 0.403, 0.362, 0.339, 0.329, 0.343, 0.406, 0.401, 0.317, 0.276,	
0.529, 0.473, 0.441, 0.404, 0.381, 0.370, 0.362, 0.410, 0.392, 0.316, 0.272,	
0.555, 0.503, 0.473, 0.440, 0.423, 0.414, 0.418, 0.446, 0.404, 0.335, 0.289,	
0.582, 0.528, 0.499, 0.469, 0.457, 0.451, 0.470, 0.487, 0.432, 0.364, 0.310,	
0.596, 0.546, 0.519, 0.492, 0.483, 0.479, 0.510, 0.516, 0.448, 0.372, 0.313
};
static const Real logTg_[lenTg_] = {
  0.5       , 0.88888889, 1.27777778, 1.66666667, 2.05555556,  2.44444444,
  2.83333333, 3.22222222, 3.61111111, 4.
};
static const Real lognH_[lennH_] = {
  0.        , 0.42857143, 0.85714286, 1.28571429, 1.71428571, 2.14285714,
  2.57142857, 3.        , 3.42857143, 3.85714286, 4.28571429, 4.71428571,
  5.14285714, 5.57142857, 6.        
};
static const Real logps_[lennH_ * lenTg_] = {
   33.60923439, 32.35048647, 31.6458604 , 31.02132235, 30.42222289 ,
   29.83261   , 29.24673384, 28.66235604, 28.0785789 , 27.49504472 ,
   33.18091039, 31.92216147, 31.21753302, 30.59298934, 29.99387726 ,
   29.40423904, 28.8183211 , 28.23389204, 27.65007024, 27.06650687 ,
   32.75300149, 31.49424423, 30.78959638, 30.16501048, 29.5658181  ,
   28.97605855, 28.39000596, 27.80546791, 27.22157705, 26.6379757  ,
   32.32619855, 31.06738031, 30.36260198, 29.73777592, 29.13823862 ,
   28.54811853, 27.96178914, 27.37708183, 26.79309991, 26.20945204 ,
   31.90230918, 30.64308622, 29.93758718, 29.3117745 , 28.71125923 ,
   28.12042115, 27.53366605, 26.94873467, 26.36464058, 25.78093706 ,
   31.48587783, 30.22433325, 29.51586354, 28.88730605, 28.28488476 ,
   27.69295575, 27.10563836, 26.52043033, 25.93620173, 25.35243225 ,
   31.0872853 , 29.81525158, 29.09829079, 28.46439304, 27.85909326 ,
   27.26572717, 26.67771529, 26.09217504, 25.50778679, 24.92393938 ,
   30.72591396, 29.41932635, 28.68506438, 28.0430074 , 27.4339034  ,
   26.83875916, 26.24991203, 25.66397699, 25.07939997, 24.49546057 ,
   30.42733764, 29.03861249, 28.27630163, 27.62323376, 27.00938204 ,
   26.41209045, 25.82224861, 25.23584618, 24.65104624, 24.06699834 ,
   30.21146872, 28.67535786, 27.87248634, 27.20529029, 26.58563564 ,
   25.9857724 , 25.39474974, 24.80779463, 24.22273154, 23.63855567 ,
   30.08004755, 28.33375429, 27.47457781, 26.78951505, 26.16280582 ,
   25.55986866, 24.96744513, 24.37983659, 23.79446288, 23.21013604 ,
   30.0136033 , 28.02095875, 27.08405993, 26.37636621, 25.74107062 ,
   25.13445646, 24.54037027, 23.95198897, 23.36624852, 22.7817436  ,
   29.98466705, 27.74775444, 26.70307122, 25.96644194, 25.32065037 ,
   24.70962891, 24.1135674 , 23.52427179, 22.93809825, 22.35338321 ,
   29.97312226, 27.52792637, 26.3346648 , 25.56052005, 24.90181737 ,
   24.28549826, 23.68708687, 23.09670875, 22.51002362, 21.92506064 ,
   29.96870043, 27.37312455, 25.98324711, 25.15962125, 24.48490972 ,
   23.86220014, 23.26098873, 22.66932799, 22.08203828, 21.49678269 ,
};

static const Real CPE_[7] = {5.22, 2.25, 0.04996, 0.00430, 0.147, 0.431,0.692};
static const Real DPE_[5] = {0.4535, 2.234, -6.266, 1.442, 0.05089};

/*----------------------------------------------------------------------------*/
/* PUBLIC FUCNTIONS                                                           */
/*----------------------------------------------------------------------------*/
static Real heating(const Real x_e, const Real x_HI, const Real x_H2,
                    const Real nH, const Real T, const Real Z,
                    const Real xi_CR, const Real G_PE, const Real G_H2);
static Real cooling(const Real x_e, const Real x_HI, const Real x_H2,
                    const Real x_Cplus, const Real x_CI,
                    const Real x_CO, const Real x_OI, 
                    const Real nH, const Real T, const Real dvdr,
                    const Real Z, const Real G_PE);

/*----------------------------------------------------------------------------*/
/* PRIVATE FUCNTIONS                                                          */
/*----------------------------------------------------------------------------*/
//Chemistry-------------------------------------------------------------------
//H2 abundance consider CR destruction and FUV
static Real fH2(const Real nH, const Real T, const Real Z_d,
                const Real xi_CR, const Real G_H2);
//C+ abundance, depending on abundances of e- and H2
static Real fCplus(const Real x_e, const Real x_H2, const Real nH, const Real T,
                   const Real Z_d, const Real Z_g, const Real xi_CR, 
                   const Real G_PE, const Real G_CI);
//H+ abundance, depending on abundances of C+, e- and H2
static Real fHplus(const Real x_e, const Real x_Cplus, const Real x_H2,
                   const Real nH, const Real T, const Real Z_d, const Real xi_CR,
                   const Real G_PE);
//ion abundances: ion = C+ + H+, depending on e- and H2 abundances.
static Real fions(const Real x_e, const Real x_H2, const Real nH, const Real T,
                  const Real Z_d, const Real Z_g, const Real xi_CR, 
                  const Real G_PE, const Real G_CI);
//e- abundances, solved from iteration, assume e- balances from C+ and H+
static Real fe(const Real x_H2, const Real nH, const Real T,
               const Real Z_d, const Real Z_g, const Real xi_CR, 
               const Real G_PE, const Real G_CI);
//CO abundance. Fit from Gong, Ostriker and Wolfire 2017
//Note: This is only tested in the case Z_d = Z_g
static Real fCO(const Real x_H2, const Real x_Cplus, const Real nH, 
                const Real Z_d, const Real Z_g, const Real xi_CR, const Real G_CO);

//heating---------------------------------------------------------------------
//cosmic ray heating
static Real heaingCR(const Real x_e, const Real x_HI, const Real x_H2,  
                     const Real nH, const Real xi_CR);
//photo electric heating on dust
static Real heatingPE(const Real x_e, const Real nH, const Real T, 
                      const Real Z_d, const Real G_PE);
//UV-pumping of H2
static Real heatingH2pump(const Real x_HI, const Real x_H2, const Real nH,
                          const Real T, const Real G_H2);

//cooling---------------------------------------------------------------------
//HI Lyman alpha cooling
static Real coolingLya(const Real x_e, const Real x_HI, const Real nH,
                       const Real T);
//OI cooling
static Real coolingOI(const Real x_e, const Real x_OI, const Real x_HI,
                      const Real x_H2, const Real nH, const Real T);
//C+ cooling
static Real coolingCII(const Real x_e, const Real x_Cplus, const Real x_HI,
                       const Real x_H2, const Real nH, const Real T);
//CI cooling
static Real coolingCI(const Real x_e, const Real x_CI, const Real x_HI,
                      const Real x_H2, const Real nH, const Real T);
//CO rotational line cooling, dvdr in cgs units
static Real coolingCO(const Real x_e, const Real x_CO, const Real x_HI,
                      const Real x_H2, const Real nH, const Real T, 
                      const Real dvdr);
//cooling by recombination of e on PAHs
static Real coolingRec(const Real x_e, const Real nH, const Real T, 
                       const Real Z_d, const Real G_PE);

//helper functions for cooling-------------------------------------------------
static Real cooling2Level_(const Real q01, const Real q10, const Real A10,
                           const Real E10, const Real xs);

static Real cooling3Level_(const Real q01, const Real q10, const Real q02,
                           const Real q20, const Real q12, const Real q21,
													 const Real A10, const Real A20, const Real A21,
                           const Real E10, const Real E20, const Real E21,
													 const Real xs);
static int linearInterpIndex_(const int len, const Real xarr[], const Real x);
static Real linearInterp_(const Real x0, const Real x1, const Real y0,
                          const Real y1, const Real x);
static Real LP1Di_(const Real *xarr, const Real *data, const int ix,
                   const Real x);
static Real LP2Di_(const Real *xarr, const Real *yarr,
                   const int lenx, const int ix, const int iy,
                   const Real *data, const Real x, const Real y);


/*----------------------------------------------------------------------------*/
/* IMPLEMENTATION of FUCNTIONS                                                */
/*----------------------------------------------------------------------------*/

Real fH2(const Real nH, const Real T, const Real Z_d, const Real xi_CR, 
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

Real fCplus(const Real x_e, const Real x_H2, const Real nH, const Real T,
            const Real Z_d, const Real Z_g, const Real xi_CR, 
            const Real G_PE, const Real G_CI) {
  const Real small_ = 1e-50;
  Real k_C_cr = 3.85 * xi_CR;
  Real k_C_photo = 3.5e-10*G_CI;
  Real k_Cplus_e = CII_rec_rate_(T);
  Real psi_gr = 1.7 * G_PE * sqrt(T)/(nH * x_e + small_);
  const Real cCp_[7] = {45.58, 6.089e-3, 1.128, 4.331e2, 4.845e-2,
                        0.8120, 1.333e-4};
  Real k_Cplus_gr = 1.0e-14 * cCp_[0] / 
		           (
			           1.0 + cCp_[1]*pow(psi_gr, cCp_[2]) * 
								   (1.0 + cCp_[3] * pow(T, cCp_[4])
										             *pow( psi_gr, -cCp_[5]-cCp_[6]*log(T) ) 
									 ) 
								) * Z_d;
  Real k_Cplus_H2 = 3.3e-13 * pow(T, -1.3) * exp(-23./T);
  Real c = (k_C_cr + k_C_photo) / nH;
  Real al = k_Cplus_e*x_e + k_Cplus_gr + k_Cplus_H2*x_H2 + c;
  Real ar = xCtot_ * c;
  Real x_Cplus = ar / al;
  return x_Cplus;
}

Real fHplus(const Real x_e, const Real x_Cplus, const Real x_H2,
            const Real nH, const Real T, const Real Z_d, const Real xi_CR, 
            const Real G_PE) {
  const Real small_ = 1e-50;
  Real x_H = 1. - 2. * x_H2 - (x_e - x_Cplus);
  x_H = MAX(x_H, 0.0);
  Real k_Hplus_e = 2.753e-14 * pow( 315614.0 / T, 1.5) * pow( 
               1.0 + pow( 115188.0 / T, 0.407) , -2.242 );
  const Real cHp_[7] = {12.25, 8.074e-6, 1.378, 5.087e2,
                               1.586e-2, 0.4723, 1.102e-5};
  Real psi_gr = 1.7 * G_PE * sqrt(T)/(nH * x_e + small_);
  Real k_Hplus_gr = 1.0e-14 * cHp_[0] / 
             (
               1.0 + cHp_[1]*pow(psi_gr, cHp_[2]) * 
                 (1.0 + cHp_[3] * pow(T, cHp_[4])
                               *pow( psi_gr, -cHp_[5]-cHp_[6]*log(T) ) 
                 ) 
              ) * Z_d;
  Real k_cr_H = xi_CR * (2.3*x_H2 + 1.5*x_H);
  Real c = k_cr_H * x_H / nH;
  Real x_Hplus = c/(k_Hplus_e *  x_e + k_Hplus_gr);
  x_Hplus = MIN(x_Hplus, 1.0);
  return x_Hplus;
}

Real fions(const Real x_e, const Real x_H2, const Real nH, const Real T,
           const Real Z_d, const Real Z_g, const Real xi_CR, 
           const Real G_PE, const Real G_CI) {
  const Real x_Cplus = fCplus(x_e, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI);
  const Real x_Hplus = fHplus(x_e, x_Cplus, x_H2, nH, T, Z_d, xi_CR, G_PE);
  return (x_Cplus + x_Hplus);
}

Real fe(const Real x_H2, const Real nH, const Real T,
        const Real Z_d, const Real Z_g, const Real xi_CR, 
        const Real G_PE, const Real G_CI) { 
  const Real rtol = 1e-2;
  const Real small_x = 1e-6;
  const Real small_f = 1e-20;
  const int maxiter = 20;
  Real x = 0;
  Real f = 0;
  Real xnew = 0;
  Real fnew = 0;
  int niter = 0;
  Real xprev = 0.5;
  Real fprev = 0;
  Real a = 0.0;
  Real fa = 0.0;
  Real b = 1.0;
  bool flag = true;
  niter = 0;
  fprev = fions(xprev, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI) - xprev;
  while (1) {
    x = fions(xprev, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI);
    if ( abs((x-xprev)/(xprev+small_x)) < rtol || abs(fprev) < small_x || 
         abs((a-b)/(b+small_x)) < rtol) {
      break;
    }
    if (niter > maxiter) {
      printf("fe(): WARNING: niter>maxiter(=%d), x=%.2e, xprev=%.2e\n", 
             maxiter, x, xprev);
      printf("a=%.2e, b=%.2e\n", a, b);
      //TODO: throw ath_error
      break;
    }
    f = fions(x, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI) - x;
    if ( abs(f - fprev) < small_f ) {
      xnew = (x + xprev)/2.;
    } else {
      xnew = x - f * (x-xprev)/(f-fprev);
    }
    if (xnew < a || xnew > b) {
      xnew = (a + b)/2.;
    }
    if (flag) {
      fa = fions(a, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI)-a;
    }
    fnew = fions(xnew, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI)-xnew;
    if (fa * fnew < 0.0) {
      b = xnew;
      flag = false;
    } else {
      a = xnew;
      flag = true;
    }
    xprev = xnew;
    fprev = fnew;
    niter = niter + 1;
  }
  return xprev;
}

Real fCO(const Real x_H2, const Real x_Cplus, const Real nH, 
         const Real Z_d, const Real Z_g, const Real xi_CR, const Real G_CO) {
  Real xCtot = xCstd * Z_g;
  Real kcr16 = xi_CR/1.0e-16;
  Real term1 = 4.0e3*Z_d/(kcr16*kcr16);
  Real ncrit2 = 0;
  Real x_CO = 0;
  Real x_CO_max1 = 0;
  Real x_CO_max2 = 0;
  term1 = MAX(term1, 1.);
  ncrit2 = 2.*pow(term1, pow(G_CO, 1./3.)) * (50.*kcr16/pow(Z_d, 1.4));
  if (nH >= ncrit2) {
    x_CO = 1.;
  } else {
    x_CO = nH/ncrit2;
  }
  x_CO *= xCtot;
  x_CO_max1 = xCtot - x_Cplus;
  x_CO_max2 = xCtot * x_H2 * 2.;
  x_CO = MIN(x_CO, x_CO_max1);
  x_CO = MIN(x_CO, x_CO_max2);
  return x_CO;
}

Real heaingCR(const Real x_e, const Real x_HI, const Real x_H2,
              const Real nH, const Real xi_CR) {
	/* ionization rate per H*/
  const Real kcr_H_fac = 1.15 * 2*x_H2 + 1.5 * x_HI;
  const Real kHI = xi_CR * kcr_H_fac;
  const Real kHe = xi_CR * 1.1;
  const Real kH2 = xi_CR * 2.0 * kcr_H_fac;
	const Real ktot = kHI*x_HI + kHe*xHetot + kH2*x_H2;
	/* heating rate per ionization in atomic region. 
	 * Draine ISM book eq (30.1)*/
  Real qHI;
  if (x_e > 1.0e-9) {
	  qHI = ( 6.5 + 26.4 * sqrt( x_e / (x_e+0.07) ) ) * eV_;
  } else { //prevent sqrt of small negative number
	  qHI =  6.5 * eV_;
  }

	/* Heating rate per ioniztion in molecular region.
	 * Despotic paper Appendix B*/
	Real qH2;
  const Real lognH = log10(nH);
  if (nH < 100.) { //prevent log of small negative number
    qH2 = 10. * eV_;
  } else if (lognH < 4) {
    qH2 = ( 10. + 3.*(lognH - 2.)/2. ) * eV_;
  } else if (lognH < 7) {
    qH2 = ( 13. + 4.*(lognH - 4.)/3. ) * eV_;
  } else if (lognH < 10) {
    qH2 = ( 17. + (lognH - 7.)/3. ) * eV_;
  } else {
    qH2 = 18. * eV_;
  }
	const Real qtot = x_HI*qHI + 2*x_H2*qH2;
	return (ktot*qtot);
}

Real heatingPE(const Real x_e, const Real nH, const Real T, 
               const Real Z_d, const Real G_PE) {
  const Real ne = x_e * nH;
  const Real x = 1.7 * G_PE * sqrt(T)/ne + 50.;
  const Real fac = ( CPE_[0] + CPE_[1]*pow(T, CPE_[4]) ) / 
    (
     1. + CPE_[2]*pow(x, CPE_[5]) * ( 1. + CPE_[3]*pow(x, CPE_[6]) ) 
     );
  const Real heating = 1.7e-26 * G_PE * Z_d * fac;
  return heating;
}

Real heatingH2pump(const Real x_HI, const Real x_H2, const Real nH,
                   const Real T, const Real G_H2) {
  const Real dot_xH2_photo = 5.7e-11 * G_H2 * x_H2;
  const Real de = 1.6 * x_HI * exp( - pow(400./T, 2) )
                    + 1.4 * x_H2 * exp( - (12000./(T + 1200.)));
  const Real ncr = 1.0e6 / sqrt(T) / de;
  const Real f = 1. / (1. + ncr/nH);
  return dot_xH2_photo * 9. * 2.2*f * eV_;
}

Real cooling2Level_(const Real q01, const Real q10, const Real A10,
                    const Real E10, const Real xs) {
	const Real f1 = q01 / (q01 + q10 + A10);
	return f1*A10*E10*xs;
}

Real cooling3Level_(const Real q01, const Real q10, const Real q02,
                    const Real q20, const Real q12, const Real q21,
										const Real A10, const Real A20, const Real A21,
                    const Real E10, const Real E20, const Real E21,
										const Real xs) {
	const Real R10 = q10 + A10;
	const Real R20 = q20 + A20;
	const Real R21 = q21 + A21;
	const Real a0 = R10*R20 + R10*R21 + q12*R20;
	const Real a1 = q01*R20 + q01*R21 + R21*q02;
	const Real a2 = q02*R10 + q02*q12 + q12*q01;
	const Real de = a0 + a1 + a2;
	const Real f1 = a1 / de;
	const Real f2 = a2 / de;
	return ( f1*A10*E10 + f2*(A20*E20 + A21*E21) )*xs;
}

int linearInterpIndex_(const int len, const Real xarr[], const Real x){
  int i = 0;
  if ( x < xarr[0]) {
    return 0;
  } else if ( x > xarr[len-1]) {
    return len-2;
  } else {
    for (i=0; x>xarr[i]; i++) {}
    return i-1;
  }
}

Real linearInterp_(const Real x0, const Real x1, const Real y0,
                   const Real y1, const Real x){
  return y0 + ( (y1-y0)/(x1-x0) ) * (x-x0);
}

Real LP1Di_(const Real *xarr, const Real *data, const int ix,
            const Real x) {
  return linearInterp_(xarr[ix], xarr[ix+1], data[ix], data[ix+1], x);
}

Real LP2Di_(const Real *xarr, const Real *yarr,
            const int lenx, const int ix, const int iy,
            const Real *data, const Real x, const Real y) {
  Real fl1, fl2;
  const Real x0 = xarr[ix];
  const Real x1 = xarr[ix+1];
  fl1 = linearInterp_(x0, x1, data[iy*lenx + ix], data[iy*lenx + ix+1], x);
  fl2 = linearInterp_(x0, x1, data[(iy+1)*lenx + ix], data[(iy+1)*lenx + ix+1], x);
  return linearInterp_(yarr[iy], yarr[iy+1], fl1, fl2, y);
}

Real coolingLya(const Real x_e, const Real x_HI, const Real nH, const Real T) {
  const Real ne = x_e * nH;
  const Real A = 6.3803e-9;
  const Real beta = 1.17;
  const Real T4 = T / 1.0e4;
  const Real fac = A * pow(T4, beta);
  const Real k01e = fac * exp(-11.84/T4);
  const Real q01 = k01e * ne;
  const Real q10 = (g0HI_/g1HI_) * fac * ne;
	return cooling2Level_(q01, q10, A10HI_, E10HI_, x_HI);
}

Real coolingOI(const Real x_e, const Real x_OI, const Real x_HI,
               const Real x_H2, const Real nH, const Real T) {
  const Real nHI = x_HI * nH;
  const Real nH2 = x_H2 * nH;
  const Real ne = x_e * nH;
  /*collisional rates from  Draine (2011) ISM book Appendix F Table F.6*/
  const Real T2 = T/100;
  const Real lnT2 = log(T2);
  /*HI*/
  const Real k10HI = 3.57e-10 * pow(T2, 0.419-0.003*lnT2);
  const Real k20HI = 3.19e-10 * pow(T2, 0.369-0.006*lnT2);
  const Real k21HI = 4.34e-10 * pow(T2, 0.755-0.160*lnT2);
  /*H2*/
  const Real k10H2p = 1.49e-10 * pow(T2, 0.264+0.025*lnT2);
  const Real k10H2o = 1.37e-10 * pow(T2, 0.296+0.043*lnT2);
  const Real k20H2p = 1.90e-10 * pow(T2, 0.203+0.041*lnT2);
  const Real k20H2o = 2.23e-10 * pow(T2, 0.237+0.058*lnT2);
  const Real k21H2p = 2.10e-12 * pow(T2, 0.889+0.043*lnT2);
  const Real k21H2o = 3.00e-12 * pow(T2, 1.198+0.525*lnT2);
  const Real k10H2 = k10H2p*fp_ + k10H2o*fo_;
	const Real k20H2 = k20H2p*fp_ + k20H2o*fo_;
	const Real k21H2 = k21H2p*fp_ + k21H2o*fo_;
  /*e*/
  /*fit from Bell+1998*/
  const Real k10e = 5.12e-10 * pow(T, -0.075);
  const Real k20e = 4.86e-10 * pow(T, -0.026);
  const Real k21e = 1.08e-14 * pow(T, 0.926);
  /*total collisional rates*/
	const Real q10 = k10HI*nHI + k10H2*nH2 + k10e * ne;
	const Real q20 = k20HI*nHI + k20H2*nH2 + k20e * ne;
	const Real q21 = k21HI*nHI + k21H2*nH2 + k21e * ne;
	const Real q01 = (g1OI_/g0OI_) * q10 * exp( -E10OI_/(kb_*T) );
	const Real q02 = (g2OI_/g0OI_) * q20 * exp( -E20OI_/(kb_*T) );
	const Real q12 = (g2OI_/g1OI_) * q21 * exp( -E21OI_/(kb_*T) );

	return cooling3Level_(q01,q10, q02, q20, q12, q21, A10OI_,A20OI_,
												A21OI_, E10OI_, E20OI_, E21OI_, x_OI);
}

Real coolingCII(const Real x_e, const Real x_Cplus, const Real x_HI,
                const Real x_H2, const Real nH, const Real T) {
  const Real nHI = x_HI * nH;
  const Real nH2 = x_H2 * nH;
  const Real ne = x_e * nH;
	const Real q10 = q10CII_(nHI, nH2, ne, T);
	const Real q01 = (g1CII_/g0CII_) * q10 * exp( -E10CII_/(kb_*T) );
	return cooling2Level_(q01, q10, A10CII_, E10CII_, x_Cplus);
}

Real coolingCI(const Real x_e, const Real x_CI, const Real x_HI,
               const Real x_H2, const Real nH, const Real T) {
  const Real nHI = x_HI * nH;
  const Real nH2 = x_H2 * nH;
  const Real ne = x_e * nH;
	/*e collisional coefficents from Johnson, Burke, & Kingston 1987, 
	 * JPhysB, 20, 2553*/
	const Real T2 = T/100.;
	const Real lnT2 = log(T2);
	const Real lnT = log(T);
	/*ke(u,l) = fac*gamma(u,l)/g(u)*/
	const Real fac = 8.629e-8 * sqrt(1.0e4/T);
	Real k10e, k20e, k21e;
	Real lngamma10e, lngamma20e, lngamma21e; /*collisional strength*/
	if (T < 1.0e3) {
		lngamma10e = (((-6.56325e-4*lnT -1.50892e-2)*lnT + 3.61184e-1)*lnT 
				          -7.73782e-1)*lnT - 9.25141;
		lngamma20e = (((0.705277e-2*lnT - 0.111338)*lnT +0.697638)*lnT
				          - 1.30743)*lnT -7.69735;
		lngamma21e = (((2.35272e-3*lnT - 4.18166e-2)*lnT +0.358264)*lnT 
				          - 0.57443)*lnT -7.4387;

	} else {
		lngamma10e = (((1.0508e-1*lnT - 3.47620)*lnT + 4.2595e1)*lnT
									- 2.27913e2)*lnT + 4.446e2;
		lngamma20e = (((9.38138e-2*lnT - 3.03283)*lnT +3.61803e1)*lnT 
									- 1.87474e2)*lnT +3.50609e2;
		lngamma21e = (((9.78573e-2*lnT - 3.19268)*lnT +3.85049e1)*lnT 
				          - 2.02193e2)*lnT +3.86186e2;
	}
	k10e = fac * exp(lngamma10e) / g1CI_;
	k20e = fac * exp(lngamma20e) / g2CI_;
	k21e = fac * exp(lngamma21e) / g2CI_;
	/*HI collisional rates, Draine (2011) ISM book Appendix F Table F.6
	 * NOTE: this is more updated than the LAMBDA database.*/
	const Real k10HI = 1.26e-10 * pow(T2, 0.115+0.057*lnT2);
	const Real k20HI = 0.89e-10 * pow(T2, 0.228+0.046*lnT2);
	const Real k21HI = 2.64e-10 * pow(T2, 0.231+0.046*lnT2);
	/*H2 collisional rates, Draine (2011) ISM book Appendix F Table F.6*/
	const Real k10H2p = 0.67e-10 * pow(T2, -0.085+0.102*lnT2);
	const Real k10H2o = 0.71e-10 * pow(T2, -0.004+0.049*lnT2);
	const Real k20H2p = 0.86e-10 * pow(T2, -0.010+0.048*lnT2);
	const Real k20H2o = 0.69e-10 * pow(T2, 0.169+0.038*lnT2);
	const Real k21H2p = 1.75e-10 * pow(T2, 0.072+0.064*lnT2);
	const Real k21H2o = 1.48e-10 * pow(T2, 0.263+0.031*lnT2);
	const Real k10H2 = k10H2p*fp_ + k10H2o*fo_;
	const Real k20H2 = k20H2p*fp_ + k20H2o*fo_;
	const Real k21H2 = k21H2p*fp_ + k21H2o*fo_;
	/* The totol collisonal rates*/
	const Real q10 = k10HI*nHI + k10H2*nH2 + k10e*ne;
	const Real q20 = k20HI*nHI + k20H2*nH2 + k20e*ne;
	const Real q21 = k21HI*nHI + k21H2*nH2 + k21e*ne;
	const Real q01 = (g1CI_/g0CI_) * q10 * exp( -E10CI_/(kb_*T) );
	const Real q02 = (g2CI_/g0CI_) * q20 * exp( -E20CI_/(kb_*T) );
	const Real q12 = (g2CI_/g1CI_) * q21 * exp( -E21CI_/(kb_*T) );
	return cooling3Level_(q01,q10, q02, q20, q12, q21, A10CI_,A20CI_,
												A21CI_, E10CI_, E20CI_, E21CI_, x_CI);
}

Real coolingCO(const Real x_e, const Real x_CO, const Real x_HI,
               const Real x_H2, const Real nH, const Real T, 
               const Real dvdr) {
  const Real nHI = x_HI * nH;
  const Real nH2 = x_H2 * nH;
  const Real ne = x_e * nH;
  const Real nCO = x_CO * nH;
  //effective column of CO
  //maximum escape probability length, in cgs unites
  const Real Leff_CO_max = 3.086e20; //100 pc 
  const Real mCO = 4.68e-23;
  const Real vth = sqrt(2. * kb_ * T / mCO);
	const Real grad_small = vth/Leff_CO_max;
  const Real gradeff = MAX(dvdr, grad_small);
  const Real NCOeff = nCO / gradeff;
  //maximum temperature above which use Tmax for cooling rate interpolation
  const Real Tmax_CO = 2000.; 
  Real T1 = 0;;
  if (T < Tmax_CO) {
    T1 = T;
  } else {
    T1 = Tmax_CO;
  }
  const Real facT = pow(1. - exp(-T1), 1.0e3);
  /*small number for a very small NCOeff*/
  const Real eps = 1.0e13;
  const Real log_NCOeff = log10(NCOeff*1.0e5 + eps); /*unit: cm^-2 / (km/s) */
  const Real Troot4 = pow(T1, 0.25);
  const Real neff = nH2 + 1.75*Troot4 * nHI + 680.1/Troot4 * ne;
  /* interpolate parameters using given T and NCOeff*/
  /* index of T and Neff*/
  const int iT0 =  linearInterpIndex_(lenTCO_, TCO_, T1);
  const int iNeff0 = linearInterpIndex_(lenNeffCO_, NeffCO_, log_NCOeff);
  /* L0 */
  const Real log_L0 = - LP1Di_(TCO_, L0CO_, iT0, T1);
  const Real L0 = pow(10, log_L0);
  /* LLTE */
  const Real log_LLTE = - LP2Di_(TCO_, NeffCO_, lenTCO_, iT0, iNeff0, 
                                           LLTECO_, T1, log_NCOeff);
  const Real LLTE = pow(10, log_LLTE);
  /* n1/2*/
  const Real log_nhalf = LP2Di_(TCO_, NeffCO_, lenTCO_, iT0, iNeff0, 
                                          nhalfCO_, T1, log_NCOeff);
  const Real nhalf = pow(10, log_nhalf);
  /* alpha*/
  const Real alpha = LP2Di_(TCO_, NeffCO_, lenTCO_, iT0, iNeff0, 
                                      alphaCO_, T1, log_NCOeff);
  const Real inv_LCO = 1./L0 + neff/LLTE 
                         + 1./L0 * pow(neff/nhalf, alpha) * (1. - nhalf*L0/LLTE);
  return (1./inv_LCO) * neff * x_CO * facT;
}

Real coolingRec(const Real x_e, const Real nH, const Real T, 
                const Real Z_d, const Real G_PE) {
  const Real ne = x_e * nH;
  const Real x = 1.7 * G_PE * sqrt(T)/ne + 50.;
  const Real lnx = log(x);
  const Real cooling = 1.0e-28 * ne * pow(T, DPE_[0] + DPE_[1]/lnx) 
                          * exp( DPE_[2] + (DPE_[3] - DPE_[4]*lnx)*lnx );
  return cooling * Z_d;
}

Real heating(const Real x_e, const Real x_HI, const Real x_H2,
             const Real nH, const Real T, const Real Z,
             const Real xi_CR, const Real G_PE, const Real G_H2) {
  const Real h_PE = heatingPE(x_e, nH, T, Z, G_PE);
  const Real h_CR = heatingCR(x_e, x_HI, x_H2, nH, xi_CR)
  const Real h_H2pump = heatingH2pump(x_HI, x_H2, nH, T, G_H2)
  return nH * (h_PE, h_CR, h_H2pump);
}

Real cooling(const Real x_e, const Real x_HI, const Real x_H2,
             const Real x_Cplus, const Real x_CI, const Real x_CO, const Real x_OI, 
             const Real nH, const Real T, const Real dvdr,
             const Real Z, const Real G_PE) {
  const Real c_Lya = coolingLya(x_e, x_HI, nH, T);
  const Real c_OI = coolingOI(x_e, x_OI, x_HI, x_H2, nH, T);
  const Real c_CII = coolingCII(x_e, x_Cplus, x_HI, x_H2, nH, T);
  const Real c_CI = coolingCI(x_e, x_CI, x_HI, x_H2, nH, T);
  const Real c_CO = coolingCO(x_e, x_CO, x_HI, x_H2, nH, T, dvdr);
  const Real c_Rec = coolingRec(x_e, nH, T, Z, G_PE);
  return nH * (c_Lya + c_OI + c_CII + c_CI + c_CO + c_Rec);
}
