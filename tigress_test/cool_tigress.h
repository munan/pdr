#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#define Real double //TODO:replace this by include defs.in file in athena

/*----------------------------------------------------------------------------*/
/* PUBLIC FUCNTIONS                                                           */
/*----------------------------------------------------------------------------*/
void get_abundances(const Real nH, const Real T, const Real dvdr, const Real Z,
                    const Real xi_CR, const Real G_PE, const Real G_CI,
                    const Real G_CO, const Real G_H2,
                    Real *px_e, Real *px_HI, Real *px_H2, Real *px_Cplus,
                    Real *px_CI, Real *px_CO, Real *px_OI);
Real heating(const Real x_e, const Real x_HI, const Real x_H2,
             const Real nH, const Real T, const Real Z,
             const Real xi_CR, const Real G_PE, const Real G_H2);
Real cooling(const Real x_e, const Real x_HI, const Real x_H2,
             const Real x_Cplus, const Real x_CI,
             const Real x_CO, const Real x_OI, 
             const Real nH, const Real T, const Real dvdr,
             const Real Z, const Real G_PE);

// TODO: move this back to private functions
//heating---------------------------------------------------------------------
//cosmic ray heating
Real heatingCR(const Real x_e, const Real x_HI, const Real x_H2,  
              const Real nH, const Real xi_CR);
//photo electric heating on dust
Real heatingPE(const Real x_e, const Real nH, const Real T, 
               const Real Z_d, const Real G_PE);
//UV-pumping of H2
Real heatingH2pump(const Real x_HI, const Real x_H2, const Real nH,
                   const Real T, const Real G_H2);

//cooling---------------------------------------------------------------------
//HI Lyman alpha cooling
Real coolingLya(const Real x_e, const Real x_HI, const Real nH,
                const Real T);
//OI cooling
Real coolingOI(const Real x_e, const Real x_OI, const Real x_HI,
               const Real x_H2, const Real nH, const Real T);
//C+ cooling
Real coolingCII(const Real x_e, const Real x_Cplus, const Real x_HI,
                const Real x_H2, const Real nH, const Real T);
//CI cooling
Real coolingCI(const Real x_e, const Real x_CI, const Real x_HI,
               const Real x_H2, const Real nH, const Real T);
//CO rotational line cooling, dvdr in cgs units
Real coolingCO(const Real x_e, const Real x_CO, const Real x_HI,
               const Real x_H2, const Real nH, const Real T, 
               const Real dvdr);
//cooling by recombination of e on PAHs
Real coolingRec(const Real x_e, const Real nH, const Real T, 
                       const Real Z_d, const Real G_PE);
