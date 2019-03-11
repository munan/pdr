/* Date: June, 19, 2015. Author: Munan Gong
 * Heating and cooling processes.
 */

#ifndef THERMO_H_
#define THERMO_H_
#include <math.h>
#include "interp.h"

class Thermo {
  friend class gow17;
  friend class RadField;
  friend class CoolingFunction;
	public:
		Thermo();
		~Thermo();
		/*----heating and cooling rate per H----*/
		/* Heating by cosmic ray ionization of H, He, and H2.
		 * Arguments:
		 * xs = ns/nH, abundances of species s.
		 * ks: cosmic ray ioniztion rate per s particle.
		 * nH: number density of H atom.
		 * Return: 
		 * cosmic ray heating in erg H^-1 s^-1*/
		static double HeatingCr(const double xe, const double nH,
										 const double xHI, const double xHe, const double xH2,
		   							 const double kHI, const double kHe, const double kH2);
		/* Heating by photo electric effect on dust, including collsional cooling. 
     * WD2001 Table 2 and 3. Use the diffuse ISM: Rv=3.1, ISRF, bc=4.
		 * Arguments:
		 * G: UV radiation scaled by solar neighbourhood value, including
     * extinction.  G=G0 * exp(-NH*sigmaPE_) at one line of sight.
		 * Zd: dust abundance scaled by solar neighbourhood value.
     * T: temperature in K
     * ne: electron number density in cm^-3.
		 * Return:
		 * Photo electric by dust heating rate in erg H^-1 s^-1.*/
		static double HeatingPE(const double G, const double Zd, const double T,
                            const double ne);
    /* Heating by H2 formation on dust grains.
     * From Hollenbach + McKee 1979
     * (0) *H + *H + gr -> H2 + gr
     * Arguments:
     * xi = ni/nH
     * T: temperature in K
     * kgr: grain reaction rates, kgr_ in gow17. 
     * Return:
     * Heating rate by H2 formation on dust grains in erg H^-1 s^-1. */
    static double HeatingH2gr(const double xHI, const double xH2, const double nH,
                       const double T, const double kgr);
    /* Heating by H2 UV pumping.
     * From Hollenbach + McKee 1979
     * Arguments:
     * xi = ni/nH
     * T: temperature in K
     * dot_xH2_photo = dxH2/dt by photo dissociation of H2 by UV light.
     * Calculated in RHS in gow17.
     * Return:
     * Heating rate by H2 UV pumping in erg H^-1 s^-1.*/
    static double HeatingH2pump(const double xHI, const double xH2, const double nH,
                         const double T, const double dot_xH2_photo);
    /* Heating by H2 photo dissiociation.
     * From Black + Dalgarno 1977, 0.4eV per reaction.
     * Arguments:
     * dot_xH2_photo = dxH2/dt by photo dissociation of H2 by UV light.
     * Calculated in RHS in gow17.
     * Return:
     * Heating rate by H2 photo dissiociation in erg H^-1 s^-1.*/
    static double HeatingH2diss(const double dot_xH2_photo);
		/* Cooling by C+ fine structure line.
		 * Collisional species: HI, H2, e.
		 * Arguments:
		 * xCII = nC+/nH.
		 * ni: number density of species i, in cm^-3.
		 * T: temperature in K
		 * Return:
		 * Cooling rate for C+ fine structure line in erg H^-1 s^-1*/
		static double CoolingCII(const double xCII, const double nHI, const double nH2,
										  const double ne, const double T);
		/* Cooling by CI fine structure line.
		 * Collisional species: HI, H2
     * Note: ignored H+.
     * Arguments:
		 * xCI = nCI/nH.
		 * ni: number density of species i, in cm^-3.
		 * T: temperature in K
		 * Return:
		 * Cooling rate for C fine structure line in erg H^-1 s^-1*/
		static double CoolingCI(const double xCI, const double nHI, const double nH2,
                     const double ne, const double T);
		/* Cooling by OI fine structure line.
		 * Collisional species: HI, H2, e
     * Note: the cooling rate is very insensitive to ne, for xe <~0.5.
     * At xe >~0.5 region, the O+ and other cooling will start to be important
     * anyway.
		 * Arguments:
		 * xOI = nOI/nH.
		 * ni: number density of species i, in cm^-3.
		 * T: temperature in K
		 * Return:
		 * Cooling rate for OI fine structure line in erg H^-1 s^-1*/
		static double CoolingOI(const double xOI, const double nHI, const double nH2,
										 const double ne, const double T);
		/* Cooling by collisional exicited lyman alphya line.
		 * Collisional species: e
		 * Arguments:
		 * xHI = nHI/nH.
		 * ni: number density of species i, in cm^-3.
		 * T: temperature in K
		 * Return:
		 * Cooling rate for Lyman alpha line in erg H^-1 s^-1*/
    static double CoolingLya(const double xHI, const double ne, const double T);
    /* Cooling by CO rotational lines.
     * Collision species: HI, H2, e
     * Note: from Omukai+2010
     * Auguments:
     * xCO = nCO/nH
     * ni: number density of species i, in cm^-3
     * T: temperature in K
     * NCOeff: effective column density of CO using LVG approximation. See
     * notes.
     * Return:
     * Cooling rate for CO rotational lines in erg H^-1 s^-1*/
    static double CoolingCOR(const double xCO, const double nHI, const double nH2, 
                      const double ne, const double T, const double NCOeff);
    /* Cooling by H2 vibration and rotation lines.
     * Collision species: HI, H2, He, H+, e
     * Note: Using Glover + Abel 2008 fitting formulas in Table 8, assuming
     * that ortho to para ration of H2 is 3:1.
     * Auguments:
     * xH2 = nH2 / nH
     * ni: number density of species i, in cm^-3
     * T: temperature in K
     * Return:
     * Cooling rate for H2 vibrational and rotational lines in erg H^-1 s^-1 */
    static double CoolingH2(const double xH2, const double nHI, const double nH2,
                     const double nHe, const double nHplus, const double ne,
                     const double T);
    /* Cooling by dust thermo emission.
     * Tabulated from Despotic.
     * dedt_dust = L_CMB - G_dust, thermo - PsiGD = 0
     * dedt_dust = L_CMB + L_dust, ISRF - PsiGD = 0
     * maxmium rate of above.
     * Ignored ISRF and IR heating. Assume Zd = 1.
     * Auguments:
     * Zd: dust metalicity compared to solar neighbourhood.
     * nH: hydrogen number density. Here implicitly assume all in H2, which
     * determines alpha_gd in despotic (see eq B8 in despotic paper).
     * Tg: gas temperature
     * GISRF: strength of ISRF = chi * exp(-sigma_{d, ISRF} * NH)
     * Return:
     * Cooling rate for dust in erg H^-1 s^-1 */
    static double CoolingDust(const double Zd, const double nH, const double Tg,
                              const double GISRF);
    /* Cooling of gas by dust coupling. Assume a constant dust temperature Td.
     * Auguments:
     * Zd: dust metalicity compared to solar neighbourhood.
     * nH: hydrogen number density. Here implicitly assume all in H2, which
     * Tg: gas temperature
     * Td: dust temperature
     * Return:
     * Cooling rate for dust in erg H^-1 s^-1 */
    static double CoolingDustTd(const double Zd, const double nH,  const double Tg, 
                                const double Td);
    /* Cooling by reconbination of e on PAHs.
     * From WD2001 Eq(45).
     * Arguments:
     * Zd: dust metalicity compared to solar neighbourhood.
     * T: temperature in K.
     * ne: number denstiy of electrons.
		 * G: UV radiation scaled by solar neighbourhood value, including
     * extinction. G=G0 * exp(-NH*sigmaPE_) at one line of sight.
     * Return:
     * Cooling rate for  recombination of e on PAHs in erg H^-1 s^-1 */
    static double CoolingRec(const double Zd, const double T, const double ne, 
                             const double G);
    /* Cooling by collisional dissociation of H2
     *  H2 + *H -> 3 *H 
     *  H2 + H2 -> H2 + 2 *H
     * reaction heat: 4.48 eV from Krome Paper
     * Arguments:
     * xi = ni/nH
     * k_H2_H, k_H2_H2: reaction rate cooefficients (k2body_ in gow17)
     * Return:
     * Cooling rate for H2 collisional dissociation in erg H^-1 s^-1.*/
    static double CoolingH2diss(const double xHI, const double xH2,
                         const double k_H2_H, const double k_H2_H2);
    /* Cooling by collisional ionization of HI
     *  *H + *e -> H+ + 2 *e
     * reaction heat: 13.6 eV from Krome Paper
     * Arguments:
     * xi = ni/nH
     * k_H_e: reaction rate cooefficients (k2body_ in gow17)
     * Return:
     * Cooling rate for collisional ionization of HI in erg H^-1 s^-1. */
    static double CoolingHIion(const double xHI, const double xe, const double k_H_e);
    /* specific heat, assume that H2 rotational and vibrational levels not
     * excited.
     * xH2, xe = nH2 or ne / nH
     * xHe_total = xHeI + xHeII = 0.1 for solar value.
     * Return: specific heat per H atom.*/
    static double CvCold(const double xH2, const double xHe_total, const double xe);
	private:
		static const double eV_; /*eV in erg*/
		static const double kb_; /*boltzmann constant in erg/K*/
    static const double ca_; /*speed of light * radiation constant, or
                               stephan-bolzmann constant*4 */
    static const double TCMB_; /*CMB temperature*/
		static const double o2p_;/*ratio of ortho to para H2*/
		static const double fo_; /*ortho H2 fraction*/
		static const double fp_; /*para H2 fraction*/
    static const double sigmaPE_; /*dust cross-section for 8-13.6eV photons in cm2*/
    static const double sigmaISRF_; /*dust cross-section for ISRF in cm2*/
    static const double sigmad10_; /*dust cross-section for IR at T=10K*/
    static const double alpha_GD_;
		/*-----C+ atomic data------*/
		/* Return Collisional rate for C+ atom in s^-1.
		 * Collisional species: HI, H2, e.
		 * T: gas temeperature.
		 * ni: number density of species i, in cm^-3.*/
		static double q10CII_(const double nHI, const double nH2, const double ne,
									 const double T);
		static const double A10CII_;
		static const double E10CII_;
		static const double g0CII_;
		static const double g1CII_;
		/*-----HI atomic data------*/
    static const double A10HI_;
    static const double E10HI_;
		static const double g0HI_;
		static const double g1HI_;
		/*-----CI atomic data------*/
		static const double g0CI_;
		static const double g1CI_;
		static const double g2CI_;
		static const double A10CI_;
		static const double A20CI_;
		static const double A21CI_;
		static const double E10CI_;
		static const double E20CI_;
		static const double E21CI_;
		/*-----OI atomic data------*/
		static const double g0OI_;
		static const double g1OI_;
		static const double g2OI_;
		static const double A10OI_;
		static const double A20OI_;
		static const double A21OI_;
		static const double E10OI_;
		static const double E20OI_;
		static const double E21OI_;
    /*-----CO cooling table data, from Omukai+2010-----*/
    static const int lenTCO_ = 11;
    static const int lenNeffCO_ = 11;
    static const double TCO_[lenTCO_];
    static const double NeffCO_[lenNeffCO_];
    static const double L0CO_[lenTCO_];
    static const double LLTECO_[lenNeffCO_*lenTCO_];
    static const double nhalfCO_[lenNeffCO_*lenTCO_];
    static const double alphaCO_[lenNeffCO_*lenTCO_];
    /*-----PE heating coefficients from WD2001 Table 2, second last line ---*/
    static const double CPE_[7];
    /*-----Collisional cooling included in PE heating, WD2001 Table3, second
     * last line---*/
    static const double DPE_[5];
    /*------dust cooling table from Despotic--------*/
    /* dedt_dust = L_CMB - L_dust, thermo - PsiGD = 0
     * Ignored ISRF heating.
     * Tabulate PsiGD as a function of Tg and nH.*/
    static const int lenTg_ = 10;
    static const int lennH_ = 15;
    static const double logTg_[lenTg_];
    static const double lognH_[lennH_];
    static const double logps_[lennH_ * lenTg_];
		/*line cooling rate per H for 2 level atom. Ignore radiative excitation
		 * and de-excitation, and assume optically thin.
		 * Arguments:
		 * q01 = \sum (nc * k_{s, 01}). Collisional excitation rate per second 
		 * for all the collider species sumed up together from level 0 to 1.
		 * q10 = \sum (nc * k_{s, 10}). Collisional de-excitation rate per second,
		 * similar to q01.
		 * A10: Enstein A coefficent for spontanious emission from level 1 to 0,
		 * in sec^-1.
		 * E10: energy difference of levels E1 - E0, in erg.
		 * xs = ns/nH, abundances of species s.
		 * Return:
		 * Line cooling rate in erg H^-1 s^-1.*/
		static double Cooling2Level_(const double q01, const double q10,
											    const double A10, const double E10,
													const double xs);
		/*line cooling rate per H for 3 level atom. Ignore radiative excitation
		 * and de-excitation, and assume optically thin.
		 * Arguments:
		 * qij = \sum (nc * k_{s, ij}). Collisional excitation rate per second 
		 * for all the collider species sumed up together from level i to j.
		 * Aij: Enstein A coefficent for spontanious emission from level i to j,
		 * in sec^-1, i > j.
		 * Eij: energy difference of levels Ei - Ej, in erg.
		 * xs = ns/nH, abundances of species s.
		 * Return:
		 * Total line cooling rate in erg H^-1 s^-1.*/
		static double Cooling3Level_(const double q01, const double q10,
													const double q02, const double q20,
													const double q12, const double q21,
													const double A10, const double A20,
													const double A21, const double E10,
													const double E20, const double E21,
													const double xs);
};
#endif /*THERMO_H_*/
