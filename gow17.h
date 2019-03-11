/* Date: May 4, 2015, Author: Munan Gong
 * chemistry network in Gong, Ostriker and Wolifre (2017)
 */

#ifndef GOW17_H_
#define GOW17_H_

#include <stdio.h>
#include <math.h> /*a^x = pow(a,x)*/
#include <string>
#include <map>
#include <algorithm> /*std::min*/
#include <sundials/sundials_types.h> /* realtype type*/
#include <nvector/nvector_serial.h> /* N_Vector type*/
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM*/
#include "ode.h" /*base class*/
#include "sundial.h" /*Ith IJth macro and CheckFlag function*/
#include "shielding.h" /*CO self shielding*/
#include "thermo.h"
#include "radfield.h"
//#include "chemReactions.h" /*chemical reactions network*/

class gow17 : public Ode {
  friend class RadField;
  friend class CoolingFunction;
  public:
		/* Map spec_list members to integer index*/
		typedef std::map<std::string, int> SpecMap;

    gow17(); 
    ~gow17();
    int Dimen() const {return kDimen;}
    /*return the number of heating and cooling processes*/
    int GetnE() const {return nE_;}

    /* right hand side of ode */
    int RHS(const realtype t, const N_Vector y,
                   N_Vector ydot);
    /*Jacobian*/
    int Jac(const realtype t,
            const N_Vector y, const N_Vector fy, 
            SUNMatrix J, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

		/* Map of species name to index*/
		static const SpecMap spec_index_map;
		void PrintSpecMap() const;
    /*return the index of species*/
    static int id(const std::string s) {
      SpecMap::const_iterator it = gow17::spec_index_map.find(s);
      if (it == gow17::spec_index_map.end()) {
        printf("spec=%s\n", s.c_str());
        throw std::runtime_error("ERROR: spec not in map");
      } else {
        return it->second;
      }
    }

		/* Print chemistry network and current rates.*/
		void PrintChemNet(FILE *pf) const;
    void Etol(const double abstol_T);
		void SetnH(const double nH);
		void SetbCO(const double bCO);
		void SetNCO(const double NCO);
		void SetNH(const double NH);
    void SetZd(const double Zd);
    void SetZg(const double Zg);
    void SetIonRate(const double ion_rate);
    void SetRadField(double *GPE, double *Gph, double *GISRF);
    /* Set total C and O abundance */
    void SetxCtot(const double xC);
    void SetxOtot(const double xO);
    /* Set temperature to be constant in ChemInit_() calculations of rates*/
    void SetConstTemp(const double T);
		/*Set velocity gradient: 1e14 cm/s/cm = 0.3km/s/pc*/
		void SetGradv(const double gradv);
    /*Set to use global NCO for NCOeff or not*/
    void SetNCOeffGlobal(const bool isNCOeff_global);
    void Setdvdr(const bool is_dvdr);
    void SetCoolingCOThin(const bool isCoolingCOThin);
    /*Set to bCO as a function of length scale or not*/
    void SetbCOL(const bool isbCO_L);
		void SetdxCell(const double dx_cell);
		/*output physical parameters to file*/
		void WriteParam(FILE *pf) const;
		/*output rates of reactions to file*/
		void WriteRates(FILE *pf);
		/*return current time*/
		double GetTime() const;
		/*return current nH*/
		double GetnH() const;
		/*return ionization rate*/
		double GetIonRate() const;
    /* Calculate internal energy */
    static double GetE(const double T, const double xH2, const double xe);
    /* Return total C abundance*/
    double GetxCtot() const;
    /* Return total O abundance*/
    double GetxOtot() const;
    /* Return total He abundance*/
    double GetxHetot() const;
    /* Return dust metallicity*/
    double GetZd() const;
    /* Set to scale dust rain reaction*/
    void SetfH2gr(double fH2gr);
    void SetfHplusgr(double fHplusgr);
    void SetfCplusgr(double fCplusgr);
    void SetfSplusgr(double fSplusgr);
    void SetfSiplusgr(double fSiplusgr);
    void SetfHeplusgr(double fHeplusgr);
    void SetfCplusCR(double fCplusCR); /* set scale for CR ionization of CI*/
    void IsH2dissHeating(const bool isH2diss_heating);
    void IsH2rvCooling(const bool isH2rv_Cooling);
    void IsH2grHeating(const bool isH2gr_heating);

		/*copy y_ to array*/
		void CopyAbd(double *y) const;
    void CopyThermoRates(double *y) const;

  private:
    static const int kDimen = 14; 
		/*number of ghost species: the abundances of which are calculated from
		 * other species or a fixed value.*/
		static const int n_ghost_ = 7;
		static const double xHe_; /*He aboundance per H*/
    static const double xC_std_; /* total C aboundance per H, std ISM value*/
    static const double xO_std_; /* total O aboundance per H, std ISM value*/
		static const double xS_std_; /*S aboundance per H, std ISM value*/
		static const double xSi_std_; /*Si aboundance per H, std ISM value*/
    /* mass of mH and mCO */
    static const double mH_;
    static const double mCO_;
		/*Array of species*/
		static const std::string spec_list_[kDimen+n_ghost_];
		/*Use to initialize the spec_index_map*/
		static SpecMap InitMap_();

		/*-------------Physical parameters-----------------*/
		/*Cosmic ray ionization rate in (s-1 H-1). This is a bit different from
		 * NL99, see DESPOTIC NL99 implementation.*/
    double Zg_; /*gas metallicity*/
    double Zd_; /*dust abundance relative to solar neighbourhood*/
    double xC_; /* total C aboundance per H*/
    double xO_; /* total O aboundance per H*/
    double xS_; /* total S aboundance per H*/
    double xSi_; /* total Si aboundance per H*/
		double ion_rate_;
		/*H atom number density = rho_gas/mH, in cm-3*/
		double nH_;
		/*column denstiy of CO, in cm^-2*/
		double NCO_;
		double NH_;
		/*dx_cell_: when not equal to zero, take into the column density of 
		 * dx_cell_/2 * nH * xCO(xH2) into account when calculating
		 * self-sheilding*/
		double dx_cell_;
    /*Veolcity dispersion of CO for calculating NCOeff if using global NCO*/
		double bCO_;
    /*Set to use global value for NCOeff or not*/
    bool isNCOeff_global_;
    bool is_dvdr_;
    /*Set to use bCO as a function of length scale or not.*/
    bool isbCO_L_;
    /*Set whether to use optically thin CO cooling (NCOeff = 0) */
    bool isCoolingCOThin_;
		
		/*----------------------Chemical reactions----------------------------*/
		/*store index for useful species*/
		const int iCO_;
		const int iH_;
		const int iH2_;
		const int ie_;
		const int iH2plus_;
		const int iC_;
		const int iHCOplus_;
		const int iCH_;
		const int iCplus_;
		const int iO_;
		const int iOH_;
		const int iHe_;
		const int iHeplus_;
		const int iH3plus_;
		const int iHplus_;
		const int iOplus_;
		const int iS_;
		const int iSplus_;
		const int iSi_;
		const int iSiplus_;
    const int iE_; /* internal energy index */
		
		/*initialize chemistry reations.
     *Arguments:
     *y: array of species and internal energy.*/
		void ChemInit_(const double *y);
    /* Caculate change of internal energy for all heating and cooling processes.
     * Processes included are listed below.
     * Heating (5): 
     *  - LCR: cosmic ray ionization of H, He, and H2
     *  - LPE: photo electric effect on dust
     *  - LH2gr: H2 formation on dust grains
     *  - LH2pump: H2 UV pumping
     *  - LH2diss: H2 photo dissiociation
     * Cooling (10): 
     *  - GCII: C+ fine structure line
     *  - GCI: CI fine structure line
     *  - GOI: OI fine structure line
     *  - GLya: collisional exicited lyman alphya line 
     *  - GCOR: CO rotational lines
     *  - GH2: H2 vibration and rotation lines
     *  - GDust: dust thermo emission
     *  - GRec: reconbination of e on PAHs
     *  - GH2diss: collisional dissociation of H2
     *  - GHIion: collisional ionization of HI
     * Arguments:
     * y: array of species and internal energy.
     * Return:
     * Rate of internal energy changing per H in erg s^-1 H^-1.
     */
    /* calculate C+ + e RR and DR rates*/
    double CII_rec_rate_(const double temp);
    double dEdt_(const double *y, const bool is_store_rates=true);
		/*copy y and set the ghost species to yghost*/
		void SetGhostSpec_(const N_Vector y, double yghost[kDimen+n_ghost_]);
		/*Cosmic ray ionizations*/
		static const int n_cr_ = 8; /*number of reactions*/
    static const int icr_H2_;
    static const int icr_He_;
    static const int icr_H_;
		static const int incr_[n_cr_]; /*reactant*/
		static const int outcr_[n_cr_]; /*product*/
		static const double kcr_base_[n_cr_]; /*coefficents of rates relative to H*/
		/*rates of reactions, updated by ChemInit() the physical parameters*/
		double kcr_[n_cr_]; 
		/*2 body reactions*/
		static const int n_2body_ = 35;
    static const int i2body_H2_H;
    static const int i2body_H2_H2;
    static const int i2body_H_e;
		static const int in2body1_[n_2body_];
		static const int in2body2_[n_2body_];
		static const int out2body1_[n_2body_];
		static const int out2body2_[n_2body_];
		static const double k2Texp_[n_2body_]; /*exponent of temp dependance*/
		static const double k2body_base_[n_2body_]; /*base rate coefficents*/
		double k2body_[n_2body_];
    static const double A_kCHx_;
    static const double n_kCHx_;
    static const double c_kCHx_[4];
    static const double Ti_kCHx_[4];
    /*(15) H2 + *H -> 3 *H  and (16) H2 + H2 -> H2 + 2 *H*/
    /*temperature above which collisional dissociation (15), (16) and (17)
     * will be importatant k>~1e-30*/
    static const double temp_coll_;
    double k9l_;
    double k9h_;
    double k10l_;
    double k10h_;
    double ncrH_;
    double ncrH2_;
		/*photo reactions*/
		static const int n_ph_ = 7;
    static const int iph_C_;
    static const int iph_CO_;
    static const int iph_H2_;
		static const int inph_[n_ph_];
		static const int outph1_[n_ph_];
		static const double kph_base_[n_ph_]; /*base rate of photo reaction*/
		static const double kph_avfac_[n_ph_];/*exponent factor in front of AV*/
		double kph_[n_ph_];
		/*Grain assisted recombination of H and H2*/
		static const int n_gr_ = 6;
    static const int igr_H_;
		static const int ingr_[n_gr_];
		static const int outgr_[n_gr_];
		/*constants for H+, C+, He+ grain recombination.
     * See Draine ISM book table 14.9 page 158, or Weingartner+Draine2001.*/
		static const double cHp_[7]; 
		static const double cCp_[7]; 
		static const double cHep_[7]; 
		static const double cSp_[7]; 
		static const double cSip_[7]; 
		/*factor to calculate psi in H+ recombination on grain*/
		double psi_gr_fac_;
		double kgr_[n_gr_];
    /* constant temperature */
    bool const_temp_;
    double temp_;
    /* absolute value of velocity gradiant */
    double gradv_;
    /* absolute value of number density gradiant*/
    double gradnH_;
    /* heating (L) and cooling (G) processes. Details see dEdt_() comments.*/
    static const int nE_ = 15; /*number of heating and cooling processes*/
    double LCR_;
    double LPE_;
    double LH2gr_;
    double LH2pump_;
    double LH2diss_;
    double GCII_;
    double GCI_;
    double GOI_;
    double GLya_;
    double GCOR_;
    double GH2_;
    double GDust_;
    double GRec_;
    double GH2diss_;
    double GHIion_;
    /* radiation field pointer*/
    double *GPE_;
    double *GISRF_;
    double *Gph_;
    /* scale formation rate on grains.*/
    double fH2gr_;
    double fHplusgr_;
    double fCplusgr_;
    double fSplusgr_;
    double fSiplusgr_;
    double fCplusCR_;
    double fHeplusgr_;
    bool isH2diss_heating_;
    bool isH2gr_heating_;
    bool isH2rv_cooling_;

};

#endif /*GOW17_H_*/
