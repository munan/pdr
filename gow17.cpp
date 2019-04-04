/* Date: May 4, 2015 Author: Munan Gong
 * Implimentation of the gow17 class
 */
#include "gow17.h"

/*------------Initialize static members of the class --------*/
/*species list and map*/
const std::string gow17::spec_list_[kDimen+n_ghost_] = 
		{"He+", "OHx", "CHx", "CO", "C+", "HCO+", "H2", "H+", "H3+",
		/*below are ghost species. The aboundances of ghost species are
		 * recalculated in RHS everytime by other species. */
		  "H2+", "S+", "Si+", "O+", "E", "*Si", "*S", "*C", "*O", "*He", "*e", "*H"}; 

const gow17::SpecMap gow17::spec_index_map = InitMap_();

gow17::SpecMap gow17::InitMap_() {
	gow17::SpecMap spec_index_map;
	for (int i=0; i<kDimen+n_ghost_; i++){
		spec_index_map[spec_list_[i]] = i;
	}
	return spec_index_map;
}

/*network specific constants*/
const double gow17::xHe_ = 0.1; 
const double gow17::xC_std_ = 1.6e-4;
const double gow17::xO_std_ = 3.2e-4; 
const double gow17::xS_std_ = 3.5e-6; 
const double gow17::xSi_std_ = 1.7e-6; 
const double gow17::temp_coll_ = 7.0e2;
const double gow17::mH_ = 1.67e-24;
const double gow17::mCO_ = 4.68e-23;

/*chemical network*/
/*cosmic ray chemistry network*/
/* (0) cr + H2 -> H2+ + *e
 * (1) cr + *He -> He+ + *e 
 * (2) cr + *H  -> H+ + *e 
 * -----added as Clark + Glover 2015---- 
 * (3) cr + *C -> C+ + *e     --including direct and cr induce photo reactions 
 * (4) crphoto + CO -> *O + *C       
 * (5) cr + CO -> HCO+ + e  --schematic for cr + CO -> CO+ + e
 * -----S, CR induced photo ionization, experimenting----
 * (6) cr + S -> S+ + e, simply use 2 times rate of C, as in UMIST12
 * -----Si, CR induced photo ionization, experimenting----
 * (7) cr + Si -> Si+ + e, UMIST12
 */
const int gow17::icr_H2_ = 0;
const int gow17::icr_He_ = 1;
const int gow17::icr_H_ = 2;
const int gow17::incr_[n_cr_] = {id("H2"), id("*He"), id("*H"), 
                                 id("*C"), id("CO"), id("CO"),
                                 id("*S"), id("*Si")};
const int gow17::outcr_[n_cr_] = {id("H2+"), id("He+"), id("H+"), 
                                  id("C+"), id("*O"), id("HCO+"),
                                  id("S+"), id("Si+")};
const double gow17::kcr_base_[n_cr_] = 
{2.0, 1.1, 1.0, 
560., 90., 6.52,
2040., 8400.}; 


/*2 body reactions*/
/*NOTE: photons from recombination are ignored*/
/* Reactions are, in order.
 * -- are equations of special rate treatment in Glover, Federrath+ 2010:
 (0) H3+ + *C -> CH + H2         --Vissapragada2016 new rates
 (1) H3+ + *O -> OH + H2        
 (2) H3+ + CO -> HCO+ + H2
 (3) He+ + H2 -> H+ + *He + *H   --fit to Schauer1989
 (4) He+ + CO -> C+ + *O + *He   
 (5) C+ + H2 -> CH + *H         -- schematic reaction for C+ + H2 -> CH2+
 (6) C+ + OH -> HCO+             -- Schematic equation for C+ + OH -> CO+ + H.
 Use rates in KIDA website.
 (7) CH + *O -> CO + *H
 (8) OH + *C -> CO + *H          --exp(0.108/T)
 (9) He+ + *e -> *He             --(17) Case B
 (10) H3+ + *e -> H2 + *H
 (11) C+ + *e -> *C              -- Include RR and DR, Badnell2003, 2006.
 (12) HCO+ + *e -> CO + *H
 ----added in GO2012--------
 (13) H2+ + H2 -> H3+ + *H       --(54) exp(-T/46600)
 (14) H+ + *e -> *H              --(12) Case B
 ---collisional dissociation, only important at high temperature T>1e3---
 (15) H2 + *H -> 3 *H            --(9) Density dependent. See Glover+MacLow2007
 (16) H2 + H2 -> H2 + 2 *H       --(10) Density dependent. See Glover+MacLow2007
 (17) *H + *e -> H+ + 2 *e       --(11) Relates to Te
 ----added for H3+ destruction in addition to (10)----
 (18) H3+ + *e -> *3H            --(111)
 ----added He+ destruction in addtion to (3), from UMIST12----
 (19) He+ + H2 -> H2+ + *He
 ----added CH reaction to match for abundances of CH---
 (20) CH + *H -> H2 + *C         
 ----added to match the Meudon code ---
 (21) OH + *O -> *O + *O + *H
 ---branching of C+ + H2 ------
 (22) C+ + H2 + *e -> *C + *H + *H
 ---S , rate from UMIST12---
 (23) S+ + *e -> *S
 (24) C+ + *S -> S+ + *C
 ---Si , rate from UMIST12---
 (25) Si+ + *e -> *Si
 (26) C+ + *Si -> Si+ + *C -- not included
 --- H2O+ + e reaction ---
 (27) H3+ + *O + *e -> H2 + *O + *H
 --- OH destruction with He+
 (28) He+ + OH -> O+ + *He + *H 
 --- H2+ charge exchange with H ---
 (29) H2+ + *H -> H+ + H2 
 --- O+ reactions ---
 (30) H+ + *O -> O+ + *H -- fitting in Stancil et al. 1999, exp(-227/T)
 (31) O+ + *H -> H+ + *O -- fitting in Stancil et al. 1999
 (32) O+ + H2 -> OH + *H     -- branching of H2O+
 (33) O+ + H2 -> *O + *H + *H  -- branching of H2O+
 --- CH+ from hot gas --
 (34) C+ + H2 -> CHx + *H
 */
const int gow17::i2body_H2_H = 15;
const int gow17::i2body_H2_H2 = 16;
const int gow17::i2body_H_e = 17;
const int gow17::in2body1_[n_2body_] = 
          {id("H3+"), id("H3+"), id("H3+"), id("He+"), id("He+"),    
           id("C+"), id("C+"), id("CHx"), id("OHx"), id("He+"),
           id("H3+"), id("C+"), id("HCO+"),
					 id("H2+"), id("H+"), id("H2"), id("H2"), id("*H"),
					 id("H3+"), id("He+"), 
           id("CHx"), id("OHx"), id("C+"), id("S+"), id("C+"),
           id("Si+"), id("C+"), id("H3+"), id("He+"), id("H2+"),
           id("H+"), id("O+"), id("O+"), id("O+"), id("C+")};
const int gow17::in2body2_[n_2body_] = 
          {id("*C"), id("*O"), id("CO"), id("H2"), id("CO"),   
           id("H2"), id("OHx"), id("*O"), id("*C"), id("*e"),   
           id("*e"), id("*e"), id("*e"),
					 id("H2"), id("*e"), id("*H"), id("H2"), id("*e"),
					 id("*e"), id("H2"), 
           id("*H"), id("*O"), id("H2"), id("*e"), id("*S"),
           id("*e"), id("*Si"), id("*O"), id("OHx"), id("*H"),
           id("*O"), id("*H"), id("H2"), id("H2"), id("H2")};
/*Note: output to ghost species doesn't matter. The abundances of ghost species
 * are updated using the other species at every timestep*/
const int gow17::out2body1_[n_2body_] = 
          {id("CHx"), id("OHx"), id("HCO+"), id("H+"), id("C+"),   
           id("CHx"), id("HCO+"), id("CO"), id("CO"), id("*He"),   
           id("H2"), id("*C"), id("CO"),
					 id("H3+"), id("*H"), id("*H"), id("H2"), id("H+"),
					 id("*H"), id("H2+"), 
           id("H2"), id("*O"), id("*C"), id("*S"), id("S+"),
           id("*Si"), id("Si+"), id("H2"), id("O+"), id("H+"),
           id("O+"), id("H+"), id("OHx"), id("*O"), id("CHx")};
const int gow17::out2body2_[n_2body_] = 
          {id("H2"), id("H2"), id("H2"), id("*He"), id("*O"),   
           id("*H"), id("*H"), id("*H"), id("*H"), id("*H"),   
           id("*H"), id("*H"), id("*H"),
					 id("*H"), id("*H"), id("*H"), id("*H"), id("*e"),
					 id("*H"), id("*He"), 
           id("*C"), id("*H"), id("*H"), id("*H"), id("*C"),
           id("*H"), id("*C"), id("*O"), id("*He"), id("H2"),
           id("*H"), id("*O"), id("*H"), id("*H"), id("*H")};
const double gow17::k2Texp_[n_2body_] = 
 {0.0, -0.190, 0.0, 0.0, 0.0, 
  -1.3, 0.0, 0.0, -0.339, -0.5, 
  -0.52, 0.0, -0.64,
  0.042, 0.0, 0.0, 0.0, 0.0,
  -0.52, 0.0,
  0.26, 0.0, -1.3, -0.59, 0.0,
  -0.62, 0.0, -0.190, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0};
const double gow17::k2body_base_[n_2body_] = 
                {1.00, 1.99e-9, 1.7e-9, 1.26e-13, 1.6e-9, 
                 3.3e-13 * 0.7, 1.00, 7.0e-11, 7.95e-10, 1.0e-11, 
                 4.54e-7, 1.00, 1.06e-5,
								 1.76e-9, 2.753e-14, 1.00, 1.00, 1.00,
								 8.46e-7, 7.20e-15, 
                 2.81e-11, 3.5e-11, 3.3e-13 * 0.3, 1.6e-10, 5e-11,
                 1.46e-10, 2.1e-9, 1.99e-9, 1.00, 6.4e-10,
                 1.00, 1.00, 1.6e-9, 1.6e-9, 3.8e-14};
//rates for H3+ + C forming CH+ and CH2+
const double gow17::A_kCHx_ = 1.04e-9;
const double gow17::n_kCHx_ = 2.31e-3;
const double gow17::c_kCHx_[4] = {3.4e-8, 6.97e-9, 1.31e-7, 1.51e-4};
const double gow17::Ti_kCHx_[4] = {7.62, 1.38, 2.66e1, 8.11e3};

/* photo reactions.
 * Reaction rates in Drain 1978 field units.
 * Reactions are, in order:
 (0) h nu + *C -> C+ + *e
 (1) h nu + CH -> *C + *H
 (2) h nu + CO -> *C + *O            --self-shielding and shielding by H2
 (3) h nu + OH -> *O + *H
 ----added in GO2012--------
 (4) h nu + H2 -> *H + *H            --self- and dust shielding
 ----S, from UMIST12
 (5) h nu + *S -> S+
 ----Si, from UMIST12
 (6) h nu + *Si -> Si+
 */
const int gow17::iph_C_ = 0;
const int gow17::iph_CO_ = 2;
const int gow17::iph_H2_ = 4;
const int gow17::inph_[n_ph_] = {
              id("*C"), id("CHx"), id("CO"),
              id("OHx"), id("H2"), id("*S"), id("*Si")};
const int gow17::outph1_[n_ph_] = {
              id("C+"), id("*C"), id("*C"),
              id("*O"), id("*H"), id("S+"), id("Si+")};
const double gow17::kph_base_[n_ph_] = {3.5e-10, 9.1e-10, 2.4e-10,   
																			  3.8e-10, 5.7e-11, 
                                        6e-10, 4.5e-9}; 
const double gow17::kph_avfac_[n_ph_] = {3.76, 2.12, 3.88,  
	                                       2.66, 4.18,
                                         3.10, 2.61};

/* Grain assisted recombination of H, H2, C+ and H+
 (0) *H + *H + gr -> H2 + gr
 (1) H+ + *e + gr -> *H + gr
 (2) C+ + *e + gr -> *C + gr
 (3) He+ + *e + gr -> *He + gr
 ------S, from WD2001-----
 (4) S+ + *e + gr -> *S + gr
 ------Si, from WD2001-----
 (5) Si+ + *e + gr -> *Si + gr
 */
const int gow17::igr_H_ = 0;
const int gow17::ingr_[n_gr_] = {id("*H"), id("H+"), id("C+"), id("He+"), 
                                 id("S+"), id("Si+")};
const int gow17::outgr_[n_gr_] = {id("H2"), id("*H"), id("*C"), id("*He"), 
                                 id("*S"), id("*Si")};
const double gow17::cHp_[7] = {12.25, 8.074e-6, 1.378, 5.087e2, 1.586e-2,
															 0.4723, 1.102e-5}; 
const double gow17::cCp_[7] = {45.58, 6.089e-3, 1.128, 4.331e2, 4.845e-2,
                               0.8120, 1.333e-4};
const double gow17::cHep_[7] = {5.572, 3.185e-7, 1.512, 5.115e3, 3.903e-7,
                                0.4956, 5.494e-7};
const double gow17::cSp_[7] = {3.064, 7.769e-5, 1.319, 1.087e2, 3.475e-1,
                               0.4790, 4.689e-2};
const double gow17::cSip_[7] = {2.166, 5.678e-8, 1.874, 4.375e4, 1.635e-6,
                               0.8964, 7.538e-5};

/*------------------Constructor and Destructor---------------------*/
gow17::gow17()
  :Ode(),
	 iCO_(id("CO")),
	 iH_(id("*H")),
	 iH2_(id("H2")),
	 ie_(id("*e")),
   iH2plus_(id("H2+")),
	 iC_(id("*C")),
	 iHCOplus_(id("HCO+")),
	 iCH_(id("CHx")),
	 iCplus_(id("C+")),
	 iO_(id("*O")),
	 iOH_(id("OHx")),
	 iHe_(id("*He")),
	 iHeplus_(id("He+")),
	 iH3plus_(id("H3+")),
	 iHplus_(id("H+")),
	 iOplus_(id("O+")),
	 iS_(id("*S")),
	 iSplus_(id("S+")),
	 iSi_(id("*Si")),
	 iSiplus_(id("Si+")),
   iE_(id("E")),
   const_temp_(false),
   isH2diss_heating_(true),
   isH2gr_heating_(true),
   isH2rv_cooling_(true),
   temp_(0.),
   LCR_(0.),
   LPE_(0.),
   LH2gr_(0.),
   LH2pump_(0.),
   LH2diss_(0.),
   GCII_(0.),
   GCI_(0.),
   GOI_(0.),
   GLya_(0.),
   GCOR_(0.),
   GH2_(0.),
   GDust_(0.),
   GRec_(0.),
   GH2diss_(0.),
   GHIion_(0.),
   GPE_(NULL),
   GISRF_(NULL),
   Gph_(NULL),
   fH2gr_(1.),
   fHplusgr_(1.),
   fCplusgr_(1.),
   fSplusgr_(1.),
   fSiplusgr_(1.),
   fCplusCR_(1.),
   fHeplusgr_(1.)
{
  y_ = N_VNew_Serial(kDimen); /* Create serial vector y*/
  CheckFlag((void *)y_, "N_VNew_Serial", 0);

	/*Default values of physical parameters*/
  Zg_ = 1.;
  Zd_ = 1.;
  xC_ = Zg_ * xC_std_;
  xO_ = Zg_ * xO_std_;
  xS_ = Zg_ * xS_std_;
  xSi_ = Zg_ * xSi_std_;
	ion_rate_ = 2.0e-16;
	nH_ = 100;
	bCO_ = 0.3e5;
	/*default negative - use the value in its own cell*/
	NCO_ = 0.;
	NH_ = 0.;
  isNCOeff_global_ = false;
  is_dvdr_ = false;
  isbCO_L_ = false;
  isCoolingCOThin_ = false;
	dx_cell_ = 0;
  gradv_ = 1.0e-14; /*0.3km/s /pc*/
  gradnH_ = 3.0e-17; /* 100 cm^-3 / pc */
}

gow17::~gow17() {
  /* destroy serial vector y*/
  N_VDestroy_Serial(y_); 
} 

/*-------------------------RHS and Jac-----------------------------*/

int gow17::RHS(const realtype t, const N_Vector y, N_Vector ydot)
{
	double rate;
	/*ydot including the ghost species*/
	double ydot_[kDimen+n_ghost_];
	/*store previous y*/
	double yprev[kDimen+n_ghost_];

	/* copy y to yprev and set ghost species*/
	SetGhostSpec_(y, yprev);
  ChemInit_(yprev);

	/*set the initial ydot to zero*/
	for(int i=0; i<kDimen+n_ghost_; i++) {
		ydot_[i] = 0.0;
	}

	/*cosmic ray reactions*/
	for (int i=0; i<n_cr_; i++) {
		rate = kcr_[i] * yprev[incr_[i]];
		ydot_[incr_[i]] -= rate;
		ydot_[outcr_[i]] += rate;
	}

	for (int i=0; i<n_2body_; i++) {
		rate =  k2body_[i] * yprev[in2body1_[i]] * yprev[in2body2_[i]];
		ydot_[in2body1_[i]] -= rate;
		ydot_[in2body2_[i]] -= rate;
		ydot_[out2body1_[i]] += rate;
		ydot_[out2body2_[i]] += rate;
	}

	/*photo reactions*/
	for (int i=0; i<n_ph_; i++) {
		rate = kph_[i] * yprev[inph_[i]];
		ydot_[inph_[i]] -= rate;
		ydot_[outph1_[i]] += rate;
	}

  /*energy equation*/
  if (!const_temp_) {
    ydot_[iE_] = dEdt_(yprev);
  }

	/*grain assisted reactions*/
	for (int i=0; i<n_gr_; i++) {
		rate = kgr_[i] * yprev[ingr_[i]];
		ydot_[ingr_[i]] -= rate;
		ydot_[outgr_[i]] += rate;
	}

	/*set ydot to return*/
	for (int i=0; i<kDimen; i++) {
		NV_Ith_S(ydot, i) = ydot_[i];
	}
  return 0;
}

int gow17::Jac(const realtype t,
               const N_Vector y, const N_Vector fy, 
               SUNMatrix J, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
	double rate_pa = 0; /* rate for partial derivative respect to species a*/
	double rate_pb = 0;
	int ia, ib, ic, id;/*index for species a, b, c, d*/
	/*store previous y*/
	double yprev[kDimen+n_ghost_];
	/*Jacobian include ghost indexes*/
	double J_[kDimen+n_ghost_][kDimen+n_ghost_];

	/* copy y to yprev and set ghost species*/
	SetGhostSpec_(y, yprev);
  /* TODO: We can might skip this, which was caluclated in RHS*/
  ChemInit_(yprev);

	/*initialize J_ to be zero*/
	for (int i=0; i<kDimen+n_ghost_; i++) {
		for (int j=0; j<kDimen+n_ghost_; j++) {
			J_[i][j] = 0;
		}
	}

	/* 2 body reactions: a+b -> c+d*/
	for (int i=0; i<n_2body_; i++) {
		ia = in2body1_[i];
		ib = in2body2_[i];
		ic = out2body1_[i];
		id = out2body2_[i];
		rate_pa = k2body_[i] * yprev[ib];
		rate_pb = k2body_[i] * yprev[ia];
		J_[ia][ia] -= rate_pa;
		J_[ib][ia] -= rate_pa;
		J_[ic][ia] += rate_pa;
		J_[id][ia] += rate_pa;
		J_[ia][ib] -= rate_pb;
		J_[ib][ib] -= rate_pb;
		J_[ic][ib] += rate_pb;
		J_[id][ib] += rate_pb;
	}
	/* photo reactions a + photon -> c+d*/
	for (int i=0; i<n_ph_; i++) {
		ia = inph_[i];
		ic = outph1_[i];
		rate_pa = kph_[i];
		J_[ia][ia] -= rate_pa;
		J_[ic][ia] += rate_pa;
	}
	/*Cosmic ray reactions a + cr -> c */
	for (int i=0; i<n_cr_; i++) {
		ia = incr_[i];
		ic = outcr_[i];
		rate_pa = kcr_[i];
		J_[ia][ia] -= rate_pa;
		J_[ic][ia] += rate_pa;
	}

	/*grain reactions a + gr -> c*/
	for (int i=0; i<n_gr_; i++) {
		ia = ingr_[i];
		ic = outgr_[i];
		rate_pa = kgr_[i];
		J_[ia][ia] -= rate_pa;
		J_[ic][ia] += rate_pa;
	}

	/*copy J to return*/
	for (int i=0; i<kDimen; i++) {
		for (int j=0; j<kDimen; j++) {
			SM_ELEMENT_D(J, i, j) = J_[i][j];
		}
	}
  return 0;
}

void gow17::PrintSpecMap() const {
	SpecMap::const_iterator iter;
	for (iter=spec_index_map.begin(); iter!=spec_index_map.end(); iter++) {
		/*put () around ghost species*/
		if (iter->second >= kDimen) {
			printf("(%7s      %d)\n", (iter->first).c_str(), iter->second);
		} else {
			printf("%8s      %d\n", (iter->first).c_str(), iter->second);
		}
	}
	return;
}

/*--------------------Chemistry reactions-----------------------*/
void gow17::PrintChemNet(FILE *pf) const {
	for (int i=0; i<n_cr_; i++) {
		fprintf(pf, "cr  + %4s -> %4s,     kcr = %.2e ksi s-1 H-1\n", 
		 spec_list_[incr_[i]].c_str(), spec_list_[outcr_[i]].c_str(), kcr_base_[i]);
	}
	for (int i=0; i<n_2body_; i++) {
		fprintf(pf, "%4s  + %4s -> %4s  + %4s,     k2body = %.1e T^%.2f cm3 s-1 H-1\n", 
		 spec_list_[in2body1_[i]].c_str(), spec_list_[in2body2_[i]].c_str(),
		 spec_list_[out2body1_[i]].c_str(), spec_list_[out2body2_[i]].c_str(),
		 k2body_base_[i], k2Texp_[i]);
	}
	for (int i=0; i<n_ph_; i++) {
		fprintf(pf, "h nu  + %4s -> %4s,     kph = %.1e G0 exp(-%.1f Av) s-1 H-1\n", 
		 spec_list_[inph_[i]].c_str(), spec_list_[outph1_[i]].c_str(),
		 kph_base_[i], kph_avfac_[i]);
	}
	for (int i=0; i<n_gr_; i++) {
		fprintf(pf, "gr  + %4s -> %4s \n", 
		 spec_list_[ingr_[i]].c_str(), spec_list_[outgr_[i]].c_str());
	}
}

void gow17::ChemInit_(const double *y) {
  double T;
  if (const_temp_) {
    T = temp_;
  } else {
    T = y[iE_] / Thermo::CvCold(y[iH2_], xHe_, y[ie_]);
    //printf("Calcuate T from E.\n");
  }
	const double logT = log10(T);
	const double logT4 = log10(T/1.0e4);
	const double lnTe = log(T * 8.6173e-5);
  double ncr, n2ncr;
	double psi; /*H+ grain recombination parameter*/
  double kcr_H_fac;//ratio of total rate to primary rate
  double t1_CHx, t2_CHx;

	/*cosmic ray reactions*/
	for (int i=0; i<n_cr_; i++) {
		kcr_[i] = kcr_base_[i] * ion_rate_;
	}
  /*cosmic ray induced photo-reactions, proportional to x(H2)*/
  /*(3) cr + *C -> C+ + *e     --including direct and cr induce photo reactions 
   *(4) crphoto + CO -> *O + *C 
   *(6) cr + S -> S+ + e, simply use 2 times rate of C, as in UMIST12
   *(7) cr + Si -> Si+ + e, UMIST12 */
 /* (0) cr + H2 -> H2+ + *e
 * (1) cr + *He -> He+ + *e 
 * (2) cr + *H  -> H+ + *e */
  kcr_H_fac = 1.15 * 2*y[iH2_] + 1.5 * y[iH_];
  kcr_[0] *= kcr_H_fac;
  kcr_[2] *= kcr_H_fac;
  kcr_[3] *= (2*y[iH2_] + 3.85/kcr_base_[3]) * fCplusCR_;
  kcr_[4] *= 2*y[iH2_];
  kcr_[6] *= 2*y[iH2_];
  kcr_[7] *= 2*y[iH2_];

	/*2 body reactions*/
	for (int i=0; i<n_2body_; i++){
		k2body_[i] = k2body_base_[i] * pow(T, k2Texp_[i]) * nH_;
	}
	/*Special treatment of rates for some equations*/
  /*(0) H3+ + *C -> CH + H2         --Vissapragada2016 new rates*/
  t1_CHx = A_kCHx_ * pow( 300./T, n_kCHx_);
  t2_CHx = c_kCHx_[0] * exp(-Ti_kCHx_[0]/T) + c_kCHx_[1] * exp(-Ti_kCHx_[1]/T)
           + c_kCHx_[2]*exp(-Ti_kCHx_[2]/T) + c_kCHx_[3] *exp(-Ti_kCHx_[3]/T);
  k2body_[0] *= t1_CHx + pow(T, -1.5) * t2_CHx;
	/*(3) He+ + H2 -> H+ + *He + *H   --fit to Schauer1989 */
	k2body_[3] *= exp(-22.5/T);
  /*(5) C+ + H2 -> CH + *H         -- schematic reaction for C+ + H2 -> CH2+*/
  k2body_[5] *= exp(-23./T);
  /* ---branching of C+ + H2 ------
 (22) C+ + H2 + *e -> *C + *H + *H*/
  k2body_[22] *= exp(-23./T);
  /* (6) C+ + OH -> HCO+             -- Schematic equation for C+ + OH -> CO+ + H.
     Use rates in KIDA website.*/
  k2body_[6] = 9.15e-10 * ( 0.62 + 45.41 / sqrt(T) ) * nH_;
  //k2body_[6] = 2.9e-9 * pow( T/300.0, -0.33 ) * nH_; //Mark's rate
  /*(8) OH + *C -> CO + *H          --exp(0.108/T)*/
  k2body_[8] *= exp(0.108/T);
	/*(9) He+ + *e -> *He             --(17) Case B */
	k2body_[9] *= 11.19 + (-1.676 + (-0.2852 + 0.04433*logT) * logT )* logT;
  /* (11) C+ + *e -> *C              -- Include RR and DR, Badnell2003, 2006. */
  k2body_[11] = CII_rec_rate_(T) * nH_;
  /* (13) H2+ + H2 -> H3+ + *H       --(54) exp(-T/46600) */
  k2body_[13] *= exp(- T/46600.);
	/* (14) H+ + *e -> *H              --(12) Case B */
	k2body_[14] *= pow( 315614.0 / T, 1.5) 
									 * pow(  1.0 + pow( 115188.0 / T, 0.407) , -2.242 );
  //k2body_[14] = 3.5e-12 * pow( 300./T, 0.75 ) * nH_;
  /* (28) He+ + OH -> *H + *He + *O(O+)*/
  k2body_[28] = 1.35e-9 *( 0.62 + 0.4767 * 5.5 * sqrt(300./T) ) * nH_;
  /*--- H2O+ + e branching--
  (1) H3+ + *O -> OH + H2        
  (27) H3+ + *O + *e -> H2 + *O + *H     */
  const double h2oplus_ratio = 
                6e-10 * y[iH2_] / ( 5.3e-6 / sqrt(T) * y[ie_] );
  const double fac_H2Oplus_H2 = h2oplus_ratio / (h2oplus_ratio + 1.);
  const double fac_H2Oplus_e = 1. / (h2oplus_ratio + 1.);
  k2body_[1] *= fac_H2Oplus_H2;
  k2body_[27] *= fac_H2Oplus_e;
  /*--- O+ reactions ---
    (30) H+ + *O -> O+ + *H -- exp(-227/T)
    (31) O+ + *H -> H+ + *O 
    (32) O+ + H2 -> OH + *H     -- branching of H2O+
    (33) O+ + H2 -> *O + *H + *H  -- branching of H2O+ */
  k2body_[30] *= ( 1.1e-11 * pow(T, 0.517) + 4.0e-10 * pow(T, 6.69e-3) )*exp(-227./T);
  k2body_[31] *= 4.99e-11* pow(T, 0.405) + 7.5e-10 * pow(T, -0.458);
  k2body_[32] *= fac_H2Oplus_H2;
  k2body_[33] *= fac_H2Oplus_e;
  // (26) C+ + *Si -> Si+ + *C -- not included
  k2body_[26] = 0;

  /*Collisional dissociation, k>~1.0e-30 at T>~5e2.*/
  if (T > temp_coll_) {
    /*(15) H2 + *H -> 3 *H   
     *(16) H2 + H2 -> H2 + 2 *H
     * --(9) Density dependent. See Glover+MacLow2007*/
  	k9l_ = 6.67e-12 * sqrt(T) * exp(-(1. + 63590./T)); 
    k9h_ = 3.52e-9 * exp(-43900.0 / T);
    k10l_ = 5.996e-30 * pow(T, 4.1881) / pow((1.0 + 6.761e-6 * T), 5.6881)  
            * exp(-54657.4 / T);
    k10h_ = 1.3e-9 * exp(-53300.0 / T); 
    ncrH_ = pow(10, (3.0 - 0.416 * logT4 - 0.327 * logT4*logT4));
    ncrH2_ = pow(10, (4.845 - 1.3 * logT4 + 1.62 * logT4*logT4));
    ncr = 1. / ( y[iH_]/ncrH_ + y[iH2_]/ncrH2_ );
    n2ncr = nH_ / ncr;
    k2body_[15] = pow(10, log10(k9h_) *  n2ncr/(1. + n2ncr) 
                         + log10(k9l_) / (1. + n2ncr)) * nH_;
    k2body_[16] = pow(10, log10(k10h_) *  n2ncr/(1. + n2ncr) 
                         + log10(k10l_) / (1. + n2ncr)) * nH_;
    /* (17) *H + *e -> H+ + 2 *e       --(11) Relates to Te */
    k2body_[17] *= exp( -3.271396786e1 + 
                      (1.35365560e1 + (- 5.73932875 + (1.56315498 
                    + (- 2.877056e-1 + (3.48255977e-2 + (- 2.63197617e-3
                    + (1.11954395e-4 + (-2.03914985e-6)
        *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe
                     );
  } else {
    k2body_[15] = 0.;
    k2body_[16] = 0.;
    k2body_[17] = 0.;
  }

	/*photo reactions*/
	for (int i=0; i<n_ph_; i++) {
    kph_[i] = kph_base_[i] * Gph_[i];
	}

	/* Grain assisted recombination of H and H2*/
	/*	 (0) *H + *H + gr -> H2 + gr , constant rate from Wolfire2008 and
   *	 Hollenbach 2012*/
	kgr_[0] = 3.0e-17 * nH_ * Zd_;
	/*	 (1) H+ + *e + gr -> *H + gr
   *	 (2) C+ + *e + gr -> *C + gr
   *   (3) He+ + *e + gr -> *He + gr
   *   (4) S+ + *e + gr -> *S + gr
   *   (5) Si+ + *e + gr -> *Si + gr
   *   , rate dependent on e aboundance. */
	psi_gr_fac_ = 1.7 * (*GPE_) * sqrt(T) / nH_; 
	psi = psi_gr_fac_ / y[ie_];
	kgr_[1] = 1.0e-14 * cHp_[0] / 
		           (
			           1.0 + cHp_[1]*pow(psi, cHp_[2]) * 
								   (1.0 + cHp_[3] * pow(T, cHp_[4])
										             *pow( psi, -cHp_[5]-cHp_[6]*log(T) ) 
									 ) 
								) * nH_ * Zd_ * fHplusgr_;
	kgr_[2] = 1.0e-14 * cCp_[0] / 
		           (
			           1.0 + cCp_[1]*pow(psi, cCp_[2]) * 
								   (1.0 + cCp_[3] * pow(T, cCp_[4])
										             *pow( psi, -cCp_[5]-cCp_[6]*log(T) ) 
									 ) 
								) * nH_ * Zd_ * fCplusgr_;
	kgr_[3] = 1.0e-14 * cHep_[0] / 
		           (
			           1.0 + cHep_[1]*pow(psi, cHep_[2]) * 
								   (1.0 + cHep_[3] * pow(T, cHep_[4])
										             *pow( psi, -cHep_[5]-cHep_[6]*log(T) ) 
									 ) 
								) * nH_ * Zd_ * fHeplusgr_;
	kgr_[4] = 1.0e-14 * cSp_[0] / 
		           (
			           1.0 + cSp_[1]*pow(psi, cSp_[2]) * 
								   (1.0 + cSp_[3] * pow(T, cSp_[4])
										             *pow( psi, -cSp_[5]-cSp_[6]*log(T) ) 
									 ) 
								) * nH_ * Zd_ * fSplusgr_;
	kgr_[5] = 1.0e-14 * cSip_[0] / 
		           (
			           1.0 + cSip_[1]*pow(psi, cSip_[2]) * 
								   (1.0 + cSip_[3] * pow(T, cSip_[4])
										             *pow( psi, -cSip_[5]-cSip_[6]*log(T) ) 
									 ) 
								) * nH_ * Zd_ * fSiplusgr_;
	return;
}

void gow17::SetGhostSpec_(const N_Vector y, double yghost[kDimen+n_ghost_]) {
	/*copy the aboundances in y to yghost*/
	for (int i=0; i<kDimen; i++) {
		yghost[i] = NV_Ith_S(y, i);
	}
	/*set the ghost species*/
 	yghost[iC_] = xC_ - yghost[iHCOplus_] -  yghost[iCH_] 
                     - yghost[iCO_] - yghost[iCplus_]; 
	yghost[iO_] = xO_ - yghost[iHCOplus_] -  yghost[iOH_] 
                     - yghost[iCO_] - yghost[iOplus_]; 
	yghost[iHe_] = xHe_ - yghost[iHeplus_]; /*HeI*/
	yghost[iS_] = xS_ - yghost[iSplus_]; /*SI*/
	yghost[iSi_] = xSi_ - yghost[iSiplus_]; /*SiI*/
  yghost[ie_] = yghost[iHeplus_] + yghost[iCplus_] + yghost[iHCOplus_]
                     + yghost[iH3plus_] + yghost[iH2plus_] + yghost[iHplus_]
                     + yghost[iSplus_] + yghost[iSiplus_] + yghost[iOplus_]; /*e*/
	yghost[iH_] = 1.0 - (yghost[iOH_] + yghost[iCH_] + yghost[iHCOplus_]
                     + 3.0*yghost[iH3plus_] + 2.0*yghost[iH2plus_] + yghost[iHplus_]
										 + 2.0*yghost[iH2_]);
	return;
}

/*--------------------------Set physical parameters----------------------*/
void gow17::SetnH(const double nH) {
	nH_ = nH;
	return;
}


void gow17::SetbCO(const double bCO) {
	bCO_ = bCO;
	return;
}

void gow17::SetZd(const double Zd) {
	Zd_ = Zd;
	return;
}

void gow17::SetZg(const double Zg) {
	Zg_ = Zg;
  xC_ = Zg_ * xC_std_;
  xO_ = Zg_ * xO_std_;
  xS_ = Zg_ * xS_std_;
  xSi_ = Zg_ * xSi_std_;
	return;
}

void gow17::SetIonRate(const double ion_rate) {
	ion_rate_ = ion_rate;
	return;
}

void gow17::SetNCO(const double NCO) {
	NCO_ = NCO;
	return;
}

void gow17::SetNH(const double NH) {
	NH_ = NH;
	return;
}

void gow17::SetdxCell(const double dx_cell) {
	if (dx_cell < 0) {
		throw std::runtime_error("SetdxCell: dx_cell < 0\n");
	}
	dx_cell_ = dx_cell;
	return;
}

void gow17::WriteParam(FILE *pf) const {
	fprintf(pf, "ionRate = %0.4e\n", ion_rate_);
	fprintf(pf, "nH = %0.4e\n", nH_);
	fprintf(pf, "xC = %0.4e\n", xC_);
	fprintf(pf, "xO = %0.4e\n", xO_);
	fprintf(pf, "xHe = %0.4e\n", xHe_);
	fprintf(pf, "NCO = %0.4e\n", NCO_);
	fprintf(pf, "dx_cell = %0.4e\n", dx_cell_);
	return;
}

void gow17::WriteRates(FILE *pf) {
	double yprev[kDimen+n_ghost_];
	double rates_cr[n_cr_];
	double rates_2body[n_2body_];
	double rates_ph[n_ph_];
	double rates_gr[n_gr_];

	SetGhostSpec_(y_, yprev);
  ChemInit_(yprev);

	/*cosmic ray reactions*/
	for (int i=0; i<n_cr_; i++) {
		rates_cr[i] = kcr_[i] * yprev[incr_[i]];
	}

	/*2 body reactions*/
	for (int i=0; i<n_2body_; i++) {
		rates_2body[i] =  k2body_[i] * yprev[in2body1_[i]] * yprev[in2body2_[i]];
	}
	/*photo reactions*/
	for (int i=0; i<n_ph_; i++) {
		rates_ph[i] = kph_[i] * yprev[inph_[i]];
  }
	/*grain assisted reactions*/
	for (int i=0; i<n_gr_; i++) {
		rates_gr[i] = kgr_[i] * yprev[ingr_[i]];
	}

	for (int i=0; i<n_cr_; i++) {
		fprintf(pf, "  %12.4e", rates_cr[i]);
	}
	for (int i=0; i<n_2body_; i++) {
		fprintf(pf, "  %12.4e", rates_2body[i]);
	}
	for (int i=0; i<n_ph_; i++) {
		fprintf(pf, "  %12.4e", rates_ph[i]);
	}
	for (int i=0; i<n_gr_; i++) {
		fprintf(pf, "  %12.4e", rates_gr[i]);
	}
	fprintf(pf, "\n");
	return;
}

double gow17::GetTime() const {
	return t_;
}

double gow17::GetnH() const {
	return nH_;
}

double gow17::GetZd() const {
	return Zd_;
}

double gow17::GetIonRate() const {
	return ion_rate_;
}

void gow17::CopyAbd(double *y) const {
	for (int i=0; i<kDimen; i++) {
		y[i] = NV_Ith_S(y_, i);
	}
	return;
}

void gow17::CopyThermoRates(double *y) const {
    y[0] =  LCR_;
    y[1] =  LPE_;
    y[2] =  LH2gr_;
    y[3] = LH2pump_;
    y[4] = LH2diss_;
    y[5] = GCII_;
    y[6] = GCI_;
    y[7] = GOI_;
    y[8] = GLya_;
    y[9] = GCOR_;
    y[10] = GH2_;
    y[11] = GDust_;
    y[12] = GRec_;
    y[13] = GH2diss_;
    y[14] = GHIion_;
    return;
}


double gow17::GetE(const double T, const double xH2, const double xe) {
  return T * Thermo::CvCold(xH2, xHe_, xe);
}

void gow17::SetConstTemp(const double T) {
  const_temp_ = true;
  temp_ = T;
  return;
}

void gow17::IsH2dissHeating(const bool isH2diss_heating) {
  isH2diss_heating_ = isH2diss_heating;
  return;
}

void gow17::IsH2rvCooling(const bool isH2rv_cooling) {
  isH2rv_cooling_ = isH2rv_cooling;
  return;
}

void gow17::IsH2grHeating(const bool isH2gr_heating) {
  isH2gr_heating_ = isH2gr_heating;
  return;
}

double gow17::dEdt_(const double *y, const bool is_store_rates) {
  double T = 0.;
  if (const_temp_) {
    printf("dEdt_: Constant temperature. return 0.\n");
    return 0.;
  } else {
    T = y[iE_] / Thermo::CvCold(y[iH2_], xHe_, y[ie_]);
  }
  double dEdt = 0.;
  /*--------------------------heating-----------------------------*/
  /*cosmic ray ionization of H, He, and H2 */
  /*NOTE: because these depends on rates, make sure ChemInit is called before.*/
  /*NOTE: the kcr_[i] assume the order of equastions are not changed*/
  const double LCR = Thermo::HeatingCr(y[ie_],  nH_,
										        y[iH_],  y[iHe_],  y[iH2_],
		   							        kcr_[icr_H_],  kcr_[icr_He_],  kcr_[icr_H2_]);
  /*photo electric effect on dust*/
  const double LPE = Thermo::HeatingPE((*GPE_), Zd_, T, nH_*y[ie_]);
  /*H2 formation on dust grains*/
  double LH2gr;
  if (isH2gr_heating_) {
    LH2gr = Thermo::HeatingH2gr(y[iH_],  y[iH2_],  nH_,
                                           T,  kgr_[igr_H_]);
  } else {
    LH2gr = 0;
  }
  /*H2 UV pumping*/
  const double dot_xH2_photo = kph_[iph_H2_] * y[iH2_];
  const double LH2pump = Thermo::HeatingH2pump(y[iH_],  y[iH2_],  nH_,
                                               T,  dot_xH2_photo);
  /*H2 photo dissiociation.*/
  double LH2diss;
  if (isH2diss_heating_) {
    LH2diss = Thermo::HeatingH2diss(dot_xH2_photo);
  } else {
    LH2diss = 0;
  }
  /*--------------------------cooling-----------------------------*/
  /* C+ fine structure line */
  const double GCII = Thermo::CoolingCII(y[iCplus_],  nH_*y[iH_],  nH_*y[iH2_],
                                         nH_*y[ie_],  T);
  /* CI fine structure line */
  const double GCI = Thermo:: CoolingCI(y[iC_],  nH_*y[iH_],  nH_*y[iH2_],
                                        nH_*y[ie_],  T);
  /* OI fine structure line */
  const double GOI = Thermo:: CoolingOI(y[iO_],  nH_*y[iH_],  nH_*y[iH2_],
                                        nH_*y[ie_],  T);
  /* collisional exicited lyman alphya line */
  const double GLya = Thermo::CoolingLya(y[iH_], nH_*y[ie_],  T);
  /* CO rotational lines */
  /* Calculate effective CO column density*/
  const double vth = sqrt(2. * Thermo::kb_ * T / mCO_);
  const double nCO = nH_ * y[iCO_];
  double NCOeff, Leff_n, Leff_v, Leff;
  if (isCoolingCOThin_) {
    NCOeff = 0;
  } else if (isNCOeff_global_) {
    const double NCO = NCO_ + y[iCO_]*nH_*dx_cell_;
    if (isbCO_L_) {
      double bCO_L = 1.0e5 * sqrt(NH_/nH_ / 3.086e18);
      NCOeff = NCO / bCO_L;
    } else {
      NCOeff = NCO / bCO_;
    }
  } else if (is_dvdr_) {
    NCOeff = nCO / gradv_;
  } else {
    /*TODO: dx_cell_ can't be too small*/
    if (gradv_ > vth / dx_cell_) {
      NCOeff = nCO / gradv_;
    } else {
      Leff_n = nH_ / gradnH_;
      Leff_v = vth / gradv_;
      Leff = std::min(Leff_n, Leff_v);
      NCOeff = nCO * Leff / vth;
    }
  }
  const double GCOR = Thermo::CoolingCOR(y[iCO_], nH_*y[iH_],  nH_*y[iH2_],
                                         nH_*y[ie_],  T,  NCOeff);
  /* H2 vibration and rotation lines */
  double GH2;
  if (isH2rv_cooling_) {
    GH2 = Thermo::CoolingH2(y[iH2_], nH_*y[iH_],  nH_*y[iH2_],
                                       nH_*y[iHe_],  nH_*y[iHplus_], nH_*y[ie_],
                                       T);
  } else {
    GH2 = 0;
  }
  /* dust thermo emission */
  //const double GDust = Thermo::CoolingDust(Zd_,  nH_, T, (*GISRF_));
  const double GDust = Thermo::CoolingDustTd(Zd_,  nH_, T, 10.);
  /* reconbination of e on PAHs */
  const double GRec = Thermo::CoolingRec(Zd_,  T,  nH_*y[ie_], (*GPE_));
  /* collisional dissociation of H2 */
  const double GH2diss = Thermo::CoolingH2diss(y[iH_],  y[iH2_], k2body_[i2body_H2_H],
                                               k2body_[i2body_H2_H2]);
  /* collisional ionization of HI */
  const double GHIion = Thermo::CoolingHIion(y[iH_],  y[ie_],
                                             k2body_[i2body_H_e]);
  /* Store heating and cooling rates in system */
  if (is_store_rates) {
     LCR_ = LCR;
     LPE_ = LPE;
     LH2gr_ = LH2gr;
     LH2pump_ = LH2pump;
     LH2diss_ = LH2diss;
     GCII_ = GCII;
     GCI_ = GCI;
     GOI_ = GOI;
     GLya_ = GLya;
     GCOR_ = GCOR;
     GH2_ = GH2;
     GDust_ = GDust;
     GRec_ = GRec;
     GH2diss_ = GH2diss;
     GHIion_ = GHIion;
  }
  dEdt = (LCR + LPE + LH2gr + LH2pump + LH2diss)
            - (GCII + GCI + GOI + GLya + GCOR 
                + GH2 + GDust + GRec + GH2diss + GHIion);
  //printf("LCR=%.2e, LPE=%.2e, LH2gr=%.2e, LH2pump=%.2e LH2diss=%.2e\n",
  //       LCR , LPE , LH2gr , LH2pump , LH2diss);
  //printf("GCII=%.2e, GCI=%.2e, GOI=%.2e, GLya=%.2e, GCOR=%.2e\n",
  //       GCII , GCI , GOI , GLya , GCOR);
  //printf("GH2=%.2e, GDust=%.2e, GRec=%.2e, GH2diss=%.2e, GHIio=%.2e\n",
  //       GH2 , GDust , GRec , GH2diss , GHIion);
  //printf("T=%.2e, dEdt=%.2e\n", T, dEdt);
  return dEdt;
}

void gow17::SetGradv(const double gradv) {
	gradv_ = gradv;
	return;
}

void gow17::SetNCOeffGlobal(const bool isNCOeff_global) {
  isNCOeff_global_ = isNCOeff_global;
  return;
}

void gow17::Setdvdr(const bool is_dvdr) {
  is_dvdr_ = is_dvdr;
  return;
}

void gow17::SetbCOL(const bool isbCO_L) {
  isbCO_L_ = isbCO_L;
  return;
}

void gow17::SetCoolingCOThin(const bool isCoolingCOThin) {
  isCoolingCOThin_ = isCoolingCOThin;
  return;
}

double gow17::CII_rec_rate_(const double temp) {
  double A, B, T0, T1, C, T2, BN, term1, term2, alpharr, alphadr;
  A = 2.995e-9;
  B = 0.7849;
  T0 =  6.670e-3;
  T1 = 1.943e6;
  C = 0.1597;
  T2 = 4.955e4;
  BN = B + C * exp(-T2/temp);
  term1 = sqrt(temp/T0);
  term2 = sqrt(temp/T1);
  alpharr = A / ( term1*pow(1.0+term1, 1.0-BN) * pow(1.0+term2, 1.0+BN) );
  alphadr = pow( temp, -3.0/2.0 ) * ( 6.346e-9 * exp(-1.217e1/temp) +
        9.793e-09 * exp(-7.38e1/temp) + 1.634e-06 * exp(-1.523e+04/temp) );
  return (alpharr+alphadr);
}

double gow17::GetxCtot() const {
  return xC_;
}

double gow17::GetxOtot() const {
  return xO_;
}

double gow17::GetxHetot() const {
  return xHe_;
}

void gow17::SetxCtot(const double xC) {
  xC_ = xC;
  return;
}

void gow17::SetxOtot(const double xO) {
  xO_ = xO;
  return;
}

void gow17::SetRadField(double *GPE, double *Gph, double *GISRF){
  GPE_ = GPE;
  GISRF_ = GISRF;
  Gph_ = Gph;
  return;
}
void gow17::SetfH2gr(double fH2gr) {
  fH2gr_ = fH2gr;
  return;
}

void gow17::SetfHplusgr(double fHplusgr) {
  fHplusgr_ = fHplusgr;
  return;
}

void gow17::SetfCplusCR(double fCplusCR) {
  fCplusCR_ = fCplusCR;
  return;
}

void gow17::SetfCplusgr(double fCplusgr) {
  fCplusgr_ = fCplusgr;
  return;
}

void gow17::SetfHeplusgr(double fHeplusgr) {
  fHeplusgr_ = fHeplusgr;
  return;
}

void gow17::SetfSplusgr(double fSplusgr) {
  fSplusgr_ = fSplusgr;
  return;
}

void gow17::SetfSiplusgr(double fSiplusgr) {
  fSiplusgr_ = fSiplusgr;
  return;
}
