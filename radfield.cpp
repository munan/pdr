/* slab radiation, beamed.
 * Date: Oct 2, 2015, Author: Munan Gong
 */

#include "radfield.h"

RadField::RadField(const long int ngrid, const double G0, const double Zd)
  :ngrid_(ngrid),
   G0_(G0),
   Zd_(Zd)
{
	bH2_ = 3.0e5;
	isfsH2_ = true;
	isfsCO_ = true;
	isfsC_ = false;
  fs_CO_ = 1.0;
  fs_H2_ = 1.0;
  fs_C_ = 1.0;
  isdust_= true;
  GPE = new double [ngrid];
  GISRF = new double [ngrid];
  Gph = new double* [ngrid];
  for (int i=0; i<ngrid; i++) {
    Gph[i] = new double [gow17::n_ph_];
  }
}

RadField::~RadField() {
  for (int i=0; i<ngrid_; i++) {
    delete [] Gph[i];
  }
  delete [] GPE;
  delete [] GISRF;
  delete [] Gph;
}

void RadField::Beamed(const long int igrid, const double NH,
                      const double NH2, const double NCO, const double NC) {
  double AV;
  if (isdust_) {
    AV = NH * Zd_ / 1.87e21;
    GPE[igrid] = (G0_/2.) * exp(-NH * Thermo::sigmaPE_ * Zd_);
    GISRF[igrid] = (G0_/2.) * exp(-NH * Thermo::sigmaISRF_ * Zd_);
  } else {
    AV = 0.;
    GPE[igrid] = (G0_/2.);
    GISRF[igrid] = (G0_/2.);
  }
  /*photo-reactions*/
	for (int i=0; i<gow17::n_ph_; i++) {
		Gph[igrid][i] = (G0_/2.) * exp( -gow17::kph_avfac_[i] * AV );
	}
	/*(2) h nu + CO -> *C + *O            --self-shielding and shielding by H2
	  (5) h nu + H2 -> *H + *H            --self- and dust shielding*/
  /*(0) h nu + *C -> C+ + *e*/
 /*self-sheild and sheildig by H2 of CO*/
  if (isfsCO_) {
    fs_CO_ = Shielding::fShield_CO_V09(NCO, NH2);
    Gph[igrid][gow17::iph_CO_] *= fs_CO_;
  }
  /*self-shielding of CI*/
  if (isfsC_) {
    fs_C_ = Shielding::fShield_C(NC, NH2);
    Gph[igrid][gow17::iph_C_] *= fs_C_;
  }
  /*self-shielding of H2*/
  if (isfsH2_) {
    fs_H2_ = Shielding::fShield_H2(NH2, bH2_);
    Gph[igrid][gow17::iph_H2_] *= fs_H2_;
  }

  /*debug*/
//  printf("igrid=%ld, GPH=%.1e\n", igrid, GPE[igrid]);
//  printf("Gph=");
//  for (long int i=0; i<gow17::n_ph_; i++) {
//    printf("%.1e     ", Gph[igrid][i]);
//  }
//  printf("\n");
//
  return;
}

void RadField::IsotropicApprox(const long int igrid, const double NH,
                         const double NH2, const double NCO, const double NC) {
  Beamed(igrid, NH*2, NH2*2, NCO*2, NC*2);
  return;
}

void RadField::Isotropic(const long int igrid, const double NH,
                         const double NH2, const double NCO, const double NC,
                         const long int nangle) {
  double sumGPE = 0.;
  double sumGISRF = 0.;
  double sumGph[gow17::n_ph_];
  double dist_fac = 0.;
  for (int i=0; i<gow17::n_ph_; i++) {
    sumGph[i] = 0.;
  }
  /*at 0 and pi/2*/
  Beamed(igrid, NH, NH2, NCO, NC);
  sumGPE += 0.5 * (2*GPE[igrid]);
  sumGISRF += 0.5 * (2*GISRF[igrid]);
  for (int i=0; i<gow17::n_ph_; i++) {
    sumGph[i] += 0.5 * (2*Gph[igrid][i]);
  }
  /*if less than 1 nangle, then average theta = 0, pi/2*/
  if (nangle <= 1) {
    GPE[igrid] = sumGPE * 0.5;
    GISRF[igrid] = sumGISRF * 0.5;
    for (int i=0; i<gow17::n_ph_; i++) {
      Gph[igrid][i] = sumGph[i] * 0.5;
    }
    return;
  }
  /*else, do the values in between*/
  const double dx = 1. / nangle;
  double x;
  for (long int i=1; i<nangle; i++){
    x = dx * i;
    dist_fac = 1. / x;
    Beamed(igrid, NH*dist_fac, NH2*dist_fac, NCO*dist_fac, NC*dist_fac);
    sumGPE += 2*GPE[igrid];
    sumGISRF += 2*GISRF[igrid];
    for (int i=0; i<gow17::n_ph_; i++) {
      sumGph[i] +=  2*Gph[igrid][i];
    }
  }
  /*assign results*/
  GPE[igrid] = sumGPE * dx / 2;
  GISRF[igrid] = sumGISRF * dx / 2;
  for (int i=0; i<gow17::n_ph_; i++) {
    Gph[igrid][i] = sumGph[i] * dx / 2;
  }
  return;
}



void RadField::IsMolSheildH2(const bool isfsH2) {
	isfsH2_ = isfsH2;
	return;
}

void RadField::IsMolSheildCO(const bool isfsCO) {
	isfsCO_ = isfsCO;
	return;
}

void RadField::IsSelfSheildC(const bool isfsC) {
	isfsC_ = isfsC;
	return;
}

void RadField::SetbH2(const double bH2) {
	bH2_ = bH2;
	return;
}

void  RadField::IsdustSheild(const bool isdust) {
  isdust_ = isdust;
  return;
}

double RadField::GetfShieldH2mol() {
  return fs_H2_;
}

double RadField::GetfShieldCOmol() {
  return fs_CO_;
}
