#include "coolingFunction.h"
#include "thermo.h"
#include <stdio.h>
#include <math.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

CoolingFunction::CoolingFunction(Slab &myslab)
  :myslab_(myslab),
	 nE_(myslab.ode_.GetnE())
{
  //input parameters
  kcr_ = myslab_.ode_.GetIonRate();
  nH_ = myslab_.ode_.GetnH();
  Zd_ = myslab_.ode_.GetZd();
  xCtot_ = myslab_.ode_.GetxCtot();
  xOtot_ = myslab_.ode_.GetxOtot();
  xHetot_ = myslab_.ode_.GetxHetot();
  ngrid_ = myslab_.ngrid_;
  //printf("kcr=%.2e, nH=%.2e, Zd=%.2e, xCtot=%.2e, xOtot=%.2e, xHetot=%.2e\n",
  //        kcr_, nH_, Zd_, xCtot_, xOtot_, xHetot_);
	T_ = new double [ngrid_];
	GPE_ = new double [ngrid_];
	GCI_ = new double [ngrid_];
	GCO_ = new double [ngrid_];
	GH2_ = new double [ngrid_];

  //output paramters
	fe_ = new double [ngrid_];
  fHplus_ = new double [ngrid_];
  fH2_ = new double [ngrid_];
  fHI_ = new double [ngrid_];
  fCplus_ = new double [ngrid_];
  fCO_ = new double [ngrid_];
  fCI_ = new double [ngrid_];
  fOI_ = new double [ngrid_];
	yE_ = new double* [ngrid_];
	for (int i=0; i<ngrid_; i++) {
		yE_[i] = new double [nE_];
	}
  
  //set values for parameters
	for (int i=0; i<ngrid_; i++) {
    //read temperature
    double xe_true = myslab_.y_[i][myslab_.ode_.id("H+")]+
      myslab_.y_[i][myslab_.ode_.id("C+")] + myslab_.y_[i][myslab_.ode_.id("S+")]
      + myslab_.y_[i][myslab_.ode_.id("Si+")] + myslab_.y_[i][myslab_.ode_.id("HCO+")]
      + myslab_.y_[i][myslab_.ode_.id("H3+")] + myslab_.y_[i][myslab_.ode_.id("H2+")]
      + myslab_.y_[i][myslab_.ode_.id("O+")] + myslab_.y_[i][myslab_.ode_.id("He+")];
    T_[i] = myslab_.y_[i][myslab_.ode_.id("E")] / Thermo::CvCold(
        myslab_.y_[i][myslab_.ode_.id("H2")], xHetot_, xe_true);
    //read radiation fields
    GPE_[i] = myslab_.prad_->GPE[i];
    GCI_[i] = myslab_.prad_->Gph[i][NL99p::iph_C_];
    GCO_[i] = myslab_.prad_->Gph[i][NL99p::iph_CO_];
    GH2_[i] = myslab_.prad_->Gph[i][NL99p::iph_H2_];
    //initial value for output parameters are zero.
    fe_[i] = 0;
    fHplus_[i] = 0;
    fH2_[i] = 0; 
    fHI_[i] = 0;
    fCplus_[i] = 0;
    fCO_[i] = 0;
    fCI_[i] = 0;
    fOI_[i] = 0;
    /*initialize yE_*/
    for (int j=0; j<nE_; j++) {
      yE_[i][j] = 0.;
    }
	}
  return;
}

CoolingFunction::~CoolingFunction(){
  delete [] fe_;
  delete [] fHplus_;
  delete [] fH2_;
  delete [] fHI_;
  delete [] fCplus_;
  delete [] fCO_;
  delete [] fCI_;
  delete [] fOI_;
	for (int i=0; i<ngrid_; i++) {
		delete [] yE_[i];
	}
	delete [] yE_;
  return;
}

void CoolingFunction::ComputeAbundances() {
  get_fH2_CR_();
  get_fe_();
  get_fCplus_();
  get_fHplus_();
  get_fCO_();
	for (int i=0; i<ngrid_; i++) {
    fHI_[i] = MAX(1. - 2.0*fH2_[i] - fHplus_[i], 0.);
    fCI_[i] = MAX(xCtot_ - fCO_[i] - fCplus_[i], 0.);
    fOI_[i] = MAX(xOtot_ - fCO_[i], 0.);
  }
  return;
}

void CoolingFunction::WriteAbundances(FILE *pf) {
	for (int i=0; i<ngrid_; i++) {
	  fprintf(pf, "%12.4e  ", fe_[i]);
	  fprintf(pf, "%12.4e  ", fHplus_[i]);
	  fprintf(pf, "%12.4e  ", fH2_[i]);
	  fprintf(pf, "%12.4e  ", fHI_[i]);
	  fprintf(pf, "%12.4e  ", fCplus_[i]);
	  fprintf(pf, "%12.4e  ", fCO_[i]);
	  fprintf(pf, "%12.4e  ", fCI_[i]);
	  fprintf(pf, "%12.4e  ", fOI_[i]);
		fprintf(pf, "\n");
	}
  return;
}

void CoolingFunction::get_fH2_CR_() {
  double kgr = 3.0e-17*Zd_;
  double R = kgr * nH_ / (2. * kcr_);
  double b = - (1.65*1.5 + 2.*R);
  double a = 1.65*0.7;
  double c = R;
  double x_H2 = (-b - sqrt(b*b - 4*a*c) )/(2.*a);
	for (int i=0; i<ngrid_; i++) {
    fH2_[i] = x_H2;
  }
  return;
}

void CoolingFunction::get_fH2_() {
  double kgr = 3.0e-17*Zd_;
  double R = kgr * nH_ / (2. * kcr_);
  double a = 1.65*0.7;
  double c = R;
  double k_FUV, b;
	for (int i=0; i<ngrid_; i++) {
    k_FUV = 5.7e-11 * GH2_[i];
    b = - (1.65*1.5 + 2.*R + k_FUV/(2.*kcr_));
    fH2_[i] = (-b - sqrt(b*b - 4*a*c) )/(2.*a);
  }
  return;
}

void CoolingFunction::get_fCO_() {
  double kcr16 = kcr_/1.0e-16;
  double term1 = 4.0e3*Zd_/(kcr16*kcr16);
  double ncrit2 = 0;
  double x_CO = 0;
  double x_CO_max1 = 0;
  double x_CO_max2 = 0;
  term1 = MAX(term1, 1.);
	for (int i=0; i<ngrid_; i++) {
    ncrit2 = 2.*pow(term1, pow(GCO_[i], 1./3.)) * (50.*kcr16/pow(Zd_, 1.4));
    if (nH_ >= ncrit2) {
      x_CO = 1.;
    } else {
      x_CO = nH_/ncrit2;
    }
    x_CO *= xCtot_;
    x_CO_max1 = xCtot_ - fCplus_[i];
    x_CO_max2 = xCtot_*fH2_[i]*2.;
    x_CO = MIN(x_CO, x_CO_max1);
    x_CO = MIN(x_CO, x_CO_max2);
    fCO_[i] = x_CO;
  }
  return;
}

double CoolingFunction::fHplus_e_(double x_e, double x_Cplus, double x_H2,
                                  double temp, double G_PE) {
  const double small_ = 1e-50;
  double x_H = 1. - 2. * x_H2 - (x_e - x_Cplus);
  x_H = MAX(x_H, 0.0);
  double k_Hplus_e = 2.753e-14 * pow( 315614.0 / temp, 1.5) * pow( 
               1.0 + pow( 115188.0 / temp, 0.407) , -2.242 );
  const double cHp_[7] = {12.25, 8.074e-6, 1.378, 5.087e2,
                               1.586e-2, 0.4723, 1.102e-5};
  double psi_gr = 1.7 * G_PE * sqrt(temp)/(nH_ * x_e + small_);
  double k_Hplus_gr = 1.0e-14 * cHp_[0] / 
             (
               1.0 + cHp_[1]*pow(psi_gr, cHp_[2]) * 
                 (1.0 + cHp_[3] * pow(temp, cHp_[4])
                               *pow( psi_gr, -cHp_[5]-cHp_[6]*log(temp) ) 
                 ) 
              ) * Zd_;
  double k_cr_H = kcr_ * (2.3*x_H2 + 1.5*x_H);
  double c = k_cr_H * x_H / nH_;
  double x_Hplus = c/(k_Hplus_e *  x_e + k_Hplus_gr);
  x_Hplus = MIN(x_Hplus, 1.0);
  return x_Hplus;
}

double CoolingFunction::CII_rec_rate_(const double temp) {
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

double CoolingFunction::fCplus_e_(double x_e, double x_H2, 
                                  double temp, double G_PE, double G_CI) {
  const double small_ = 1e-50;
  double k_C_cr = 3.85 * kcr_;
  double k_C_photo = 3.5e-10*G_CI;
  double k_Cplus_e = CII_rec_rate_(temp);
  double psi_gr = 1.7 * G_PE * sqrt(temp)/(nH_ * x_e + small_);
  const double cCp_[7] = {45.58, 6.089e-3, 1.128, 4.331e2, 4.845e-2,
                          0.8120, 1.333e-4};
  double k_Cplus_gr = 1.0e-14 * cCp_[0] / 
		           (
			           1.0 + cCp_[1]*pow(psi_gr, cCp_[2]) * 
								   (1.0 + cCp_[3] * pow(temp, cCp_[4])
										             *pow( psi_gr, -cCp_[5]-cCp_[6]*log(temp) ) 
									 ) 
								) * Zd_;
  double k_Cplus_H2 = 3.3e-13 * pow(temp, -1.3) * exp(-23./temp);
  double c = (k_C_cr + k_C_photo) / nH_;
  double al = k_Cplus_e*x_e + k_Cplus_gr + k_Cplus_H2*x_H2 + c;
  double ar = xCtot_ * c;
  double x_Cplus = ar / al;
  return x_Cplus;
}

double CoolingFunction::fe_e_(double x_e, double x_H2, double temp,
                              double G_PE, double G_CI) {
  double xCplus = fCplus_e_(x_e, x_H2, temp, G_PE, G_CI);
  double xHplus = fHplus_e_(x_e, xCplus, x_H2, temp, G_PE);
  return (xCplus + xHplus);
}

void CoolingFunction::get_fe_() {
  const double rtol = 1e-2;
  const double small_x = 1e-6;
  const double small_f = 1e-20;
  const int maxiter = 20;
  double x = 0;
  double f = 0;
  double xnew = 0;
  double fnew = 0;
  int niter = 0;
  double xprev = 0.5;
  double fprev = 0;
  double a = 0.0;
  double fa = 0.0;
  double b = 1.0;
  bool flag = true;
	for (int i=0; i<ngrid_; i++) {
    niter = 0;
    fprev = fe_e_(xprev, fH2_[i], T_[i], GPE_[i], GCI_[i]) - xprev;
    while (1) {
      x = fe_e_(xprev, fH2_[i], T_[i], GPE_[i], GCI_[i]);
      if ( abs((x-xprev)/(xprev+small_x)) < rtol || abs(fprev) < small_x || 
           abs((a-b)/(b+small_x)) < rtol) {
        break;
      }
      if (niter > maxiter) {
        printf("get_fe_: WARNING: niter>maxiter(=%d), x=%.2e, xprev=%.2e\n", 
               maxiter, x, xprev);
        printf("a=%.2e, b=%.2e\n", a, b);
        //TODO: throw run time error?
        break;
      }
      f = fe_e_(x, fH2_[i], T_[i], GPE_[i], GCI_[i]) - x;
      if ( abs(f - fprev) < small_f ) {
        xnew = (x + xprev)/2.;
      } else {
        xnew = x - f * (x-xprev)/(f-fprev);
      }
      if (xnew < a || xnew > b) {
        xnew = (a + b)/2.;
      }
      if (flag) {
        fa = fe_e_(a, fH2_[i], T_[i], GPE_[i], GCI_[i])-a;
      }
      fnew = fe_e_(xnew, fH2_[i], T_[i], GPE_[i], GCI_[i])-xnew;
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
    fe_[i] = xprev;
  }
  return;
}

void CoolingFunction::get_fHplus_() {
	for (int i=0; i<ngrid_; i++) {
    fHplus_[i] = fHplus_e_(fe_[i], fCplus_[i], fH2_[i], T_[i], GPE_[i]);
  }
  return;
}

void CoolingFunction::get_fCplus_() {
	for (int i=0; i<ngrid_; i++) {
    fCplus_[i] = fCplus_e_(fe_[i], fH2_[i], T_[i], GPE_[i], GCI_[i]);
  }
  return;
}

void CoolingFunction::ComputeThermoRates() {
  const double mCO = 4.68e-23;
  double vth, nCO, NCOeff, Leff_n, Leff_v, Leff, kcr_H_fac;
  double GLya, GCII, GCI, GOI, GCOR, GRec;
  double LCR, LPE;
	for (int i=0; i<ngrid_; i++) {
    //cooling
    GLya = Thermo::CoolingLya(fHI_[i], nH_*fe_[i],  T_[i]);
    GCII = Thermo::CoolingCII(fCplus_[i],  nH_*fHI_[i],  nH_*fH2_[i], 
                              nH_*fe_[i],  T_[i]);
    GCI = Thermo:: CoolingCI(fCI_[i],  nH_*fHI_[i],  nH_*fH2_[i],
                             nH_*fe_[i],  T_[i]);
    GOI = Thermo:: CoolingOI(fOI_[i],  nH_*fHI_[i],  nH_*fH2_[i],
                             nH_*fe_[i],  T_[i]);
    //calculate effective column for CO
    vth = sqrt(2. * Thermo::kb_ * T_[i] / mCO);
    nCO = nH_ * fCO_[i];
    if (myslab_.ode_.gradv_ > vth / myslab_.ode_.dx_cell_) {
      NCOeff = nCO / myslab_.ode_.gradv_;
    } else {
      Leff_n = nH_ / myslab_.ode_.gradnH_;
      Leff_v = vth / myslab_.ode_.gradv_;
      Leff = std::min(Leff_n, Leff_v);
      NCOeff = nCO * Leff / vth;
    }
    GCOR = Thermo::CoolingCOR(fCO_[i], nH_*fHI_[i],  nH_*fH2_[i],
                              nH_*fe_[i],  T_[i],  NCOeff);
    GRec = Thermo::CoolingRec(Zd_,  T_[i],  nH_*fe_[i], GPE_[i]);
    //heating
    kcr_H_fac = 1.15 * 2*fH2_[i] + 1.5 * fHI_[i];
    LCR = Thermo::HeatingCr(fe_[i],  nH_,
										        fHI_[i],  xHetot_,  fH2_[i],
		   							        kcr_*kcr_H_fac,  kcr_*1.1, kcr_*2.0*kcr_H_fac);
    LPE = Thermo::HeatingPE(GPE_[i], Zd_, T_[i], nH_*fe_[i]);
    yE_[i][0] =  LCR;
    yE_[i][1] =  LPE;
    yE_[i][2] =  0.;
    yE_[i][3] = 0.;
    yE_[i][4] = 0.;
    yE_[i][5] = GCII;
    yE_[i][6] = GCI;
    yE_[i][7] = GOI;
    yE_[i][8] = GLya;
    yE_[i][9] = GCOR;
    yE_[i][10] = 0.;
    yE_[i][11] = 0.;
    yE_[i][12] = GRec;
    yE_[i][13] = 0.;
    yE_[i][14] = 0.;
  }
  return;
}

void CoolingFunction::WriteThermoRates(FILE *pf) {
	for (int i=0; i<ngrid_; i++) {
		for (int j=0; j<nE_; j++) {
			fprintf(pf, "%12.4e  ", yE_[i][j]);
		}
		fprintf(pf, "\n");
	}
  return;
}
